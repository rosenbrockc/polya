#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
AUTHORS: Conrad W. Rosenbrock, Wiley S. Morgan (June 2015)

Classes to support the calculation of coefficients for specific terms in a product
of multinomials. Construct a product class by specifying the exponent and target term
and then add multinomials using the Product instance's append(). The coefficient is
then available from the coeff().
"""
from functools import reduce
class Sequence(object):
    """Represents an exponent-limited sequence with a single root. Here, sequence represents a
    sequence of integer values k_1, k_2...k_j that are the exponents of a single term in a multinomial.
    The root of the sequence is one of the k_i; its children become sets of sequences including
    variables to the right of i.
    """
    def __init__(self, root, possibles, i, powersum, targets, parent=None):
        """Initializes a sequence collector for a variable. 'Term' refers to a product
        of variables like x^i.y^j.z^k, then x has index 0, y has 1, etc.

        :arg root: the exponent of the variable to the left in the multinomial term.
        :arg possibles: a list of possible values for each variable in the multinomial.
        :arg i: the index of the variable being sequenced in the term.
        :arg powersum: the maximum value that the sum of exponents in the sequence is allowed to have.
        :arg parent: a Sequence instance for the variable the *left* of this one
          (i.e. has index i-1).
        """
        self._root = root
        self.used = root + (0 if parent is None else parent.used)
        self.parent = parent

        #We only keep recursively defining sequences until we run out of variables in
        #the term. Possibles is a list of possible exponents for each variable in the
        #term and has the same number of items as variables in the term.
        if i < len(targets):
            #Filter the possible values for the variable being considered based on the
            #exponent of the multinomial. When multinomials are expanded, the sum of
            #the exponents in any term must be less than the exponent on the multinomial
            #times the maximum power of any of its (unexpanded) terms.

            #We find all the possible values for this variable by ensuring that:
            # 1) it's exponent is compatible with the exponents of all variables to the left of it.
            # 2) the exponent we are suggesting is in the list of possible values for the variable.
            # 3) the exponent remains positive.
            self.kids = [Sequence(p-root, possibles, i+1, powersum, targets, self) 
                         for p in possibles if p-root >= 0
                         and p-root <= targets[i]
                         and abs(p - root) <= powersum-self.used 
                         and abs(p-self.used) % possibles[1] == 0]
        else:
            self.kids = []

        self._kidcount = None
        self.varcount = len(targets)

    @property
    def kidcount(self): #pragma: no cover
        """Returns the number of children and grandchildren to the last generation."""
        if self._kidcount is None:
            _kidcount = sum([k.kidcount for k in self.kids])
            if _kidcount == 0:
                _kidcount = len(self.kids)
        return _kidcount

    def expand(self, depth=0):
        """Recursively generates a list of all relevant sequences for this multinomial term."""
        #Iterate through the child sequences and add their variable root values if
        #the total sequence sums to the target.
        sequences = []
        for kid in self.kids:
            for seq in kid.expand(depth+1):
                #Here is where the recursion happens; we add the sequence of this variable's
                #children to the right of this root.
                sequences.append((self._root,) + seq)

        if len(self.kids) == 0:
            if depth == self.varcount-1:
                return [(self._root,)]
            else: #pragma: no cover
                return [(self._root,) + (0,)*(self.varcount-(depth+1))]
        else:
            return sequences

    def expand_noappend(self, sequences, start, varindex): #pragma: no cover
        """Implements an expansion that doesn't use python's append."""
        if len(sequences) == 0:
            if self.kidcount == 0:
                raise ValueError("This can't happen!")
            sequences = [[0]*self.varcount for i in range(self.kidcount)]

        cursor = start
        for kid in self.kids:
            kid.expand_noappend(sequences, cursor, varindex+1)
            cursor += kid.kidcount

            if varindex == self.varcount-1:
                cursor += 1

            for k in range(start, cursor):
                sequences[k][varindex-1] = self._root
            start = cursor

        if self.kidcount == 0:
            sequences[cursor][varindex-1] = self._root
            
        return sequences

class Product(object):
    """Represents a product of multinomials for which only a single term is interesting."""
    def __init__(self, coefficient, targets):
        """Initializes the empty product of multinomials.

        :arg coefficient: the scalar integer multiplying this product of multinomials.
        :arg targets: a list of exponents for the only interesting term in the product. The
          list is in the order that the variables appear in each multinomial.
        """
        self.coefficient = coefficient
        self.targets = targets
        self.multinoms = []

    def coeff(self):
        """Returns the coefficient of the term with the target exponents if all the multinomials
        in the product were expanded and had their terms collected.
        """
        #If this is an isolated multinomial, we only need to check the coefficient of the target
        #term.
        if len(self.multinoms) == 1:
            if all([self.multinoms[0].power-t>0 for t in self.targets]):
                return 0
            else:
                return self.multinoms[0].nchoosekm(self.targets)*self.coefficient
        
        from itertools import product
        #Get a list of the possible exponents for each variable in each of the multinomials.
        #We start with the first variable and choose only those combinations of exponents
        #across *all* the multinomials that give the correct target exponent for that variable.
        possibles = [n.possible_powers for n in self.multinoms]
        
        seq0 = [s for s in product(*possibles) if sum(s) == self.targets[0]]
        #Next, we construct Sequence instances for each of the first variable compatible
        #possibilities and follow them through to the other variables.
        coeffs = 0
        for seq in seq0:
            mnseq = []
            #Each sequence calculated from the first variable has an entry for each multinomial
            #in this product. The Sequence instances construct smart sequences for the remaining
            #variables in each multinomial separately
            for i in range(len(seq)):
                varseq = Sequence(seq[i], possibles[i], 1, self.multinoms[i].powersum, self.targets)
                mnseq.append(varseq.expand())
            coeffs += self._sum_sequences(mnseq)

        return int(coeffs)*self.coefficient

    def _sum_sequences(self, mnseq):
        """Sums all the possible combinations of relevant sequences based of the variable sequence
        lists specified.

        :arg mnseq: a list of possible variable sequences in each multinomial (one for each multinomial)
          that might contribute to the correct target variable.
        """
        from itertools import product
        from operator import mul
        from functools import reduce
        #We can also filter the sequences by enforcing the constraints that the exponents
        #correctly reproduce the target across all the multionmials in the product. Get hold
        #of all the combinations of sequences across the multinomials and check each for
        #conformance to the targets.
        coeffs = 0
        for seq in product(*mnseq):
            expsum = [sum(zs) for zs in zip(*seq)]
            if expsum == self.targets:
                coeffs += reduce(mul, [m.nchoosekm(s) for m, s in zip(self.multinoms, seq)])

        return coeffs

    def __str__(self):
        #First we need to sort the multinomials by their exponent.
        sortedmns = sorted(self.multinoms, key=(lambda m: (m.exponent,m.power)), reverse=True)
        return str(self.coefficient) + ''.join([str(mn) for mn in sortedmns])            

class Multinomial(object):
    """Represents a multinomial expansion."""
    def __init__(self, power, coeff, arrowings, exponent=1):
        """Sets up the multinomial.

        :arg powers: the power on each of the *unexpanded* variables in the multinomial;
          of the form (x^2+y^2+z^2) => 2.
        :arg coeff: the coefficient on each of the variables in the multinomial; e.g.
          (x^2 + a y^2) => a. For more than two variables, each color that *has* arrowing
          gets the coefficient while the others get 1.
        :arg exponent: the exponent of the entire multinomial.
        """
        self.power = power
        self.coeff = coeff
        self.arrowings = arrowings
        self.exponent = exponent
        self.powersum = power*exponent
        """Returns the integer value that all term exponents in the multinomial should
        sum to (or be less than)."""
        self.possible_powers = list(range(0,power*exponent+1, power))
        """For each variable being considered, determines the possible powers based
        on the exponent in the multinomial."""

    def __str__(self):
        #We want to print the multinomial out in a nice, readable way, similar to how
        #they are presented in Mathematica.
        contents = ' + '.join(["{}^{}".format(self.coeff if a else 1, self.power) for a in self.arrowings])
        return "({})^{}".format(contents, self.exponent)

    def normed_seq(self, seq):
        """Normalizes the specified sequence using the powers of unexpanded terms in the multinomial.
        
        :arg seq: a list of exponents in an *expanded* term.
        """
        return [int(ai/self.power) for ai in seq]

    def nchoosekm(self, sequence):
        """Returns the number of different ways to partition an n-element
        set into disjoint subsets of sizes k1, ..., km.

        :arg sequence: an un-normed tuple of form (k1, k2, k3).
        """
        prod = 1
        if not all([seq%self.power == 0 for seq in sequence]):
            return 0
        else:
            from operator import mul
            normseq = self.normed_seq(sequence)
            for i in range(len(sequence)):
                nsum = sum(normseq[0:i+1])
                prod *= Multinomial.nchoosek(nsum, normseq[i])

            #Add the contribution from the coefficients of the variable *inside*
            #the multionomial.
            pcoeff = 1
            for iseq, arrow in enumerate(self.arrowings):
                if arrow:
                    pcoeff *= self.coeff**normseq[iseq]

            return prod*pcoeff
        
    @staticmethod
    def nchoosek(n, k):
        """This implementation was taken from "Binomial Coefﬁcient Computation: Recursion 
        or Iteration?" by Yannis Manolopoulos, ACM SIGCSE Bulletin InRoads, Vol.34, No.4, 
        December 2002. http://delab.csd.auth.gr/papers/SBI02m.pdf It is supposed to be robust 
        against large, intermediate values and to have optimal complexity.
        """
        if k < 0 or k > n: #pragma: no cover
            return 0
        if k==0 and n == 0:
            return 1
        t = 1
        if k < n-k:
            for i in range(n, n-k, -1):
                t = t*i/(n-i+1)
        else:
            for i in range(n, k, -1):
                t = t*i/(n-i+1)

        return t

def group(gen): #pragma: no cover
    """Generates an entire group using the specified generators by applying generators
    to each other recursively.

    :arg gen: a list of generators as integers.
    """
    from operator import itemgetter as iget
    def g_apply(operations, source, groupi=None):
        """Applies the specified group operations to the source list of elements and then
        appends it to the group if it is unique.
        """
        result = list(iget(*operations)(source))
        if groupi is not None and result not in groupi:
            groupi.append(result)
        return result

    #Make sure the group is zero-based for python.
    if not 0 in gen[0]:
        ngens = [list([e-1 for e in g]) for g in gen]
    else:
        ngens = gen

    groupi = []
    for i in ngens:
        for j in ngens: #filter(lambda k: k!=i, ngens):
            c = g_apply(i, j, groupi)
            d = g_apply(i, c, groupi)
            while d != c:
                d = g_apply(i, d, groupi)

    while True:
        group2 = []
        for i in ngens:
            for h in groupi:
                d = g_apply(i, h)
                if d not in groupi and d not in group2:
                    group2.append(d)

        groupi.extend(group2)
        if len(group2) == 0:
            break
    return(groupi)

def _group_to_cyclic(group, limit=None):
    """Determines the degeneracy of each r-cycle in the specified group operations."""
    result = []
    #We allow filtering so that the unit testing can access the cyclic form of the group.
    if 0 not in group[0][0]: #pragma: no cover
        group = [[[j - 1 for j in i] for i in t] for t in group]
    
    if limit is not None: #pragma: no cover
        filtered = group[limit[0]:limit[1]]
    else:
        filtered = group

    for operation in filtered:
        #visitedp has the same # of elements as the site group operation and
        #is used to make sure each element in the array is visited as
        #we loop through in a *non-sequential* order.
        #visitedc has the same number of elements as the arrow group
        #opertaion and is used to make sure each element of that array
        #is visited as we loop through it in a "non-sequential" order.
        visitedp = [0]*len(operation[0])
        visitedc = [0]*len(operation[1])
        polynomials = {}
        
        arrow_cycle_len = []
        #We start by finding the cylce lengths of the arrow group
        #operations so that we can later see how many of them are
        #divosors of the site group operations.
        while 0 in visitedc:
            #Start with the first elemenet in the arrow group that
            #hasn't been visited yet. All cycles have length > 0.
            cursor = vindex = visitedc.index(0)
            cyc_len = 1
            visitedc[cursor] = 1
            cursor = operation[1][cursor]
            while cursor != vindex:
                visitedc[cursor] = 1
                cyc_len += 1
                cursor = operation[1][cursor]
            arrow_cycle_len.append(cyc_len)

        while 0 in visitedp:
            #Start with the first element in the site group that hasn't
            #been visited yet. The first non-trivial polynomials have
            #powers > 0.
            cursor = vindex = visitedp.index(0)
            powers = 1 
            #change the current position to having been visited; move the cursor.
            visitedp[cursor] = 1
            cursor = operation[0][cursor]
            #The power of the variables in the polynomials is equal to
            #the number of group operations separating the cursor's
            #current position from its *value* in the group operations
            #list.
            while cursor != vindex:
                visitedp[cursor] = 1
                powers += 1
                cursor = operation[0][cursor]
            #We now have everything need to construct part of the
            #polynomial. This is done by taking powers and using it to
            #construct an array of length equal to the number of
            #elements in the system each entry in the array is set to
            #be equal to powers.

            #To get the coefficients right we need to see how many of
            #the arrow cycles are divisors in length of the current
            #cycle.
            coef = sum([l for l in arrow_cycle_len if powers%l == 0])
            polykey = (powers, coef)
            if polykey not in polynomials:
                polynomials[polykey] = 1
            else:
                polynomials[polykey] += 1
        result.append(polynomials)
        
    return result

def polya(concentrations, group, arrowings=None, debug=False):
    """Uses a group and concentrations to find the number of unique arrangements as described by 
    polya.
    
    :arg concentrations: specify a list of integers specifying how many of each coloring should
      be present in each of the enumerated lists.
    :arg group: group operations for permuting the colorings.
    """

    if arrowings is None:
        arrowings = [False]*len(concentrations)
    elif isinstance(arrowings, int):
        arrowings = [False]*(len(concentrations)-arrowings) + [True]*arrowings
        
    #This is to check that the concentrations sum to the number of sites the group is
    #operating on
    
    if sum(concentrations) != len(group[0][0]):
        raise ValueError("The concentrations don't sum to the number of things the group is acting on!")
            
    for k in range(2):
        for g in range(len(group)):
            if 0 not in group[g][k]:
                group[g][k] = [j-1 for j in group[g][k]]

    polyndict = {}
    #The operations in the group are used to construct the unique polynomials for each operation.
    for polynomials in _group_to_cyclic(group):
        #Construct a product of multinomials for this group operation.
        p = Product(1,concentrations)
        for exp, coeff in polynomials:
            p.multinoms.append(Multinomial(exp, coeff, arrowings, polynomials[(exp, coeff)]))

        key = str(p)
        if key not in polyndict:
            polyndict[key] = p
        else:
            polyndict[key].coefficient += 1

    if debug: #pragma: no cover
        for key in polyndict:
            print((str(polyndict[key]), " => ", polyndict[key].coeff()))
            
    pp = 0
    for p in list(polyndict.values()):
        pp += p.coeff()
        
    rad = sum([p.coeff() for p in list(polyndict.values())])
    return int(rad/float(len(group)))             

def _examples():
    """Print some examples on how to use this python version of the code."""
    helptext = ("For all the examples below, it is assumed that you know the fixed concentration "
                "term T in advance. This term is the first, *positional* argument to the script. "
                "In addition to the term T, you need to specify the group operations as permutation "
                "lists. They can be either zero- or one-based. Group operations can be specified "
                "with the group generators or as a 2D matrix with all the group operations; if the "
                "lists of values were saved directly from python using a __repr__ or __str__ then "
                "use the '-parse' argument to specify that.")
    egs = [("Find the Polya Coefficient with Group Generators",
            "The code below finds the number of unique ways to color a square with 4 corners using "
            "2 different colors such that there are 2 corners with each color. "
            "The group is specified using generators in a file called 'generators.in.paper'. The "
            "contents of the generators file are:\n  4 3 2 1\n  2 3 4 1\nand are the generators "
            "for the dihedral group of degree 4.", "./polya.py 2 2 -generators generators.in.paper"),
           ("Find the Polya Coefficient with an Entire Group",
            "This code also finds the coefficient, but for a larger group with 144 operations acting "
            "on a finite set with 20 elements. The term T is specified as [4,4,4,2,2,2,2] so that we "
            "want 4 of the first 3 colors and 2 of the last 4 colors with 7, the total number of "
            "colors in the enumeration. The group file 'group.out.cr6' can be viewed in the repo at "
            "'polya/fortran/tests/'.", "./polya.py 4 4 4 2 2 2 2 -group group.out.cr6")]

    print("POLYA ENUMERATION THEOREM SOLVER\n")
    for eg in egs:
        title, desc, code = eg
        print("--" + title + '--\n')
        print(desc + '\n')
        print('  ' + code + '\n')

def _parser_options():
    """Parses the options and arguments from the command line."""
    import argparse
    parser = argparse.ArgumentParser(description="Polya Coefficient Calculator")
    parser.add_argument("-generators",
                        help=("Specify the name/path to a file that lists the generators for "
                              "the symmetry group defining uniqueness on the lattice."))
    parser.add_argument("-group",
                        help=("Specify the name/path to a file listing the *entire* set of group "
                              "symmetry and arrow permutation operations defining uniqueness on the lattice."))
    parser.add_argument("-parse", choices=["python", "text"], default="text",
                        help=("Choose how the group files will be interpreted by the script:\n"
                              "- 'python': the text is assumed to be a valid python expression, \n"
                              "\tsuch as a list, and is interpreted using eval(). \n"
                              "- 'text': text values are split on whitespace and converted to \n"
                              "\tintegers. One group operation/generator per line."))
    parser.add_argument("concentrations", type=int, nargs="*", default=[0],
                        help=("The number of each type of coloring in the concentration restricted "
                              "enumeration on a lattice."))
    parser.add_argument("-debug", action="store_true",
                        help="Print verbose polya polynomial information for debugging.")
    parser.add_argument("-examples", action="store_true",
                        help="Print some examples for how to use the Polya solver.")
    parser.add_argument("-arrows", type=int, default=None
                        help="The numer of elements with 'arrows' on them in the system.")

    vardict = vars(parser.parse_args())
    if vardict["examples"]:
        _examples()
        exit(0)
    return vardict

def _read_file(args, filepath):
    """Parses the contents of the specified file using the 'parse' arguments from script args."""
    from os import path
    contents = []
    with open(path.expanduser(filepath)) as f:
        if args["parse"] == "text":
            for line in f:
                contents.append(map(int, line.split()))
        else:
            contents = eval(f.read())

    return contents

def script_polya(args):
    """Calculates the number of unique ways to enumerate a fixed set of colorings on a lattice
    subject to a set of symmetry operations.
    """
    if not args["generators"] and not args["group"]:
        raise ValueError("You must specify either the generators or the group.")

    if args["generators"]:
        gens = _read_file(args, args["generators"])
        grpops = group(gens)
    elif args["group"]:
        grpops = _read_file(args, args["group"])

    if args["arrows"]:
        coeff = polya(args["concentrations"], grpops, arrowings=args["arrows"], args["debug"])
    else:
        coeff = polya(args["concentrations"], grpops, args["debug"])
    
    print(coeff)
    return coeff

if __name__ == '__main__':
    script_polya(_parser_options())
