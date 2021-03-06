#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
AUTHORS: Conrad W. Rosenbrock, Wiley S. Morgan (October 2014)

Classes to support the calculation of coefficients for specific terms in a product
of multinomials. Construct a product class by specifying the exponent and target term
and then add multinomials using the Product instance's append(). The coefficient is
then available from the coeff().

EXAMPLE: find the coefficient of x^4.y^4.z^4 in 3*(x+y+z)^3*(x^2+y^2+z^2)^3*(x^3+y^3+z^3).
ANSWER: 162 from Mathematica

>> p = Product(3, [4,4,4])
>> p.append(Multinomial([1,1,1],3))
>> p.append(Multinomial([2,2,2],3))
>> p.append(Multinomial([3,3,3],1))
>> print(p.coeff())
162
"""
from numpy import array, sum as nsum

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
        if i < len(possibles):
            #Filter the possible values for the variable being considered based on the
            #exponent of the multinomial. When multinomials are expanded, the sum of
            #the exponents in any term must be less than the exponent on the multinomial
            #times the maximum power of any of its (unexpanded) terms.

            #We find all the possible values for this variable by ensuring that:
            # 1) it's exponent is compatible with the exponents of all variables to the left of it.
            # 2) the exponent we are suggesting is in the list of possible values for the variable.
            # 3) the exponent remains positive.
            self.kids = [Sequence(p-root, possibles, i+1, powersum, targets, self) 
                         for p in possibles[i] if p-root >= 0
                         and p-root <= targets[i]
                         and abs(p - root) <= powersum-self.used 
                         and abs(p-self.used) % possibles[i][1] == 0]
        else:
            self.kids = []

    def expand(self):
        """Recursively generates a list of all relevant sequences for this multinomial term."""
        #Iterate through the child sequences and add their variable root values if
        #the total sequence sums to the target.
        sequences = []
        for kid in self.kids:
            for seq in kid.expand():
                #Here is where the recursion happens; we add the sequence of this variable's
                #children to the right of this root.
                sequences.append((self._root,) + seq)

        if len(self.kids) == 0:
            return [(self._root,)]
        else:
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
            if all([p-t>0 for p,t in zip(self.multinoms[0].powers, self.targets)]):
                return 0
            else:
                return self.multinoms[0].nchoosekm(self.targets)*self.coefficient
        
        from itertools import product
        #Get a list of the possible exponents for each variable in each of the multinomials.
        #We start with the first variable and choose only those combinations of exponents
        #across *all* the multinomials that give the correct target exponent for that variable.
        possibles = [n.possible_powers for n in self.multinoms]
        possfirst = [p[0] for p in possibles]
        seq0 = [s for s in product(*possfirst) if sum(s) == self.targets[0]]

        #Next, we construct Sequence instances for each of the first variable compatible
        #possibilities and follow them through to the other variables.
        sequences = []
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
        sortedmns = sorted(self.multinoms, key=(lambda m: (m.exponent,m.powers[0])), reverse=True)
        return str(self.coefficient) + ''.join([str(mn) for mn in sortedmns])            

class Multinomial(object):
    """Represents a multinomial expansion."""
    def __init__(self, powers, exponent=1):
        """Sets up the multinomial.

        :arg powers: the powers on each of the *unexpanded* variables in the multinomial;
          of the form (x^2+y^2+z^4) => [2,2,4]. There should be an item in the list for each
          of the variables in the multinomial. No support is provided for coefficients on the
          variables inside the brackets (only Polya implementation).
        :arg exponent: the exponent of the entire multinomial.
        """
        self.powers = powers
        self.exponent = exponent
        self.powersum = max(powers)*exponent
        """Returns the integer value that all term exponents in the multinomial should
        sum to (or be less than)."""
        self.possible_powers = [list(range(0,p*exponent+1, p)) for p in powers]
        """For each variable being considered, determines the possible powers based
        on the exponent in the multinomial."""

    def __str__(self):
        #We want to print the multinomial out in a nice, readable way, similar to how
        #they are presented in Mathematica.
        contents = ' + '.join(["{}".format(p) for p in self.powers])
        return "({})^{}".format(contents, self.exponent)

    def normed_seq(self, seq):
        """Normalizes the specified sequence using the powers of unexpanded terms in the multinomial.
        
        :arg seq: a list of exponents in an *expanded* term.
        """
        return [int(ai/bi) for ai, bi in zip(seq, self.powers)]

    def nchoosekm(self, sequence):
        """Returns the number of different ways to partition an n-element
        set into disjoint subsets of sizes k1, ..., km.

        :arg sequence: an un-normed tuple of form (k1, k2, k3).
        """
        prod = 1
        if not all([sequence[i]%self.powers[i] == 0 for i in range(len(sequence))]):
            return 0
        else:
            normseq = self.normed_seq(sequence)
            for i in range(len(sequence)):
                nsum = sum(normseq[0:i+1])
                prod *= Multinomial.nchoosek(nsum, normseq[i])

            return prod
        
    @staticmethod
    def nchoosek(n, k):
        """This implementation was taken from "Binomial Coefﬁcient Computation: Recursion 
        or Iteration?" by Yannis Manolopoulos, ACM SIGCSE Bulletin InRoads, Vol.34, No.4, 
        December 2002. http://delab.csd.auth.gr/papers/SBI02m.pdf It is supposed to be robust 
        against large, intermediate values and to have optimal complexity.
        """
        if k < 0 or k > n:
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

def group(gen):
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
        ngens = [list(map(lambda e: e-1, g)) for g in gen]
    else:
        ngens = gen

    groupi = []
    for i in ngens:
        for j in filter(lambda k: k!=i, ngens):
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

def polya(concentrations, group, debug=False):
    """Uses a group and concentrations to find the number of unique arrangements as described by 
    polya.
    
    :arg concentrations: specify a list of integers specifying how many of each coloring should
      be present in each of the enumerated lists.
    :arg group: group operations for permuting the colorings.
    """
    polyndict = {}
    #The operations in the group are used to construct the unique polynomials for each operation.
    for operation in group:
        #visited has the same # of elements as the group operation and is used to make sure each
        #element in the array is visited as we loop through in a *non-sequential* order.
        visited = [0]*len(operation)
        polynomials = {}

        while 0 in visited:
            #Start with the first element in the group that hasn't been visited yet. The first
            #non-trivial polynomials have powers > 0.
            cursor = vindex = visited.index(0)
            powers = 1 
            #change the current position to having been visited; move the cursor.
            visited[cursor] = 1
            cursor = operation[cursor]

            #The power of the variables in the polynomials is equal to the number of group operations
            #separating the cursor's current position from its *value* in the group operations list.
            while cursor != vindex:
                visited[cursor] = 1
                powers += 1
                cursor = operation[cursor]

            #We now have everything need to construct part of the polynomial. This is done by taking powers
            #and using it to construct an array of length equal to the number of elements in the system 
            #each entry in the array is set to be equal to powers.
            partpolyn = (powers,)*len(concentrations)
            if partpolyn not in polynomials:
                polynomials[partpolyn] = 1
            else:
                polynomials[partpolyn] += 1

        #Construct a product of multinomials for this group operation.
        p = Product(1,concentrations)
        for explist in polynomials:
            p.multinoms.append(Multinomial(explist, polynomials[explist]))

        key = str(p)
        if key not in polyndict:
            polyndict[key] = p
        else:
            polyndict[key].coefficient += 1

    if debug:
        for key in polyndict:
            print(str(polyndict[key]), " => ", polyndict[key].coeff())

    rad = sum([p.coeff() for p in polyndict.values()])
    return int(rad/float(len(group)))             

def _parser_options():
    """Parses the options and arguments from the command line."""
    import argparse
    parser = argparse.ArgumentParser(description="Polya Coefficient Calculator")
    parser.add_argument("-generators",
                        help=("Specify the name/path to a file that lists the generators for "
                              "the symmetry group defining uniqueness on the lattice."))
    parser.add_argument("-group",
                        help=("Specify the name/path to a file listing the *entire* set of group "
                              "symmetry operations defining uniqueness on the lattice."))
    parser.add_argument("-parse", choices=["python", "text"], default="text",
                        help=("Choose how the group files will be interpreted by the script:\n"
                              "- 'python': the text is assumed to be a valid python expression, \n"
                              "\tsuch as a list, and is interpreted using eval(). \n"
                              "- 'text': text values are split on whitespace and converted to \n"
                              "\tintegers. One group operation/generator per line."))
    parser.add_argument("concentrations", type=int, nargs="+",
                        help=("The number of each type of coloring in the concentration restricted "
                              "enumeration on a lattice."))
    parser.add_argument("-debug", action="store_true",
                        help="Print verbose polya polynomial information for debugging.")
    return vars(parser.parse_args())

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

    coeff = polya(args["concentrations"], grpops, args["debug"])
    print(coeff)
    return coeff

if __name__ == '__main__':
    script_polya(_parser_options())
