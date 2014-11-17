"""This program produces input files for testing the enumeration code. Due to the shortness of 
this project there is some code duplication that has occured from the main code and here. For 
some of the tests array elements are choosen usinchg a random number so if new tests are desired 
or required please be sure to remove the previously generated input files before generating the 
new ones or the unit test may fail.

Output files include:
coeffs.out.*
concs.in.*
generators.in.*
group.out.*
kidcount.*.out.*
mnseq.in.*
mnseq.len.in.*
multinomial.in.*
possibles.in.*
prod.out.*
prod.seq.*
result.out.*
seq.in.*
seqexp.*.out.*
sum_seq.out.*
where the *s refer to the case, or trest group, number.
"""
import unittest as ut
from polya import *
from polya import _group_to_cyclic
import cProfile
import time
from pstats import Stats
from testgen import save, load

class Timer:    
"""This class impliments a timer that is used to time and profile the code as it runs.
"""
    def __enter__(self):
        """Sets the start time of the process being timed"""
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        """Sets the end time of the code after it is done running then returns the difference
        between start and stop times as the runtime for the code.
        """
        self.end = time.clock()
        self.interval = self.end - self.start

class TestPolya(ut.TestCase):
    """Runs the polya python code on each identified test group and builds appropriate input and
    output files for the fortran code."""
    def setUp(self):
        """Sets up the test groups and the folder where the resultant testing files should be 
        saved for the fortran tests.
        self.folder is the directory that the files will be saved in.
        self.generators stores the test groups to be used in the following format:
           {"gens": [the generators for the group]
            "concs": [the concentrations the different colors, objects, being enumerated]
            "size": the number of operations in the group
            "result": the number of unique arrangements after enumeration.}
        Also initializes an array for profiling the code and dictionaries for the groups and 
        products to be stored in.
        """
        self.folder = "~/Documents/research/enum4/polya/fortran/tests"
        self.pr = cProfile.Profile()
        self.generators = [{"gens": [[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1],
                                     [2, 3, 4, 1, 6, 7, 8, 5, 10, 11, 12, 9, 14, 15, 16, 13]],
                            "concs": [4,4,4,2,2],
                            "size": 1024,
                            "result": 418812},
                           {"gens": [[2, 3, 4, 5, 6, 7, 8, 1],[8, 7, 6, 5, 4, 3, 2, 1]],
                            "concs": [4,4],
                            "size": 16,
                            "result": 8},
                           {"gens": [[2, 1, 4, 3, 5, 6, 7, 8, 9, 10, 11, 12],
                                     [4, 3, 2, 1, 5, 6, 7, 8, 9, 10, 11, 12],
                                     [1, 2, 3, 4, 6, 7, 8, 5, 9, 10, 11, 12],
                                     [1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 12],
                                     [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 9, 12],
                                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 10]],
                            "concs": [2,2,2,2,2,2],
                            "size": 1152,
                            "result": 11070},
                           {"gens": [[1, 2, 3, 4, 5, 6, 7, 8, 9],[7, 4, 1, 8, 5, 2, 9, 6, 3],
                                     [9, 8, 7, 6, 5, 4, 3, 2, 1],[3, 6, 9, 2, 5, 8, 1, 4, 7],
                                     [7, 8, 9, 4, 5, 6, 1, 2, 3],[3, 2, 1, 6, 5, 4, 9, 8, 7],
                                     [1, 4, 7, 2, 5, 8, 3, 6, 9],[9, 6, 3, 8, 5, 2, 7, 4, 1]],
                            "concs": [3,3,3],
                            "size": 8,
                            "result": 228},
                           {"gens": [[2, 3, 4, 5, 1, 6, 7, 8, 9 ,10],
                                     [6, 7, 8, 9, 10, 1, 2, 3, 4, 5]],
                            "concs": [3,3,3,1],
                            "size": 50,
                            "result": 336},
                           {"gens": [[2, 3, 4, 5, 1, 7, 8, 9, 10, 6],
                                     [2, 1, 3, 4, 5, 7, 6, 8, 9, 10]],
                            "concs": [2,2,2,2,2],
                            "size": 120,
                            "result": 1110},
                           {"gens": [[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
                                      13, 14, 15, 16, 17, 18, 19, 20, 1],
                                     [2, 3, 4, 1, 6, 7, 8, 5, 10, 11, 12, 
                                      9, 14, 15, 16, 13, 18, 19, 20, 17]],
                            "concs": [5,8,7],
                            "size": 2500,
                            "result": 43332},
                           {"gens": [[2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16],
                                     [1,2,3,4,6,7,8,5,9,10,11,12,13,14,15,16],
                                     [1,2,3,4,5,6,7,8,10,11,12,9,13,14,15,16],
                                     [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,13]],
                            "concs": [8,8],
                            "size": 256,
                            "result": 196},
                           {"gens": [[2,3,4,5,1,6,7,8,9,10],[2,1,4,3,6,5,8,7,10,9]],
                            "concs": [5,2,3],
                            "size": 120,
                            "result": 55},
                           {"gens": [[3,4,5,6,7,8,9,10,1,2],[3,4,1,2,5,6,7,8,9,10]],
                            "concs": [5,5],
                            "size": 120,
                            "result": 12},
                           {"gens": [[2,3,4,5,1,7,8,9,10,6],[2,1,3,4,5,7,6,8,9,10]],
                            "concs": [5,5],
                            "size": 120,
                            "result": 12},
                           {"gens": [[2,3,4,5,6,7,8,1],[8,7,6,5,4,3,2,1]],
                            "concs": [4,4],
                            "size": 16,
                            "result": 8},
                           {"gens": [[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1],
                                     [2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19]],
                            "concs": [5,5,10],
                            "size": 200,
                            "result": 234732},
                           {"gens": [[2,3,4,5,6,7,8,9,10,1],[10,9,8,7,6,5,4,3,2,1]],
                            "concs": [4,3,3],
                            "size": 20,
                            "result": 216},
                           {"gens": [[2,3,4,5,6,1,7,8,9,10,11,12],
                                     [1,2,3,4,5,6,8,9,7,10,11,12],
                                     [1,2,3,4,5,6,8,7,9,10,11,12],
                                     [1,2,3,4,5,6,7,8,9,11,12,10],
                                     [1,2,3,4,5,6,7,8,9,11,10,12]],
                            "concs": [4,4,4],
                            "size": 216,
                            "result": 634},
                           {"gens": [[2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16],
                                     [1,2,3,4,6,7,8,5,9,10,11,12,13,14,15,16],
                                     [1,2,3,4,5,6,7,8,10,11,12,9,13,14,15,16],
                                     [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,13]],
                            "concs": [3,3,4,6],
                            "size": 256,
                            "result": 170376},
                           {"gens": [[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1],
                                     [16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]],
                            "concs": [4,4,4,4],
                            "size": 32,
                            "result": 1972059},
                           {"gens": [[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1],
                                     [20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]],
                            "concs": [8,8,4],
                            "size": 40,
                            "result": 1560534},
                           {"gens": [[2,3,4,5,1,6,7,8,9,10],
                                     [6,7,8,9,10,1,2,3,4,5]],
                            "concs": [3,3,3,1],
                            "size": 50,
                            "result": 336},
                           {"gens": [[2,3,4,1,5,6,7,8],
                                     [2,1,3,4,5,6,7,8],
                                     [1,2,3,4,6,7,8,5],
                                     [1,2,3,4,6,5,7,8]],
                            "concs": [2,2,4],
                            "size": 576,
                            "result": 9},
                           {"gens": [[2,3,4,5,6,7,1,8,9,10,11,12,13,14,15,16,17,18,19,20],
                                     [1,2,3,4,5,6,7,9,10,11,12,13,14,8,15,16,17,18,19,20],
                                     [1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,15,18,19,20],
                                     [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,18]],
                            "concs": [6,7,7],
                            "size": 441,
                            "result": 405720},
                           {"gens": [[2,3,4,5,1,6],
                                     [1,2,3,5,6,4]],
                            "concs": [1,1,1,1,1,1],
                            "size": 360,
                            "result": 2}]

        self.groups = {}
        self.products = {}
        self.total = 0

    def gen_group(self, index):
        """Tests the group generation from the generators. The method outputs generators.in.* 
        and group.out.*"""
        if index not in self.groups:
            single = self.generators[index]
            gens = single["gens"]
            save([gens],["generators"],self.folder,case=index)
            grpops = group(gens)
            self.assertEqual(len(grpops), single["size"], "{} != {}".format(len(grpops), 
                                                                            single["size"]))
            #Save the input/output files for the fortran unit tests.
            save([grpops],["group"],self.folder,infile=False,case=index)
            self.groups[index] = grpops

    def gen_product(self, index):
        """Tests the construction of a Product of multinomials. The method outputs 
        multinomial.in.* """
        self.gen_group(index)
        #We only need to test the first operation in each group since we are testing the class
        #initializer.
        cyclic = _group_to_cyclic(self.groups[index],(0,1))[0]        
        p = Product(1, self.generators[index]["concs"])
        marray = []
        for exp in cyclic:
            #turns the product into a multinomial class
            p.multinoms.append(Multinomial(exp, cyclic[exp]))
            marray.append([exp, cyclic[exp]])
        self.products[index] = p
        #Save the input/output files for the fortran unit tests.
        save([marray], ["multinomial"], self.folder, case=index)

    def gen_sequence(self, index):
        """Tests the Sequence class of Polya's initialization and expand method. The method
        outputs possibles.in.*, seqexp.*.out.*, kidcount.*.out.*, mnseq.in.*, and coeffs.out.*"""
        from itertools import product
        self.gen_group(index)
        #We only need to test the first operation in each group since later unit tests
        #will confirm the totals for the entire group.
        cyclic = _group_to_cyclic(self.groups[index],(0,1))[0]
        keys = list(cyclic.keys())
        possibles = [list(range(0,power*cyclic[power]+1, power)) for power in keys]
        seqbase = [s for s in product(*possibles) if sum(s) == self.generators[index]["concs"][0]]
        powersums = [k*cyclic[k] for k in keys]
        #Save the input/output files for the fortran unit tests.
        save([possibles],["possibles"],self.folder,case=index)

        #We choose a random number here to ensure that the sequences that are selected are not
        #all trivial, ie all zeros, producing a better test group for the expand method.
        #Here we also have some code duplication from the polya code to accomadate the radom
        #number implementation.
        from random import randint
        if len(seqbase) > 0:
            irand = randint(0, len(seqbase)-1)
            #choose a random sequence for the unittest
            seq = seqbase[irand]
            key = keys[0]
            mnseq = []
            save([seq],["seq"],self.folder,case=index)
            for i in range(len(seq)):
                #For each sequence we calculate a powersum. The sequence intances are then 
                #contsructed and expanded.
                powersum = key*cyclic[key]
                varseq = Sequence(seq[i], possibles[i], 1, powersum, 
                                  self.generators[index]["concs"])
                varexp = varseq.expand()
                #append each expanded sequence to the mnseq array
                mnseq.append(varexp)
                #Save the input/output files for the fortran unit tests.
                save([varexp, varseq.kidcount], ["seqexp.{}".format(i), "kidcount.{}".format(i)],
                     self.folder,case=index,infile=False)

            #Save the input/output files for the fortran unit tests.
            save(mnseq,["mnseq.in.{}".format(i+1) for i in range(len(mnseq))],self.folder, 
                 infile=None, case=index)
            save([len(mnseq)], ["mnseq.len"], self.folder, case=index)

            #add this to our dictionary of products if it's not there yet
            if index not in self.products:
                self.gen_product(index)
            p = self.products[index]
            sumseq = int(p._sum_sequences(mnseq))
            
            #Save the input/output files for the fortran unit tests.
            save([sumseq], ["sum_seq"], self.folder, infile=False, case=index)

        coeffs = 0
        if index not in self.products:
            self.gen_product(index)
        p = self.products[index]
        for seq in seqbase:
            mnseq = []
            #Each sequence calculated from the first variable has an entry for each multinomial
            #in this product. The Sequence instances construct smart sequences for the remaining
            #variables in each multinomial separately
            for i in range(len(seq)):
                varseq = Sequence(seq[i], possibles[i], 1, powersums[i], 
                                  self.generators[index]["concs"])
                mnseq.append(varseq.expand())
            coeffs += int(p._sum_sequences(mnseq))

        #Save the input/output files for the fortran unit tests.
        save([coeffs],["coeffs"], self.folder, infile=False, case=index)

    def gen_test(self, index):
        """ Tests the polya code for each of the given test cases. The method outputs concs.in.*
        and results.out.*"""
        single = self.generators[index]
        self.gen_group(index)
        grpops = self.groups[index]

        #Save the input/output files for the fortran unit tests.
        save([single["concs"]],["concs"],self.folder,case=index)
        save([single["result"]],["result"],self.folder,case=index,infile=False)

        #Time the polya code while getting the output for testing
        with Timer() as t:
            self.pr.enable()
            coeff = polya(single["concs"], grpops)
            self.pr.disable()
        self.total += t.interval
        tstr = 'in %.05f sec.' % t.interval

        #print the result and the generators to the screen.
        print("{}x{} => {} {}".format(len(single["gens"]), len(single["gens"][0]), 
                                      single["result"], tstr))
        self.assertEqual(coeff, single["result"], 
                         "{} => {} != {}".format(single["gens"], single["result"], coeff))

    def test_groups(self):
        """Builds an array of the group for testing"""
        for i in range(len(self.generators)):
            self.gen_group(i)

    def test_sequences(self):
        """Builds an array of the sequences for testing."""
        for i in range(len(self.generators)):
            self.gen_sequence(i)

    def test_products(self):
        """Builds an array the products for testing."""
        for i in range(len(self.generators)):
            self.gen_product(i)

    def test_coefficients(self):
        """Runs a test of the code and checks the results then prints the accurracy, as 
        percent correct, to the screen"""
        for i in range(len(self.generators)):
           self.gen_test(i)
        print("TOTAL: %.06f" % self.total)

