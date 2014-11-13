import unittest as ut
from polya import *
from polya import _group_to_cyclic
import cProfile
import time
from pstats import Stats
from testgen import save, load

class Timer:    
    def __enter__(self):
        self.start = time.clock()
        return self

    def __exit__(self, *args):
        self.end = time.clock()
        self.interval = self.end - self.start

class TestPolya(ut.TestCase):
    def setUp(self):
        self.folder = "~/codes/projects/polya/fortran/tests"
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
                            "result": 196}]
        self.groups = {}
        self.products = {}
        self.total = 0

    def gen_group(self, index):
        """Tests the group generation from the generators."""
        if index not in self.groups:
            single = self.generators[index]
            gens = single["gens"]
            save([gens],["generators"],self.folder,case=index)
            grpops = group(gens)
            self.assertEqual(len(grpops), single["size"], "{} != {}".format(len(grpops), single["size"]))
            save([grpops],["group"],self.folder,infile=False,case=index)
            self.groups[index] = grpops

    def gen_product(self, index):
        """Tests the construction of a Product of multinomials."""
        self.gen_group(index)
        #We only need to test the first operation in each group since we are testing the class
        #initializer.
        cyclic = _group_to_cyclic(self.groups[index],(0,1))[0]        
        p = Product(1, self.generators[index]["concs"])
        marray = []
        for exp in cyclic:
            p.multinoms.append(Multinomial(exp, cyclic[exp]))
            marray.append([exp, cyclic[exp]])
        self.products[index] = p
        save([marray], ["multinomial"], self.folder, case=index)

    def gen_sequence(self, index):
        """Tests the Sequence class of Polya's initialization and expand method."""
        from itertools import product
        self.gen_group(index)
        #We only need to test the first operation in each group since later unit tests
        #will confirm the totals for the entire group.
        cyclic = _group_to_cyclic(self.groups[index],(0,1))[0]
        keys = list(cyclic.keys())
        possibles = [list(range(0,power*cyclic[power]+1, power)) for power in keys]
        seqbase = [s for s in product(*possibles) if sum(s) == self.generators[index]["concs"][0]]
        powersums = [k*cyclic[k] for k in keys]
        save([possibles],["possibles"],self.folder,case=index)

        from random import randint
        if len(seqbase) > 0:
            irand = randint(0, len(seqbase)-1)
            seq = seqbase[irand]
            key = keys[0]
            mnseq = []
            save([seq],["seq"],self.folder,case=index)
            for i in range(len(seq)):
                powersum = key*cyclic[key]
                varseq = Sequence(seq[i], possibles[i], 1, powersum, self.generators[index]["concs"])
                varexp = varseq.expand()
                mnseq.append(varexp)
                save([varexp, varseq.kidcount], ["seqexp.{}".format(i), "kidcount.{}".format(i)],
                     self.folder,case=index,infile=False)

            save(mnseq,["mnseq.in.{}".format(i+1) for i in range(len(mnseq))],self.folder, 
                 infile=None, case=index)
            save([len(mnseq)], ["mnseq.len"], self.folder, case=index)
            if index not in self.products:
                self.gen_product(index)
            p = self.products[index]
            sumseq = int(p._sum_sequences(mnseq))
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
                varseq = Sequence(seq[i], possibles[i], 1, powersums[i], self.generators[index]["concs"])
                mnseq.append(varseq.expand())
            coeffs += int(p._sum_sequences(mnseq))

        save([coeffs],["coeffs"], self.folder, infile=False, case=index)

    def gen_test(self, index):
        single = self.generators[index]
        self.gen_group(index)
        grpops = self.groups[index]

        #Save the input/output files for the fortran unit tests.
        save([single["concs"]],["concs"],self.folder,case=index)
        save([single["result"]],["result"],self.folder,case=index,infile=False)

        with Timer() as t:
            self.pr.enable()
            coeff = polya(single["concs"], grpops)
            self.pr.disable()
        self.total += t.interval
        tstr = 'in %.05f sec.' % t.interval
        print("{}x{} => {} {}".format(len(single["gens"]), len(single["gens"][0]), single["result"], tstr))
        self.assertEqual(coeff, single["result"], 
                         "{} => {}".format(single["gens"], single["result"]))

    @ut.skip("Once-off Okay!")
    def test_groups(self):
        for i in range(len(self.generators)):
            self.gen_group(i)

    def test_sequences(self):
        for i in range(len(self.generators)):
            self.gen_sequence(i)

    @ut.skip("Once-off Okay!")
    def test_products(self):
        for i in range(len(self.generators)):
            self.gen_product(i)

    def test_coefficients(self):
        for i in range(len(self.generators)):
           self.gen_test(i)
        print("TOTAL: %.06f" % self.total)

        # p = Stats(self.pr)
        # p.strip_dirs()
        # p.sort_stats('cumtime')
        # p.print_stats()
