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
    """This class impliments a timer that is used to time and profile the code as it runs."""
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
        self.generators = [{"fname": 0,
                            "gens": [[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1],
                                     [2, 3, 4, 1, 6, 7, 8, 5, 10, 11, 12, 9, 14, 15, 16, 13]],
                            "concs": [4,4,4,2,2],
                            "size": 1024,
                            "result": 418812},
                           {"fname": 1,
                            "gens": [[2, 3, 4, 5, 6, 7, 8, 1],[8, 7, 6, 5, 4, 3, 2, 1]],
                            "concs": [4,4],
                            "size": 16,
                            "result": 8},
                           {"fname": 2,
                            "gens": [[2, 1, 4, 3, 5, 6, 7, 8, 9, 10, 11, 12],
                                     [4, 3, 2, 1, 5, 6, 7, 8, 9, 10, 11, 12],
                                     [1, 2, 3, 4, 6, 7, 8, 5, 9, 10, 11, 12],
                                     [1, 2, 3, 4, 6, 5, 7, 8, 9, 10, 11, 12],
                                     [1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 9, 12],
                                     [1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 10]],
                            "concs": [2,2,2,2,2,2],
                            "size": 1152,
                            "result": 11070},
                           {"fname": 3,
                            "gens": [[1, 2, 3, 4, 5, 6, 7, 8, 9],[7, 4, 1, 8, 5, 2, 9, 6, 3],
                                     [9, 8, 7, 6, 5, 4, 3, 2, 1],[3, 6, 9, 2, 5, 8, 1, 4, 7],
                                     [7, 8, 9, 4, 5, 6, 1, 2, 3],[3, 2, 1, 6, 5, 4, 9, 8, 7],
                                     [1, 4, 7, 2, 5, 8, 3, 6, 9],[9, 6, 3, 8, 5, 2, 7, 4, 1]],
                            "concs": [3,3,3],
                            "size": 8,
                            "result": 228},
                           {"fname": 4,
                            "gens": [[2, 3, 4, 5, 1, 6, 7, 8, 9 ,10],
                                     [6, 7, 8, 9, 10, 1, 2, 3, 4, 5]],
                            "concs": [3,3,3,1],
                            "size": 50,
                            "result": 336},
                           {"fname": 5,
                            "gens": [[2, 3, 4, 5, 1, 7, 8, 9, 10, 6],
                                     [2, 1, 3, 4, 5, 7, 6, 8, 9, 10]],
                            "concs": [2,2,2,2,2],
                            "size": 120,
                            "result": 1110},
                           {"fname": 6,
                            "gens": [[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
                                      13, 14, 15, 16, 17, 18, 19, 20, 1],
                                     [2, 3, 4, 1, 6, 7, 8, 5, 10, 11, 12, 
                                      9, 14, 15, 16, 13, 18, 19, 20, 17]],
                            "concs": [5,8,7],
                            "size": 2500,
                            "result": 43332},
                           {"fname": 7,
                            "gens": [[2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16],
                                     [1,2,3,4,6,7,8,5,9,10,11,12,13,14,15,16],
                                     [1,2,3,4,5,6,7,8,10,11,12,9,13,14,15,16],
                                     [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,13]],
                            "concs": [8,8],
                            "size": 256,
                            "result": 196},
                           {"fname": 8,
                            "gens": [[2,3,4,5,1,6,7,8,9,10],[2,1,4,3,6,5,8,7,10,9]],
                            "concs": [5,2,3],
                            "size": 120,
                            "result": 55},
                           {"fname": 9,
                            "gens": [[3,4,5,6,7,8,9,10,1,2],[3,4,1,2,5,6,7,8,9,10]],
                            "concs": [5,5],
                            "size": 120,
                            "result": 12},
                           {"fname": 10,
                            "gens": [[2,3,4,5,1,7,8,9,10,6],[2,1,3,4,5,7,6,8,9,10]],
                            "concs": [5,5],
                            "size": 120,
                            "result": 12},
                           {"fname": 11,
                            "gens": [[2,3,4,5,6,7,8,1],[8,7,6,5,4,3,2,1]],
                            "concs": [4,4],
                            "size": 16,
                            "result": 8},
                           {"fname": 12,
                            "gens": [[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1],
                                     [2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,20,19]],
                            "concs": [5,5,10],
                            "size": 200,
                            "result": 234732},
                           {"fname": 13,
                            "gens": [[2,3,4,5,6,7,8,9,10,1],[10,9,8,7,6,5,4,3,2,1]],
                            "concs": [4,3,3],
                            "size": 20,
                            "result": 216},
                           {"fname": 14,
                            "gens": [[2,3,4,5,6,1,7,8,9,10,11,12],
                                     [1,2,3,4,5,6,8,9,7,10,11,12],
                                     [1,2,3,4,5,6,8,7,9,10,11,12],
                                     [1,2,3,4,5,6,7,8,9,11,12,10],
                                     [1,2,3,4,5,6,7,8,9,11,10,12]],
                            "concs": [4,4,4],
                            "size": 216,
                            "result": 634},
                           {"fname": 15,
                            "gens": [[2,3,4,1,5,6,7,8,9,10,11,12,13,14,15,16],
                                     [1,2,3,4,6,7,8,5,9,10,11,12,13,14,15,16],
                                     [1,2,3,4,5,6,7,8,10,11,12,9,13,14,15,16],
                                     [1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,13]],
                            "concs": [3,3,4,6],
                            "size": 256,
                            "result": 170376},
                           {"fname": 16,
                            "gens": [[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,1],
                                     [16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]],
                            "concs": [4,4,4,4],
                            "size": 32,
                            "result": 1972059},
                           {"fname": 17,
                            "gens": [[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1],
                                     [20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1]],
                            "concs": [8,8,4],
                            "size": 40,
                            "result": 1560534},
                           {"fname": 18,
                            "gens": [[2,3,4,5,1,6,7,8,9,10],
                                     [6,7,8,9,10,1,2,3,4,5]],
                            "concs": [3,3,3,1],
                            "size": 50,
                            "result": 336},
                           {"fname": 19,
                            "gens": [[2,3,4,1,5,6,7,8],
                                     [2,1,3,4,5,6,7,8],
                                     [1,2,3,4,6,7,8,5],
                                     [1,2,3,4,6,5,7,8]],
                            "concs": [2,2,4],
                            "size": 576,
                            "result": 9},
                           {"fname": 20,
                            "gens": [[2,3,4,5,6,7,1,8,9,10,11,12,13,14,15,16,17,18,19,20],
                                     [1,2,3,4,5,6,7,9,10,11,12,13,14,8,15,16,17,18,19,20],
                                     [1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,15,18,19,20],
                                     [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,18]],
                            "concs": [6,7,7],
                            "size": 441,
                            "result": 405720},
                           {"fname": 21,
                            "gens": [[2,3,4,5,1,6],
                                     [1,2,3,5,6,4]],
                            "concs": [1,1,1,1,1,1],
                            "size": 360,
                            "result": 2},
                           {"fname": 'fg1',
                              "gens": [[6, 7, 8, 5, 10, 11, 12, 9, 2, 3, 4, 1], 
                                     [6, 7, 8, 5, 2, 3, 4, 1, 10, 11, 12, 9], 
                                     [6, 5, 7, 8, 10, 9, 11, 12, 2, 1, 3, 4], 
                                     [6, 5, 7, 8, 2, 1, 3, 4, 10, 9, 11, 12]],
                            "concs": [6,6],
                            "size": 144,
                            "result": 17},
                           {"fname": 'fg2',
                            "gens": 
                            [[8, 9, 7, 11, 12, 10, 14, 15, 13, 17, 18, 16, 2, 3, 1, 5, 6, 4], 
                             [8, 9, 7, 11, 12, 10, 2, 3, 1, 5, 6, 4, 14, 15, 13, 17, 18, 16], 
                             [7, 9, 11, 10, 12, 8, 13, 15, 17, 16, 18, 14, 1, 3, 5, 4, 6, 2], 
                             [7, 9, 11, 10, 12, 8, 1, 3, 5, 4, 6, 2, 13, 15, 17, 16, 18, 14]],
                            "concs": [9,9],
                            "size": 144,
                            "result": 378},
                           {"fname": 'fg3',
                            "gens": [[6, 7, 8, 5, 10, 11, 12, 9, 2, 3, 4, 1, 14, 15, 16, 13, 18, 19, 20, 17], [6, 7, 8, 5, 2, 3, 4, 1, 10, 11, 12, 9, 18, 19, 20, 17, 14, 15, 16, 13], [6, 5, 7, 8, 10, 9, 11, 12, 2, 1, 3, 4, 14, 13, 15, 16, 18, 17, 19, 20], [6, 5, 7, 8, 2, 1, 3, 4, 10, 9, 11, 12, 18, 17, 19, 20, 14, 13, 15, 16]],
                            "concs": [10,10],
                            "size": 144,
                            "result": 1562},
                           {"fname": 'fg4',
                            "gens": [[10, 11, 12, 9, 14, 15, 16, 13, 18, 19, 20, 17, 22, 23, 24, 21, 2, 3, 4, 1, 6, 7, 8, 5], [10, 11, 12, 9, 14, 15, 16, 13, 2, 3, 4, 1, 6, 7, 8, 5, 18, 19, 20, 17, 22, 23, 24, 21], [9, 12, 16, 13, 10, 11, 15, 14, 17, 20, 24, 21, 18, 19, 23, 22, 1, 4, 8, 5, 2, 3, 7, 6], [9, 12, 16, 13, 10, 11, 15, 14, 1, 4, 8, 5, 2, 3, 7, 6, 17, 20, 24, 21, 18, 19, 23, 22]],
                            "concs": [12,12],
                            "size": 144,
                            "result": 19219},
                           {"fname": 'fg5',
                            "gens": [[8, 9, 7, 11, 12, 10, 14, 15, 13, 17, 18, 16, 2, 3, 1, 5, 6, 4, 20, 21, 19, 23, 24, 22, 26, 27, 25, 29, 30, 28], [8, 9, 7, 11, 12, 10, 2, 3, 1, 5, 6, 4, 14, 15, 13, 17, 18, 16, 26, 27, 25, 29, 30, 28, 20, 21, 19, 23, 24, 22], [7, 9, 11, 10, 12, 8, 13, 15, 17, 16, 18, 14, 1, 3, 5, 4, 6, 2, 19, 21, 23, 22, 24, 20, 25, 27, 29, 28, 30, 26], [7, 9, 11, 10, 12, 8, 1, 3, 5, 4, 6, 2, 13, 15, 17, 16, 18, 14, 25, 27, 29, 28, 30, 26, 19, 21, 23, 22, 24, 20]],
                            "concs": [15,15],
                            "size": 144,
                            "result": 1081166},
                           {"fname": 'fg6',
                            "gens": [[10, 11, 12, 9, 14, 15, 16, 13, 18, 19, 20, 17, 22, 23, 24, 21, 2, 3, 4, 1, 6, 7, 8, 5, 26, 27, 28, 25, 30, 31, 32, 29, 34, 35, 36, 33, 38, 39, 40, 37], [10, 11, 12, 9, 14, 15, 16, 13, 2, 3, 4, 1, 6, 7, 8, 5, 18, 19, 20, 17, 22, 23, 24, 21, 34, 35, 36, 33, 38, 39, 40, 37, 26, 27, 28, 25, 30, 31, 32, 29], [9, 12, 16, 13, 10, 11, 15, 14, 17, 20, 24, 21, 18, 19, 23, 22, 1, 4, 8, 5, 2, 3, 7, 6, 25, 28, 32, 29, 26, 27, 31, 30, 33, 36, 40, 37, 34, 35, 39, 38], [9, 12, 16, 13, 10, 11, 15, 14, 1, 4, 8, 5, 2, 3, 7, 6, 17, 20, 24, 21, 18, 19, 23, 22, 33, 36, 40, 37, 34, 35, 39, 38, 25, 28, 32, 29, 26, 27, 31, 30]],
                            "concs": [20,20],
                            "size": 144,
                            "result": 957370938},
                           {"fname": 'fg7',
                            "gens": [[26, 27, 28, 25, 45, 46, 47, 48, 36, 33, 34, 35, 41, 42, 43, 44, 31, 32, 29, 30, 39, 40, 37, 38, 50, 51, 52, 49, 69, 70, 71, 72, 60, 57, 58, 59, 65, 66, 67, 68, 55, 56, 53, 54, 63, 64, 61, 62, 2, 3, 4, 1, 21, 22, 23, 24, 12, 9, 10, 11, 17, 18, 19, 20, 7, 8, 5, 6, 15, 16, 13, 14], [26, 27, 28, 25, 45, 46, 47, 48, 36, 33, 34, 35, 41, 42, 43, 44, 31, 32, 29, 30, 39, 40, 37, 38, 2, 3, 4, 1, 21, 22, 23, 24, 12, 9, 10, 11, 17, 18, 19, 20, 7, 8, 5, 6, 15, 16, 13, 14, 50, 51, 52, 49, 69, 70, 71, 72, 60, 57, 58, 59, 65, 66, 67, 68, 55, 56, 53, 54, 63, 64, 61, 62], [41, 42, 43, 44, 26, 27, 28, 25, 45, 46, 47, 48, 36, 33, 34, 35, 32, 29, 30, 31, 38, 39, 40, 37, 65, 66, 67, 68, 50, 51, 52, 49, 69, 70, 71, 72, 60, 57, 58, 59, 56, 53, 54, 55, 62, 63, 64, 61, 17, 18, 19, 20, 2, 3, 4, 1, 21, 22, 23, 24, 12, 9, 10, 11, 8, 5, 6, 7, 14, 15, 16, 13], [41, 42, 43, 44, 26, 27, 28, 25, 45, 46, 47, 48, 36, 33, 34, 35, 32, 29, 30, 31, 38, 39, 40, 37, 17, 18, 19, 20, 2, 3, 4, 1, 21, 22, 23, 24, 12, 9, 10, 11, 8, 5, 6, 7, 14, 15, 16, 13, 65, 66, 67, 68, 50, 51, 52, 49, 69, 70, 71, 72, 60, 57, 58, 59, 56, 53, 54, 55, 62, 63, 64, 61]],
                            "concs": [36,36],
                            "size": 144,
                            "result": 3073004178281331200},
                           {"fname": 'fg8',
                            "gens": [[26, 27, 28, 25, 45, 46, 47, 48, 36, 33, 34, 35, 41, 42, 43, 44, 31, 32, 29, 30, 39, 40, 37, 38, 50, 51, 52, 49, 69, 70, 71, 72, 60, 57, 58, 59, 65, 66, 67, 68, 55, 56, 53, 54, 63, 64, 61, 62, 2, 3, 4, 1, 21, 22, 23, 24, 12, 9, 10, 11, 17, 18, 19, 20, 7, 8, 5, 6, 15, 16, 13, 14, 74, 75, 76, 73, 93, 94, 95, 96, 84, 81, 82, 83, 89, 90, 91, 92, 79, 80, 77, 78, 87, 88, 85, 86, 98, 99, 100, 97, 117, 118, 119, 120, 108, 105, 106, 107, 113, 114, 115, 116, 103, 104, 101, 102, 111, 112, 109, 110], [26, 27, 28, 25, 45, 46, 47, 48, 36, 33, 34, 35, 41, 42, 43, 44, 31, 32, 29, 30, 39, 40, 37, 38, 2, 3, 4, 1, 21, 22, 23, 24, 12, 9, 10, 11, 17, 18, 19, 20, 7, 8, 5, 6, 15, 16, 13, 14, 50, 51, 52, 49, 69, 70, 71, 72, 60, 57, 58, 59, 65, 66, 67, 68, 55, 56, 53, 54, 63, 64, 61, 62, 98, 99, 100, 97, 117, 118, 119, 120, 108, 105, 106, 107, 113, 114, 115, 116, 103, 104, 101, 102, 111, 112, 109, 110, 74, 75, 76, 73, 93, 94, 95, 96, 84, 81, 82, 83, 89, 90, 91, 92, 79, 80, 77, 78, 87, 88, 85, 86], [41, 42, 43, 44, 26, 27, 28, 25, 45, 46, 47, 48, 36, 33, 34, 35, 32, 29, 30, 31, 38, 39, 40, 37, 65, 66, 67, 68, 50, 51, 52, 49, 69, 70, 71, 72, 60, 57, 58, 59, 56, 53, 54, 55, 62, 63, 64, 61, 17, 18, 19, 20, 2, 3, 4, 1, 21, 22, 23, 24, 12, 9, 10, 11, 8, 5, 6, 7, 14, 15, 16, 13, 89, 90, 91, 92, 74, 75, 76, 73, 93, 94, 95, 96, 84, 81, 82, 83, 80, 77, 78, 79, 86, 87, 88, 85, 113, 114, 115, 116, 98, 99, 100, 97, 117, 118, 119, 120, 108, 105, 106, 107, 104, 101, 102, 103, 110, 111, 112, 109], [41, 42, 43, 44, 26, 27, 28, 25, 45, 46, 47, 48, 36, 33, 34, 35, 32, 29, 30, 31, 38, 39, 40, 37, 17, 18, 19, 20, 2, 3, 4, 1, 21, 22, 23, 24, 12, 9, 10, 11, 8, 5, 6, 7, 14, 15, 16, 13, 65, 66, 67, 68, 50, 51, 52, 49, 69, 70, 71, 72, 60, 57, 58, 59, 56, 53, 54, 55, 62, 63, 64, 61, 113, 114, 115, 116, 98, 99, 100, 97, 117, 118, 119, 120, 108, 105, 106, 107, 104, 101, 102, 103, 110, 111, 112, 109, 89, 90, 91, 92, 74, 75, 76, 73, 93, 94, 95, 96, 84, 81, 82, 83, 80, 77, 78, 79, 86, 87, 88, 85]],
                            "concs": [60,60],
                            "size": 144,
                            "result": 670936866946976150457457926209536},
                           {"fname": 'fc1',
                              "gens": [[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,1]],
                            "concs": [10,10,20],
                            "size": 40,
                            "result": 636699333130704},
                           {"fname": 'fc2',
                            "gens": [[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,1],[1,40,39,38,37,36,35,34,33,32,31,30,29,28,27,26,25,24,23,22,21,20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2]],
                            "concs": [10,10,20],
                            "size": 80,
                            "result": 318349689844608},
                           {"fname": 'fc3',
                            "gens": [[10, 11, 12, 9, 14, 15, 16, 13, 18, 19, 20, 17, 22, 23, 24, 21, 2, 3, 4, 1, 6, 7, 8, 5, 26, 27, 28, 25, 30, 31, 32, 29, 34, 35, 36, 33, 38, 39, 40, 37], [10, 11, 12, 9, 14, 15, 16, 13, 2, 3, 4, 1, 6, 7, 8, 5, 18, 19, 20, 17, 22, 23, 24, 21, 34, 35, 36, 33, 38, 39, 40, 37, 26, 27, 28, 25, 30, 31, 32, 29], [9, 12, 16, 13, 10, 11, 15, 14, 17, 20, 24, 21, 18, 19, 23, 22, 1, 4, 8, 5, 2, 3, 7, 6, 25, 28, 32, 29, 26, 27, 31, 30, 33, 36, 40, 37, 34, 35, 39, 38], [9, 12, 16, 13, 10, 11, 15, 14, 1, 4, 8, 5, 2, 3, 7, 6, 17, 20, 24, 21, 18, 19, 23, 22, 33, 36, 40, 37, 34, 35, 39, 38, 25, 28, 32, 29, 26, 27, 31, 30]],
                            "concs": [10,10,20],
                            "size": 144,
                            "result": 176860970824040},
                           {"fname": 'fc4',
                            "gens": [[2, 3, 4, 5, 6, 7, 8, 1, 10, 11, 12, 13, 14, 15, 16, 9, 18, 19, 20, 21, 22, 23, 24, 17, 26, 27, 28, 29, 30, 31, 32, 25, 34, 35, 36, 37, 38, 39, 40, 33], [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 1, 2, 3, 4, 5, 6, 7, 8], [1, 8, 7, 6, 5, 4, 3, 2, 9, 16, 15, 14, 13, 12, 11, 10, 17, 24, 23, 22, 21, 20, 19, 18, 25, 32, 31, 30, 29, 28, 27, 26, 33, 40, 39, 38, 37, 36, 35, 34], [1,2,3,4,5,6,7,8,33,34,35,36,37,38,39,40,25,26,27,28,29,30,31,32,17,18,19,20,21,22,23,24,9,10,11,12,13,14,15,16]],
                            "concs": [10,10,20],
                            "size": 160,
                            "result": 159174945930814},
                           {"fname": 'fc5',
                            "gens": [[2,3,4,1,6,7,8,5, 10,11,12,9, 14,15,16,13, 18,19,20,17, 22,23,24,21, 26,27,28,25, 30,31,32,29, 34,35,36,33, 38,39,40,37 ], [2,1,3,4, 6,5,7,8, 10,9,11,12, 14,13,15,16, 18,17,19,20, 22,21,23,24, 26,25,27,28, 30,29,31,32, 34,33,35,36, 38,37,39,40 ], [ 5,6,7,8, 9,10,11,12, 13,14,15,16, 17,18,19,20, 21,22,23,24, 25,26,27,28, 29,30,31,32, 33,34,35,36, 37,38,39,40, 1,2,3,4],  [ 1,2,3,4, 37,38,39,40, 33,34,35,36, 29,30,31,32, 25,26,27,28, 21,22,23,24, 17,18,19,20, 13,14,15,16, 9,10,11,12, 5,6,7,8 ]],
                            "concs": [10,10,20],
                            "size": 480,
                            "result": 53066285026564}]

        self.groups = {}
        self.products = {}
        self.total = 0

    def gen_group(self, index):
        """Tests the group generation from the generators. The method outputs generators.in.* 
        and group.out.*"""
        if index not in self.groups:
            single = self.generators[index]
            gens = single["gens"]
            save([gens],["generators"],self.folder,case=single["fname"])
            grpops = group(gens)
            self.assertEqual(len(grpops), single["size"], "{} != {}".format(len(grpops), 
                                                                            single["size"]))
            #Save the input/output files for the fortran unit tests.
            save([grpops],["group"],self.folder,infile=False,case=single["fname"])
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
        save([marray], ["multinomial"], self.folder, case=self.generators[index]["fname"])

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
        save([possibles],["possibles"],self.folder,case=self.generators[index]["fname"])

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
            save([seq],["seq"],self.folder,case=self.generators[index]["fname"])
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
                     self.folder,case=self.generators[index]["fname"],infile=False)

            #Save the input/output files for the fortran unit tests.
            save(mnseq,["mnseq.in.{}".format(i+1) for i in range(len(mnseq))],self.folder, 
                 infile=None, case=self.generators[index]["fname"])
            save([len(mnseq)], ["mnseq.len"], self.folder, case=self.generators[index]["fname"])

            #add this to our dictionary of products if it's not there yet
            if index not in self.products:
                self.gen_product(index)
            p = self.products[index]
            sumseq = int(p._sum_sequences(mnseq))
            
            #Save the input/output files for the fortran unit tests.
            save([sumseq], ["sum_seq"], self.folder, infile=False, case=self.generators[index]["fname"])

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
        save([coeffs],["coeffs"], self.folder, infile=False, case=self.generators[index]["fname"])

    def gen_test(self, index):
        """ Tests the polya code for each of the given test cases. The method outputs concs.in.*
        and results.out.*"""
        single = self.generators[index]
        self.gen_group(index)
        grpops = self.groups[index]

        #Save the input/output files for the fortran unit tests.
        save([single["concs"]],["concs"],self.folder,case=single["fname"])
        save([single["result"]],["result"],self.folder,case=single["fname"],infile=False)

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
        """Builds an array ofproducts for testing."""
        for i in range(len(self.generators)):
            self.gen_product(i)

    def test_coefficients(self):
        """Runs a test of the code and checks the results then prints the accurracy, as 
        percent correct, to the screen"""
        for i in range(len(self.generators)):
           self.gen_test(i)
        print("TOTAL: %.06f" % self.total)

