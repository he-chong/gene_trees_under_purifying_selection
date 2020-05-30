# distutils: language = c++

import re
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "../src/for_ancestral_selection_graph.h":
    cdef cppclass Node:
        string getName()
        int getGenotype()

cdef extern from "../src/one_population_ancestral_selection_graph.h":
    cdef cppclass OnePopulationAncestralSelectionGraph:
        OnePopulationAncestralSelectionGraph(double mut_rate, double sel_coef, int pop_size, int samp_num)
        string pick_embedded_genealogy()
        void getExtantLineages(vector[Node*]& leaves)


cdef class CythonOnePopulationAncestralSelectionGraph:
    cdef OnePopulationAncestralSelectionGraph* asg
    def __cinit__(self, double mut_rate, double sel_coef, int pop_size, int samp_num):
        self.asg = new OnePopulationAncestralSelectionGraph(mut_rate, sel_coef, pop_size, samp_num)
    def __dealloc__(self):
        del self.asg

    def pick_embedded_genealogy(self):
        cdef string newick_tree = self.asg.pick_embedded_genealogy()
        return newick_tree

    @property
    def genotypes(self):
        cdef vector[Node*] leaves
        self.asg.getExtantLineages(leaves)
        genotypes = []
        for i from 0 <= i < leaves.size():         
            genotype = leaves[i].getGenotype()
            genotypes.append(genotype)
        leaves.clear()
        return genotypes
