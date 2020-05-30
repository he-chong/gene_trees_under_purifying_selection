import re

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "for_ancestral_selection_graph.h":
    cdef cppclass Node:
        string getName()
        int getGenotype()

cdef extern from "ancestral_selection_graph.h":
    cdef cppclass AncestralSelectionGraph:
        AncestralSelectionGraph(int population_size, double mutation_rate, double selection_coefficient, int sample_num)
        void getLeaves(vector[Node*]& leaves)
        string pick_embedded_genealogy()


cdef class CythonAncestralSelectionGraph:
    cdef AncestralSelectionGraph* asg
    def __cinit__(self, int population_size, double mutation_rate, double selection_coefficient, int sample_num):
        self.asg = new AncestralSelectionGraph(population_size, mutation_rate, selection_coefficient, sample_num)
    def __dealloc__(self):
        del self.asg

    def pick_embedded_genealogy(self):
        cdef string newick_tree = self.asg.pick_embedded_genealogy()
        return newick_tree

    @property
    def genotypes(self):
        cdef vector[Node*] leaves
        self.asg.getLeaves(leaves)
        genotypes = []
        for i from 0 <= i < leaves.size():         
            genotype = leaves[i].getGenotype()
            genotypes.append(genotype)
        leaves.clear()
        return genotypes
