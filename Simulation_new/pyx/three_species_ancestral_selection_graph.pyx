import re

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "../src/for_ancestral_selection_graph.h":
    cdef cppclass Node:
        string getName()
        int getGenotype()


cdef extern from "../src/three_species_ancestral_selection_graph.h":
    cdef cppclass ThreeSpeciesAncestralSelectionGraph:
        ThreeSpeciesAncestralSelectionGraph(double mut_rate, double sel_coef, int pop_size_1, int pop_size_2, int pop_size_3,\
            int pop_size_23, int pop_size_1_23, double t_1, double dt)
        string pick_embedded_genealogy()
        void getExtantLineages(vector[Node*]& leaves)


cdef class CythonThreeSpeciesAncestralSelectionGraph:
    cdef ThreeSpeciesAncestralSelectionGraph* three_species_asg
    def __cinit__(self, double mut_rate, double sel_coef, int pop_size_1, int pop_size_2, int pop_size_3,\
            int pop_size_23, int pop_size_1_23, double t_1, double dt):
        self.three_species_asg = new ThreeSpeciesAncestralSelectionGraph(mut_rate, sel_coef, pop_size_1, pop_size_2, pop_size_3, \
            pop_size_23, pop_size_1_23, t_1, dt)

    def __dealloc__(self):
        del self.three_species_asg

    def pick_embedded_genealogy(self):
        cdef string newick_tree = self.three_species_asg.pick_embedded_genealogy()
        return newick_tree

    @property
    def genotypes(self):
        cdef vector[Node*] leaves
        self.three_species_asg.getExtantLineages(leaves)
        genotypes = []
        for i from 0 <= i < leaves.size():
            # print leaves[0].getName(), leaves[0].getGenotype()
            genotype = leaves[i].getGenotype()
            genotypes.append(genotype)
        leaves.clear()
        return genotypes
