import re

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "for_ancestral_selection_graph.h":
    cdef cppclass Node:
        string getName()
        int getGenotype()

cdef extern from "terminal_ancestral_selection_graph.h":
    cdef cppclass TerminalAncestralSelectionGraph:
        void getLeaves(vector[Node*]& leaves)

cdef extern from "three_species_ancestral_selection_graph.h":
    cdef cppclass ThreeSpeciesAncestralSelectionGraph:
        ThreeSpeciesAncestralSelectionGraph(double mutation_rate, double selection_coefficient, int terminal_population_size_1, int terminal_population_size_2, int terminal_population_size_3,\
            int internal_population_size_23, int internal_population_size_1_23, double divergence_1, double divergence_2)
        string pick_embedded_genealogy()
        void getTerminals(vector[TerminalAncestralSelectionGraph*]& terminals);


cdef class CythonThreeSpeciesAncestralSelectionGraph:
    cdef ThreeSpeciesAncestralSelectionGraph* three_species_asg
    def __cinit__(self, double mutation_rate, double selection_coefficient, int terminal_population_size_1, int terminal_population_size_2, int terminal_population_size_3,\
            int internal_population_size_23, int internal_population_size_1_23, double divergence_1, double divergence_2):
        self.three_species_asg = new ThreeSpeciesAncestralSelectionGraph(mutation_rate, selection_coefficient, terminal_population_size_1, terminal_population_size_2, terminal_population_size_3, \
            internal_population_size_23, internal_population_size_1_23, divergence_1, divergence_2)
    def __dealloc__(self):
        del self.three_species_asg

    def pick_embedded_genealogy(self):
        cdef string newick_tree = self.three_species_asg.pick_embedded_genealogy()
        return newick_tree

    @property
    def genotypes(self):
        cdef vector[TerminalAncestralSelectionGraph*] terminals_vector
        cdef vector[Node*] leaves
        self.three_species_asg.getTerminals(terminals_vector)
        cdef TerminalAncestralSelectionGraph *terminal
        genotypes = []
        for i from 0 <= i < terminals_vector.size():
            terminal = terminals_vector[i]
            terminal.getLeaves(leaves)
            # print leaves[0].getName(), leaves[0].getGenotype()
            genotype = leaves[0].getGenotype()
            genotypes.append(genotype)
            leaves.clear()
        terminals_vector.clear()
        return genotypes
