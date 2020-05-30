/*****************************************************************************
*                                                                            *
*  This file is contains codes that defines the structure of the species     *
*  tree, and the parameters of each population in the tree.                  *
*                                                                            *
*  @author   何充 He Chong                                                   *
*  @email    biohe@foxmail.com                                               *
*                                                                            *
*****************************************************************************/

#include <iostream>
#include <deque>
#include "base_ancestral_selection_graph.h"
#include "three_species_ancestral_selection_graph.h"

using namespace std;


ThreeSpeciesAncestralSelectionGraph::ThreeSpeciesAncestralSelectionGraph(double mutation_rate, double selection_coefficient, int terminal_population_size_1, int terminal_population_size_2, int terminal_population_size_3,
        int internal_population_size_23, int internal_population_size_1_23, double divergence_1, double divergence_2): BaseAncestralSelectionGraph(),
            terminal_1("1", terminal_population_size_1, mutation_rate, selection_coefficient, divergence_1 + divergence_2, sample_num=1),
            terminal_2("2", terminal_population_size_2, mutation_rate, selection_coefficient, divergence_1, sample_num=1),
            terminal_3("3", terminal_population_size_3, mutation_rate, selection_coefficient, divergence_1, sample_num=1),
            internal_23("(2,3)", internal_population_size_23, mutation_rate, selection_coefficient, divergence_2, terminal_2, terminal_3),
            internal_1_23("(1,(2,3))", internal_population_size_1_23, mutation_rate, selection_coefficient, NONE, terminal_1, internal_23)
{
    this->mutation_rate = mutation_rate;
    this->selection_coefficient = selection_coefficient;
    this->terminal_population_size_1 = terminal_population_size_1;
    this->terminal_population_size_2 = terminal_population_size_2;
    this->terminal_population_size_3 = terminal_population_size_3;
    this->internal_population_size_23 = internal_population_size_23;
    this->internal_population_size_1_23 = internal_population_size_1_23;
    this->divergence_1 = divergence_1;
    this->divergence_2 = divergence_2; 

    ultimate_ancestor = internal_1_23.getUltimateAncestor();
}


void ThreeSpeciesAncestralSelectionGraph::getTerminals(vector<TerminalAncestralSelectionGraph*>& terminals) {
    terminals.push_back(&terminal_1);
    terminals.push_back(&terminal_2);
    terminals.push_back(&terminal_3);
}


