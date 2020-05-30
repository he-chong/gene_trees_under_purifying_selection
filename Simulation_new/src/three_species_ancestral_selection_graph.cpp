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
#include "three_species_ancestral_selection_graph.h"

using namespace std;

class unexpected_ua : public exception  
{  
public:  
    const char* what()  
    {  
        return "unexpected UA";  
    }  
}; 

ThreeSpeciesAncestralSelectionGraph::ThreeSpeciesAncestralSelectionGraph(double mut_rate, double sel_coef, int pop_size_1, int pop_size_2, int pop_size_3,
        int pop_size_23, int pop_size_1_23, double t_1, double dt): AncestralSelectionGraphComposer(),
            global_mutation_rate(mut_rate),
            global_selection_coefficient(sel_coef),
            population_size_1(pop_size_1),
            population_size_2(pop_size_2),
            population_size_3(pop_size_3),
            population_size_23(pop_size_23),
            population_size_1_23(pop_size_1_23),
            tau_1(t_1),
            delta_tau(dt),
            terminal_1("1", mut_rate, sel_coef, pop_size_1, 0, t_1+dt, 1, composer_random_engine),
            terminal_2("2", mut_rate, sel_coef, pop_size_2, 0, t_1, 1, composer_random_engine),
            terminal_3("3", mut_rate, sel_coef, pop_size_3, 0, t_1, 1, composer_random_engine),
            internal_23("(2,3)", mut_rate, sel_coef, pop_size_23, t_1, t_1+dt, terminal_2, terminal_3, composer_random_engine),
            internal_1_23("(1,(2,3))", mut_rate, sel_coef, pop_size_1_23, t_1, NONE, terminal_1, internal_23, composer_random_engine)
{
    
}


void ThreeSpeciesAncestralSelectionGraph::__setComponents() {
    components.push_back(&internal_1_23);
    components.push_back(&internal_23);
    components.push_back(&terminal_3);
    components.push_back(&terminal_2);
    components.push_back(&terminal_1);
}

void ThreeSpeciesAncestralSelectionGraph::__setUltimateAncestor() {
    if(internal_1_23.getOutNodes().size() != 1)
    {
        throw unexpected_ua();
    }
    ultimate_ancestor = internal_1_23.getOutNodes()[0];
}


void ThreeSpeciesAncestralSelectionGraph::getTerminals(vector<AncestralSelectionGraphComponent*>& terminals) {
    terminals.push_back(&terminal_1);
    terminals.push_back(&terminal_2);
    terminals.push_back(&terminal_3);
}

void ThreeSpeciesAncestralSelectionGraph::getExtantLineages(vector<Node*>& extantLineages) {
    extantLineages.push_back(terminal_1.getInNodes()[0]);
    extantLineages.push_back(terminal_2.getInNodes()[0]);
    extantLineages.push_back(terminal_3.getInNodes()[0]);
}


