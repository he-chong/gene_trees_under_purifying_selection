#include <iostream>
#include <algorithm>
#include <random>
#include "one_population_ancestral_selection_graph.h"

/* 
 * @brief The initialization of The Ancestral Selection Graph within
 *        a population needs to assign a number of extant lineages 
 *        in the parameter sample_num. The parameter "height_upper" is set to be "NONE" 
 *        which means that the backward-in-time process stops when there exist only
 *        one lineage, i.e. the ultimate ancestor, UA. 
 */
OnePopulationAncestralSelectionGraph::OnePopulationAncestralSelectionGraph (double mut_rate, double sel_coef, int pop_size, 
            int samp_num):AncestralSelectionGraphComposer(),
    population("#", mut_rate, sel_coef, pop_size, 0, NONE, samp_num, composer_random_engine)
{

}


void OnePopulationAncestralSelectionGraph::__setComponents() {
    components.push_back(&population);
}


void OnePopulationAncestralSelectionGraph::__setUltimateAncestor() {
    ultimate_ancestor = population.getOutNodes()[0];
}


void OnePopulationAncestralSelectionGraph::getExtantLineages(vector<Node*>& extantLineages) {
    for(int i=0; i<population.getSampleNum(); i++) {
        extantLineages.push_back(population.getInNodes()[i]);
    }
}
