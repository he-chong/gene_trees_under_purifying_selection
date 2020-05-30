#include <iostream>
#include <algorithm>
#include <random>
#include "ancestral_selection_graph.h"

/* 
 * @brief The initialization of The Ancestral Selection Graph within
 *        a population needs to assign a number of extant lineages 
 *        in the parameter sample_num. time_upbound=None means that
 *        its backward-in-time process stops when there exist only
 *        one lineage, i.e. the ultimate ancestor, UA. 
 */
AncestralSelectionGraph::AncestralSelectionGraph (int population_size, double mutation_rate, double selection_coefficient, 
            int sample_num):BaseAncestralSelectionGraph(UNNAMED, population_size, mutation_rate, selection_coefficient, NONE)
{
    this->sample_num = sample_num;

    for (int i=0; i<sample_num; i++)
    {
        string used_name = (sample_num == 1) ? name : (name != UNNAMED ? (name + "_" + to_string(i+1)) : to_string(i+1));
        Node* leaf = new Node(used_name);
        Branch* branch = new Branch;
        leaf->setCBranch(branch);
        branch->setDescendant(leaf);
        all_nodes.push_back(leaf);
        all_branches.push_back(branch);
        temp_branches.push_back(branch);
    }
    
    _buildAncestralSelectionGraph();
}

void AncestralSelectionGraph::getLeaves(vector<Node*>& leaves) {
    for(int i=0; i<sample_num; i++) {
        leaves.push_back(all_nodes[i]);
    }
}
