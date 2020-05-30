#include <iostream>
#include <algorithm>
#include <chrono>
#include <random>
#include "terminal_ancestral_selection_graph.h"

/* 
 * @brief The constructor of the population of the Terminal branch. 
 *        The initialization of the population of the terminal branch
 *        needs to assign a number of extant lineages in the parameter
 *        sample_num. time_upbound=t means that the backward-in-time 
 *        process stops at the time t.
 *
 */
TerminalAncestralSelectionGraph::TerminalAncestralSelectionGraph (string name, int population_size, double mutation_rate, double selection_coefficient, double time_upbound, 
            int sample_num):BaseAncestralSelectionGraph(name, population_size, mutation_rate, selection_coefficient, time_upbound)
{
    this->sample_num = sample_num;
    for (int i=0; i<sample_num; i++)
    {
        string used_name = (sample_num == 1) ? name : (name != "" ? (name + "_" + to_string(i+1)) : to_string(i+1));
        Node* leaf = new Node(used_name);
        Branch* branch = new Branch;
        leaf->setCBranch(branch);
        branch->setDescendant(leaf);
        all_nodes.push_back(leaf);
        all_branches.push_back(branch);
        temp_branches.push_back(branch);
    }

    _buildAncestralSelectionGraph();
    // cout << temp_branches.size() << endl;
}


void TerminalAncestralSelectionGraph::getLeaves(vector<Node*>& leaves) {
    for(int i=0; i<sample_num; i++) {
        leaves.push_back(all_nodes[i]);
    }
}