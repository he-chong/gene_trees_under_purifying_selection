#include <iostream>
#include <algorithm>
#include <random>
#include "internal_ancestral_selection_graph.h"

/* 
 * @brief The constructor of the population of the internal branch. 
 *        The initialization of the population of the internal branch
 *        needs the ancestral lineages in the descendant populations.
 *        time_upbound=None means that the backward-in-time process 
 *        stops when there exist only one lineage, i.e. the ultimate 
 *        ancestor, UA. time_upbound=t means that the backward-
 *        in-time process stops at the time t.
 *
 */
InternalAncestralSelectionGraph::InternalAncestralSelectionGraph (string name, int population_size, double mutation_rate, double selection_coefficient, double time_upbound, 
            BaseAncestralSelectionGraph& asg_1, BaseAncestralSelectionGraph& asg_2):
            BaseAncestralSelectionGraph(name, population_size, mutation_rate, selection_coefficient, time_upbound)
{
    this->sample_num = asg_1.getTempBranches().size() + asg_2.getTempBranches().size();
    temp_branches.insert(temp_branches.end(), asg_1.getTempBranches().begin(), asg_1.getTempBranches().end());
    temp_branches.insert(temp_branches.end(), asg_2.getTempBranches().begin(), asg_2.getTempBranches().end());

    _buildAncestralSelectionGraph();
    
}


