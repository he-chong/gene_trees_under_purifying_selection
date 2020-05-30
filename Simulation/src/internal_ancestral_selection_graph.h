#ifndef INTERNAL_ANCESTRAL_SELECTION_GRAPH_H
#define INTERNAL_ANCESTRAL_SELECTION_GRAPH_H
#include "base_ancestral_selection_graph.h"

using namespace std;


class InternalAncestralSelectionGraph: public BaseAncestralSelectionGraph
{
    public:
        InternalAncestralSelectionGraph(string name, int population_size, double mutation_rate, double selection_coefficient, double time_upbound,
            BaseAncestralSelectionGraph& asg_1, BaseAncestralSelectionGraph& asg_2);
};
#endif