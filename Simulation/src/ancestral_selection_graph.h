#ifndef ANCESTRAL_SELECTION_GRAPH_H
#define ANCESTRAL_SELECTION_GRAPH_H
#include "base_ancestral_selection_graph.h"

using namespace std;


class AncestralSelectionGraph: public BaseAncestralSelectionGraph
{
    public:
        AncestralSelectionGraph(int population_size, double mutation_rate, double selection_coefficient, int sample_num=3);
        void getLeaves(vector<Node*>& leaves);
};
#endif