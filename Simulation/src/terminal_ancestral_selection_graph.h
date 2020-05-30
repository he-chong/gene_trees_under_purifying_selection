#ifndef TERMINAL_ANCESTRAL_SELECTION_GRAPH_H
#define TERMINAL_ANCESTRAL_SELECTION_GRAPH_H
#include "base_ancestral_selection_graph.h"

using namespace std;


class TerminalAncestralSelectionGraph: public BaseAncestralSelectionGraph
{
    public:
        TerminalAncestralSelectionGraph(string name, int population_size, double mutation_rate, double selection_coefficient, double time_upbound, int sample_num=1);
        void getLeaves(vector<Node*>& leaves);


};
#endif