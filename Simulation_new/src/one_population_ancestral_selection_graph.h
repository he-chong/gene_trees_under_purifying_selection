#ifndef ONE_POPULATION_ANCESTRAL_SELECTION_GRAPH_H
#define ONE_POPULATION_ANCESTRAL_SELECTION_GRAPH_H
#include "ancestral_selection_graph_composer.h"

using namespace std;


class OnePopulationAncestralSelectionGraph: public AncestralSelectionGraphComposer
{
    private:
        AncestralSelectionGraphComponent population;

    protected:
    	void __setUltimateAncestor();
    	void __setComponents();

    public:
        OnePopulationAncestralSelectionGraph(double mut_rate, double sel_coef, int pop_size, int sample_num=3);
        void getExtantLineages(vector<Node*>& extantLineages);
};
#endif