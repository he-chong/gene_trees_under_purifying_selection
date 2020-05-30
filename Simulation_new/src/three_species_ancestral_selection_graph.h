#include "ancestral_selection_graph_composer.h"

class ThreeSpeciesAncestralSelectionGraph: public AncestralSelectionGraphComposer
{
    private:
        AncestralSelectionGraphComponent terminal_1;
        AncestralSelectionGraphComponent terminal_2;
        AncestralSelectionGraphComponent terminal_3;
        AncestralSelectionGraphComponent internal_23;
        AncestralSelectionGraphComponent internal_1_23;

        double global_mutation_rate;
        double global_selection_coefficient;

        double population_size_1;
        double population_size_2;
        double population_size_3;
        double population_size_23;
        double population_size_1_23;

        double tau_1;
        double delta_tau;

    protected:
        void __setUltimateAncestor();
        void __setComponents();

    public:
        ThreeSpeciesAncestralSelectionGraph(double mut_rate, double sel_coef, int terminal_pop_size_1, int terminal_pop_size_2, int terminal_pop_size_3,
        int internal_pop_size_23, int internal_pop_size_1_23, double t_1, double dt);
        void getTerminals(vector<AncestralSelectionGraphComponent*>& terminals);
        void getExtantLineages(vector<Node*>& extantLineages);
};