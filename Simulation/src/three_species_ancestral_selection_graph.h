#include "terminal_ancestral_selection_graph.h"
#include "internal_ancestral_selection_graph.h"

class ThreeSpeciesAncestralSelectionGraph: public BaseAncestralSelectionGraph
{
    private:
        TerminalAncestralSelectionGraph terminal_1;
        TerminalAncestralSelectionGraph terminal_2;
        TerminalAncestralSelectionGraph terminal_3;
        InternalAncestralSelectionGraph internal_23;
        InternalAncestralSelectionGraph internal_1_23;

        double terminal_population_size_1;
        double terminal_population_size_2;
        double terminal_population_size_3;
        double internal_population_size_23;
        double internal_population_size_1_23;

        double divergence_1;
        double divergence_2;



    public:
        ThreeSpeciesAncestralSelectionGraph(double mutation_rate, double selection_coefficient, int terminal_population_size_1, int terminal_population_size_2, int terminal_population_size_3,
            int internal_population_size_23, int internal_population_size_1_23, double divergence_1, double divergence_2);
        void getTerminals(vector<TerminalAncestralSelectionGraph*>& terminals);
};