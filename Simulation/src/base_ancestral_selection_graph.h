#ifndef BASE_ANCESTRAL_SELECTION_GRAPH_H
#define BASE_ANCESTRAL_SELECTION_GRAPH_H
#define UNKNOWN -1
#define NONE 0
#define UNNAMED ""
#define ADVANTAGEOUS 1
#define DISADVANTAGEOUS 2
#include <vector>
#include <deque>
#include <string>
#include <random>
#include <chrono>
#include "for_ancestral_selection_graph.h"

using namespace std;


class BaseAncestralSelectionGraph
{
    protected:
        string name;
        int population_size;
        double mutation_rate;
        double selection_coefficient;
        double theta;
        double sigma;
        double time_upbound;
        int sample_num;
        vector<Node*> all_nodes;
        vector<Branch*> all_branches;
        Node* ultimate_ancestor;
        deque<Branch*> temp_branches;
        mt19937 random_engine;

        void _coalesce();
        void _branch();
        void _mutate();
        string _event();

        int _assignGenotype(Node* node);
        void _connect(Node* node);
        string _output(Node* node);
        void _trim(Node* node);
        void _buildAncestralSelectionGraph();

    public:
        BaseAncestralSelectionGraph(string name, int population_size, double mutation_rate, double selection_coefficient, double time_upbound);
        BaseAncestralSelectionGraph();
        ~BaseAncestralSelectionGraph();
        string pick_embedded_genealogy();
        
        Node* getUltimateAncestor() {
            return ultimate_ancestor;
        }

        deque<Branch*>& getTempBranches() {
            return temp_branches;
        }

};
#endif