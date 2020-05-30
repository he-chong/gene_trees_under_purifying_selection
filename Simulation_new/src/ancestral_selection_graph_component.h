#ifndef ANCESTRAL_SELECTION_GRAPH_COMPONENT_H
#define ANCESTRAL_SELECTION_GRAPH_COMPONENT_H
#define UNKNOWN -1
#define NONE -1
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


class AncestralSelectionGraphComponent
{
    protected:
        string name;
        int population_size;
        double mutation_rate;
        double selection_coefficient;
        double theta;
        double sigma;
        double height_lower;
        double height_upper;
        int sample_num;
        vector<Node*> all_nodes;
        vector<Branch*> all_branches;
        vector<Node*> in_nodes;
        vector<Node*> out_nodes;
        deque<Branch*> temp_branches;
        mt19937& random_engine;

        void _coalesce(double height);
        void _branch(double height);
        void _mutate();
        void _aftercare(double height);
        string _event(double height);
        
        void _buildAncestralSelectionGraph();

    public:
        AncestralSelectionGraphComponent(string n, double mut_rate, double sel_coef, int pop_size, double time_up, double h_base, int samp_num, mt19937& rand_eng);
        AncestralSelectionGraphComponent(string n, double mut_rate, double sel_coef, int pop_size, double time_up, double h_base,
        AncestralSelectionGraphComponent& component_1, AncestralSelectionGraphComponent& component_2, mt19937& rand_eng);
        ~AncestralSelectionGraphComponent();

        string getName() {
            return name;
        }

        int getPopulationSize() {
            return population_size;
        }

        double getMutationRate() {
            return mutation_rate;
        }

        double getSelectionCoefficient() {
            return selection_coefficient;
        }

        int getSampleNum() {
            return sample_num;
        }

        vector<Node*>& getInNodes() {
            return in_nodes;
        }

        vector<Node*>& getOutNodes() {
            return out_nodes;
        }

        vector<Node*>& getAllNodes() {
            return all_nodes;
        }

        deque<Branch*>& getTempBranches() {
            return temp_branches;
        }

};
#endif