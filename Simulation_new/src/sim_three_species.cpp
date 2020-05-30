#include <iostream>
#include <iterator>
#include <ctime>
#include <exception>
#include "three_species_ancestral_selection_graph.h"
using namespace std;

int main()
{
    double mutation_rate;
    double selection_coefficient;
    double terminal_population_size_1;
    double terminal_population_size_2;
    double terminal_population_size_3;
    double internal_population_size_23;
    double internal_population_size_1_23;
    double divergence_1;
    double divergence_2;

    cout << "# This is for demostrating the result of three-species Ancestral Selection Graph model. The speices tree is given as (1, (2, 3))" << endl;
    cout << endl;
    cout << "# Please type into mutation_rate, selection_coefficient, terminal_population_size_1, terminal_population_size_2, terminal_population_size_3, internal_population_size_23, internal_population_size_1_23, divergence_1, divergence_2:" << endl;
    cout << endl;
    cout << "# For example: 3e-8 -7.5e-6 1e6 2e5 4e4 2e5 2e5 2e6 1.4e4" << endl;
    cin >> mutation_rate >> selection_coefficient >> terminal_population_size_1 >> terminal_population_size_2 >> terminal_population_size_3 >> internal_population_size_23 >> internal_population_size_1_23 >> divergence_1 >> divergence_2;

    time_t begin = time(0);
    char buffer[9] = {0};
    strftime(buffer, 9, "%H:%M:%S", localtime(&begin));
    cout << "Begin: " << string(buffer) << endl;
    cout << endl;
    for (int i=0; i<10; i++)
    {
        cout << "# Repeat: " << i+1 << endl;
        cout << "Building ancestral selection graph..." << endl;
        ThreeSpeciesAncestralSelectionGraph inter_asg(mutation_rate, selection_coefficient, terminal_population_size_1, terminal_population_size_2, terminal_population_size_3,
            internal_population_size_23, internal_population_size_1_23, divergence_1, divergence_2);
        cout << "Picking embedded genealogy..." << endl;
        string newick_tree = inter_asg.pick_embedded_genealogy();
        // string newick_tree = "#";
        cout << "Sampled genealogy: " << newick_tree << endl;
        vector<Node*> extantLineages;
        inter_asg.getExtantLineages(extantLineages);
        for(auto& leaf : extantLineages) {
            cout << "Genotype_" << leaf->getName() << ": ";
            cout << leaf->getGenotype() << endl;
        }
        cout << endl;
    }
    time_t end = time(0);
    strftime(buffer, 9, "%H:%M:%S", localtime(&end));
    cout << "End: " << string(buffer) << endl;
    return 0;
}