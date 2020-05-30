#include <iostream>
#include <ctime>
#include <exception>
#include <vector>
#include "ancestral_selection_graph.h"
using namespace std;

int main()
{
    double mutation_rate;
    double selection_coefficient;
    double population_size;
    int sample_num;

    cout << "# This is for demostrating the result of Ancestral Selection Graph model simulating the genealogy within a population" << endl;
    cout << endl;
    cout << "# Please type into mutation_rate, selection_coefficient, population_size_1, sample_num" << endl;
    cout << endl;
    cout << "# For example: 3e-8 -7.5e-6 2e5 3" << endl;
    cin >> mutation_rate >> selection_coefficient >> population_size >> sample_num;

    time_t begin = time(0);
    char buffer[9] = {0};
    strftime(buffer, 9, "%H:%M:%S", localtime(&begin));
    cout << "Begin: " << string(buffer) << endl;
    cout << endl;
    for (int i=0; i<10; i++)
    {
        cout << "# Repeat: " << i+1 << endl;
        cout << "Building ancestral selection graph..." << endl;
        AncestralSelectionGraph asg(population_size, mutation_rate, selection_coefficient, sample_num);
        cout << "Picking embedded genealogy..." << endl;
        string newick_tree = asg.pick_embedded_genealogy();
        cout << "Sampled genealogy: " << newick_tree << endl;
        vector<Node*> leaves;
        asg.getLeaves(leaves);
        for(auto& leaf : leaves) {
            cout << "Genotype_" << leaf->getName() << ": ";
            cout << leaf->getGenotype() << endl;
        }
        cout << endl;
        leaves.clear();
    }
    time_t end = time(0);
    strftime(buffer, 9, "%H:%M:%S", localtime(&end));
    cout << "End: " << string(buffer) << endl;

    return 0;
}