/*********************************************************************************
*                                                                                *
*  This file is contains codes that implements the backward-in-time process of   *
*  Ancestral Selection Graph                                                     *
*                                                                                *
*  @author   何充 He Chong                                                       *
*  @email    biohe@foxmail.com                                                   *
*                                                                                *
**********************************************************************************/

#include <iostream>
#include <algorithm>
#include <sstream>
#include <exception>
#include <iomanip>
#include "ancestral_selection_graph_component.h"

using namespace std;


/* 
 * @brief The initialization of the Ancestral Selection Graph in whatever situation
 *        needs to assign population size, mutation rate, selection_coefficient.
 *        height_upper=None means that the backward-in-time process 
 *        stops when there exist only one lineage, i.e. the ultimate 
 *        ancestor, UA. height_upper=t means that the backward-
 *        in-time process stops at the time t.
 *
 */
AncestralSelectionGraphComponent::AncestralSelectionGraphComponent (string n, double mut_rate, double sel_coef, int pop_size, double h_low, double h_up,
    int samp_num, mt19937& rand_eng):
    name(n), population_size(pop_size), mutation_rate(mut_rate), selection_coefficient(sel_coef), height_lower(h_low), height_upper(h_up), sample_num(samp_num),
    random_engine(rand_eng)
{
    double fitness_advantage = selection_coefficient < 0 ? fabs(selection_coefficient)/(1+selection_coefficient) : selection_coefficient;
    this->theta = 2*mutation_rate*population_size; // theta is the expected number of mutation events in a coalescent time unit.     
    this->sigma = 2*(fitness_advantage)*population_size; // sigma is the expected number of branching events in a coalescent time unit. 
    for (int i=0; i<sample_num; i++)
    {
        string used_name = (sample_num == 1) ? name : (name != UNNAMED ? (name + "_" + to_string(i+1)) : to_string(i+1));
        Node* leaf = new Node(used_name, height_lower);
        // cout << leaf->getGenotype() << endl;
        Branch* branch = new Branch;
        leaf->setCBranch(branch);
        branch->setDescendant(leaf);
        all_nodes.push_back(leaf);
        all_branches.push_back(branch);
        in_nodes.push_back(leaf);
        temp_branches.push_back(branch);
    }
    _buildAncestralSelectionGraph();   
}


AncestralSelectionGraphComponent::AncestralSelectionGraphComponent(string n, double mut_rate, double sel_coef, int pop_size, double h_low, double h_up,
    AncestralSelectionGraphComponent& component_1, AncestralSelectionGraphComponent& component_2, mt19937& rand_eng):
    name(n), population_size(pop_size), mutation_rate(mut_rate), selection_coefficient(sel_coef), height_lower(h_low), height_upper(h_up),
    random_engine(rand_eng)
{
    double fitness_advantage = selection_coefficient < 0 ? fabs(selection_coefficient)/(1+selection_coefficient) : selection_coefficient;
    theta = 2*mutation_rate*population_size; // theta is the expected number of mutation events in a coalescent time unit.
    sigma = 2*(fitness_advantage)*population_size; // sigma is the expected number of branching events in a coalescent time unit. 

    sample_num = component_1.getTempBranches().size() + component_2.getTempBranches().size();
    in_nodes.insert(in_nodes.end(), component_1.getOutNodes().begin(), component_1.getOutNodes().end());
    in_nodes.insert(in_nodes.end(), component_2.getOutNodes().begin(), component_2.getOutNodes().end());
    temp_branches.insert(temp_branches.end(), component_1.getTempBranches().begin(), component_1.getTempBranches().end());
    temp_branches.insert(temp_branches.end(), component_2.getTempBranches().begin(), component_2.getTempBranches().end());
    _buildAncestralSelectionGraph();    
}


void AncestralSelectionGraphComponent::_coalesce(double height)
{
    Node* newNode = new Node(height);
    Branch* newBranch = new Branch; 

    shuffle(temp_branches.begin(), temp_branches.end(), random_engine);
    Branch* sampled1 = temp_branches[0];
    Branch* sampled2 = temp_branches[1];
    // cout << "coalescence: " << endl;

    sampled1->setAncestor(newNode);
    sampled2->setAncestor(newNode);
    newBranch->setDescendant(newNode);
    newNode->setCBranch(newBranch);
    newNode->setLBranch(sampled1);
    newNode->setRBranch(sampled2);    

    temp_branches.pop_front();
    temp_branches.pop_front();
    temp_branches.push_back(newBranch);

    all_nodes.push_back(newNode);
    all_branches.push_back(newBranch);
}


void AncestralSelectionGraphComponent::_branch(double height)
{
    Node* newNode = new Node(height);
    Branch* newBranch1 = new Branch;
    Branch* newBranch2 = new Branch;

    uniform_int_distribution<int> u (0, temp_branches.size()-1);
    int index = u(random_engine);
    Branch* sampled = temp_branches[index];
    // cout << "branching: " + to_string(index) << endl;
    
    sampled->setAncestor(newNode);
    newNode->setLBranch(sampled);
    newNode->setCBranch(newBranch1);
    newNode->setIBranch(newBranch2);
    newBranch1->setDescendant(newNode);
    newBranch2->setDescendant(newNode);

    temp_branches.erase(temp_branches.begin()+index);
    temp_branches.push_back(newBranch1);
    temp_branches.push_back(newBranch2);

    all_nodes.push_back(newNode);
    all_branches.push_back(newBranch1);
    all_branches.push_back(newBranch2);
}


void AncestralSelectionGraphComponent::_mutate()
{
    uniform_int_distribution<int> u (0, temp_branches.size()-1);
    int index = u(random_engine);
    // cout << "mutation: " + to_string(index) << endl;
    temp_branches[index]->addMutation();

}


string AncestralSelectionGraphComponent::_event(double height)
{
    uniform_real_distribution<double> u (0, 1);
    double k = temp_branches.size();
    double prob_coalesce = (k - 1) / (k - 1 + sigma + theta);
    double prob_branch = sigma / (k - 1 + sigma + theta);
    double prob_mutate = theta / (k - 1 + sigma + theta);
    double d = u(random_engine);
    if (d < prob_mutate)
    {
        _mutate();
        return "mutation";
    }

    else if (d >= prob_mutate && d < prob_coalesce + prob_mutate)
    {
         _coalesce(height);
        return "coalescence";
    }

    else
    {
        _branch(height);
        return "branching";
    }
}


void AncestralSelectionGraphComponent::_aftercare(double height)
{
    for(int i=0; i<temp_branches.size(); i++)
    {
        Node* newNode = new Node(height);
        Branch* newBranch = new Branch;
        Branch* branch = temp_branches[0];
        branch->setAncestor(newNode);
        newBranch->setDescendant(newNode);
        newNode->setCBranch(newBranch);
        newNode->setLBranch(branch);

        temp_branches.pop_front();
        temp_branches.push_back(newBranch);

        out_nodes.push_back(newNode);
        all_nodes.push_back(newNode);
        all_branches.push_back(newBranch);
    }
}


void AncestralSelectionGraphComponent::_buildAncestralSelectionGraph()
{
    double height = height_lower;
    int k = temp_branches.size();
    if(k * (k - 1 + theta + sigma) < 1e-300 && k * (k - 1 + theta + sigma) > -1e-300 )
    {
        for (int i=0; i<k; i++)
        {
            temp_branches[i]->addLength(height_upper - height_lower);
        }
    }
    else
    {
        if (height_upper == NONE)
        {
            while (temp_branches.size() > 1)
            {
                k = temp_branches.size();
                exponential_distribution<double> exp_distr(k * (k - 1 + theta + sigma)/((double)2 * population_size));
                double time = exp_distr(random_engine);
                for (int i=0; i<k; i++)
                {
                    temp_branches[i]->addLength(time);
                }
                height += time;
                string event = _event(height);
                // cout << event<< endl;
            }
            out_nodes.push_back(temp_branches[0]->getDescendant());
        }
        else
        {
            while (1)
            {
                k = temp_branches.size();
                exponential_distribution<double> exp_distr(k * (k - 1 + theta + sigma)/((double)2 * population_size));
                double time = exp_distr(random_engine);
                if (height + time < height_upper)
                {
                    for (int i=0; i<k; i++)
                    {
                        temp_branches[i]->addLength(time);
                    }
                    height += time;
                    string event = _event(height);
                    // cout << event<< endl;
                }
                else
                {
                    for(int i=0; i<k; i++)
                    {
                        temp_branches[i]->addLength(height_upper - height);
                    }
                    height = height_upper;
                    _aftercare(height);
                    break;
                }          

            }        
        }       
    }
}


AncestralSelectionGraphComponent::~AncestralSelectionGraphComponent() {
    // cout << "deconstruction" << endl;
    temp_branches.clear();
    in_nodes.clear();
    out_nodes.clear();

    for(vector<Node*>::iterator it=all_nodes.begin(); it!=all_nodes.end(); it++) {
        delete (*it);
        (*it) = NULL;
    }

    for(vector<Branch*>::iterator it=all_branches.begin(); it!=all_branches.end(); it++) {
        delete (*it);
        (*it) = NULL;
    }

    vector<Node*>().swap(all_nodes);
    vector<Branch*>().swap(all_branches);
}
