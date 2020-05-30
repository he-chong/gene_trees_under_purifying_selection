/*****************************************************************************
*                                                                            *
*  This file is contains codes that implements the core of the Ancestral     *
*  Selection Graph                                                           *
*                                                                            *
*  @author   何充 He Chong                                                   *
*  @email    biohe@foxmail.com                                               *
*                                                                            *
*****************************************************************************/

#include <iostream>
#include <algorithm>
#include <sstream>
#include <exception>
#include <iomanip>
#include "base_ancestral_selection_graph.h"

using namespace std;

class bad_relation : public std::exception  
{  
public:  
    const char* what()  
    {  
        return "bad relation";  
    }   
}; 

class bad_branch_pattern : public std::exception  
{  
public:  
    const char* what()  
    {  
        return "bad branch pattern";  
    }   
}; 

class no_parent : public std::exception  
{  
public:  
    const char* what()  
    {  
        return "no parent";  
    }  
}; 

class no_child : public std::exception  
{  
public:  
    const char* what()  
    {  
        return "bad relation";  
    }  
}; 

class not_bifurcation : public std::exception  
{  
public:  
    const char* what()  
    {  
        return "not bifurcation";  
    }  
}; 

/* 
 * @brief The initialization of the Ancestral Selection Graph in whatever situation
 *        needs to assign population size, mutation rate, selection_coefficient.
 *        time_upbound=None means that the backward-in-time process 
 *        stops when there exist only one lineage, i.e. the ultimate 
 *        ancestor, UA. time_upbound=t means that the backward-
 *        in-time process stops at the time t.
 *
 */
BaseAncestralSelectionGraph::BaseAncestralSelectionGraph(string name, int population_size, double mutation_rate, double selection_coefficient, double time_upbound):
    random_engine(2280714801L^chrono::high_resolution_clock::now().time_since_epoch().count())
{
    this->name = name;
    this->population_size = population_size;
    this->mutation_rate = mutation_rate;
    this->selection_coefficient = selection_coefficient;
    this->theta = 2*mutation_rate*population_size; // theta is the expected number of mutation events in a coalescent time unit. 
    double fitness_advantage = selection_coefficient < 0 ? fabs(selection_coefficient)/(1+selection_coefficient) : selection_coefficient;
    this->sigma = 2*(fitness_advantage)*population_size; // sigma is the expected number of branching events in a coalescent time unit. 
    this->time_upbound = time_upbound;
}


BaseAncestralSelectionGraph::BaseAncestralSelectionGraph():
    random_engine(2280714801L^chrono::high_resolution_clock::now().time_since_epoch().count())
{
}


void BaseAncestralSelectionGraph::_buildAncestralSelectionGraph()
{
    double total_time = 0;
    int k = temp_branches.size();
    if(k * (k - 1 + theta + sigma) < 1e-300 && k * (k - 1 + theta + sigma) > -1e-300 )
    {
        for (int i=0; i<k; i++)
        {
            temp_branches[i]->addLength(time_upbound);
        }
    }
    else
    {
        if (time_upbound == NONE)
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
                string event = _event();
                // cout << event<< endl;
            }
            Branch * ultimate_branch = temp_branches[0];
            ultimate_ancestor = ultimate_branch->getDescendant();
            ultimate_ancestor->setRealBranch(ultimate_branch);
        }
        else
        {
            while (1)
            {
                k = temp_branches.size();
                exponential_distribution<double> exp_distr(k * (k - 1 + theta + sigma)/((double)2 * population_size));
                double time = exp_distr(random_engine);
                total_time += time;
                if (total_time <= time_upbound)
                {
                    for (int i=0; i<k; i++)
                    {
                        temp_branches[i]->addLength(time);
                    }
                    string event = _event();
                    // cout << event<< endl;
                }
                else
                {
                    for (int i=0; i<k; i++)
                    {
                        temp_branches[i]->addLength(time - total_time + time_upbound);
                    }
                    break;
                }          

            }        
            ultimate_ancestor = NULL;
        }       
    }

}


string BaseAncestralSelectionGraph::pick_embedded_genealogy()
{
    Branch * ultimate_branch = ultimate_ancestor->getRealBranch();
    uniform_real_distribution<double> u (0, 1);
    double x = u(random_engine);
    // cout << selection_coefficient << endl;
    if(selection_coefficient == 0)
    {
        // cout << "neutral" << endl;
        ultimate_ancestor->setGenotype(ADVANTAGEOUS);
    }
    else
    {
    	double fitness_disadvantage = selection_coefficient > 0 ? selection_coefficient/(1+selection_coefficient) : fabs(selection_coefficient);
        // cout << "else: selection_coefficient: " << selection_coefficient << endl;
        // cout << "fitness_disadvantage: " << fitness_disadvantage << endl;
        // cout << "mutation_rate/fitness_disadvantage: " << mutation_rate/fitness_disadvantage << endl;
        if (x >= mutation_rate/fitness_disadvantage) // Sampling the genotype of the UA from the mutation-selection balance
        {
            ultimate_ancestor->setGenotype(ADVANTAGEOUS);
            // cout << "ADVANTAGEOUS" << endl;
        }
        else
        {
            ultimate_ancestor->setGenotype(DISADVANTAGEOUS);
            // cout << "DISADVANTAGEOUS" << endl;
        }
    }

    deque<Node*> queue;
    queue.push_back(ultimate_ancestor);

    while (queue.size() > 0)
    {
        Node* node = queue[0];
        queue.pop_front();
        _assignGenotype(node);
        if (node->getLBranch() != NULL)
        {
            queue.push_back(node->getLBranch()->getDescendant());
        }
        if (node->getRBranch() != NULL)
        {
            queue.push_back(node->getRBranch()->getDescendant());
        }
    }

    queue.push_back(ultimate_ancestor);
    while (queue.size() > 0)
    {
        Node* node = queue[0];
        queue.pop_front();
        _trim(node);
        if (node->getLBranch() != NULL)
        {
            queue.push_back(node->getLBranch()->getDescendant());
        }
        if (node->getRBranch() != NULL)
        {
            queue.push_back(node->getRBranch()->getDescendant());
        }
    }

    _connect(ultimate_ancestor);
    ultimate_branch->setLength(0);
    string newick_tree = _output(ultimate_branch->getDescendant())+";";
    return newick_tree;
}


void BaseAncestralSelectionGraph::_coalesce()
{
    Node* newNode = new Node;
    // cout << newNode->getGenotype() << endl;
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


void BaseAncestralSelectionGraph::_branch()
{
    Node* newNode = new Node;
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


void BaseAncestralSelectionGraph::_mutate()
{
    uniform_int_distribution<int> u (0, temp_branches.size()-1);
    int index = u(random_engine);
    // cout << "mutation: " + to_string(index) << endl;
    temp_branches[index]->addMutation();

}


string BaseAncestralSelectionGraph::_event()
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
         _coalesce();
        return "coalescence";
    }

    else
    {
        _branch();
        return "branching";
    }
}

int BaseAncestralSelectionGraph::_assignGenotype(Node * node)
{
    if (node->getGenotype() != UNKNOWN)
    {
        return node->getGenotype();
    }
    else
    {
        Branch * c_branch = node->getCBranch();
        Branch * i_branch = node->getIBranch();
        if (c_branch != NULL && i_branch == NULL)
        {
            int c_genotype = _assignGenotype(c_branch->getAncestor());
            node->setRealBranch(c_branch);
            node->setGenotype(c_genotype);
        }
        else if (c_branch != NULL && i_branch != NULL)
        {
            int c_genotype = _assignGenotype(c_branch->getAncestor());
            int i_genotype = _assignGenotype(i_branch->getAncestor());
            Branch *virtual_branch;
            if (i_genotype == ADVANTAGEOUS)
            {
                node->setRealBranch(i_branch);
                node->setGenotype(i_genotype);
                virtual_branch = c_branch;
            }
            else
            {
                node->setRealBranch(c_branch);
                node->setGenotype(c_genotype);
                virtual_branch = i_branch;
            }

            Node *virtual_parent = virtual_branch->getAncestor();
            if (virtual_parent->getLBranch() == virtual_branch)
            {
                virtual_parent->setLBranch(NULL);
            }
            else if (virtual_parent->getRBranch() == virtual_branch)
            {
                virtual_parent->setRBranch(NULL);
            }
            else
            {
                throw bad_relation();
            }
        }
        else
        {
            throw bad_branch_pattern();
        }

        for (int i=0; i<node->getRealBranch()->getMutation(); i++)
        {
            if (node->getGenotype() == ADVANTAGEOUS)
            {
                node->setGenotype(DISADVANTAGEOUS);
            }
            else
            {
                node->setGenotype(ADVANTAGEOUS);
            }
        }

        return node->getGenotype();
    }
}


void BaseAncestralSelectionGraph::_trim(Node * node)
{
    if (node->getName() == UNNAMED && node->getLBranch() == NULL && node->getRBranch() == NULL)
    {
        Node *real_parent = node->getRealBranch()->getAncestor();
        if (node->getRealBranch() == real_parent->getLBranch())
        {
            real_parent->setLBranch(NULL);
        }
        else if (node->getRealBranch() == real_parent->getRBranch())
        {
            real_parent->setRBranch(NULL);
        }
        else
        {
            throw no_parent();
        }
        _trim(real_parent);
    }
}


void BaseAncestralSelectionGraph::_connect(Node * node)
{
    Branch * l_branch = node->getLBranch();
    Branch * r_branch = node->getRBranch();
    if (l_branch != NULL && r_branch != NULL)
    { 
        Node * l_child = l_branch->getDescendant();
        Node * r_child = r_branch->getDescendant();
        _connect(l_child);
        _connect(r_child);
        
    }
    else if (l_branch != NULL && r_branch == NULL)
    {
        Node * l_child = l_branch->getDescendant();
        Branch * real_branch = node->getRealBranch();
        // Node * parent = real_branch->getAncestor();
        real_branch->setDescendant(l_child);
        real_branch->addLength(l_branch->getLength());
        l_child->setRealBranch(real_branch);
        node->setLBranch(NULL);
        l_branch->setDescendant(NULL);
        _connect(l_child);
    }
    else if (l_branch == NULL && r_branch != NULL)
    {
        Node * r_child = r_branch->getDescendant();
        Branch * real_branch = node->getRealBranch();
        // Node * parent = real_branch->getAncestor();
        real_branch->setDescendant(r_child);
        real_branch->addLength(r_branch->getLength());
        r_child->setRealBranch(real_branch);
        node->setRBranch(NULL);
        r_branch->setDescendant(NULL);
        _connect(r_child);
    }
    else
    {

    }
}


string BaseAncestralSelectionGraph::_output(Node * node)
{
    stringstream ss;
    double length = node->getRealBranch()->getLength();
    string length_str;
    ss << length;
    ss >> length_str; 

    if (node->getName() != UNNAMED)
    {
        return node->getName() + ": " + length_str;
    }
    else
    {
        Branch * l_branch = node->getLBranch();
        Branch * r_branch = node->getRBranch();
        if (l_branch != NULL && r_branch != NULL)
        { 
            if(length == 0)
            {
                return("(" + _output(l_branch->getDescendant()) + "," + _output(r_branch->getDescendant()) +")");
            }
            else
            {
                return("(" + _output(l_branch->getDescendant()) + "," + _output(r_branch->getDescendant()) +"): " + length_str);
            }

        }
        else if (l_branch != NULL && r_branch == NULL)
        {
            throw not_bifurcation();
        }
        else if (l_branch == NULL && r_branch != NULL)
        {
            throw not_bifurcation();
        }
        else
        {
            throw no_child();
        }
    }
}


BaseAncestralSelectionGraph::~BaseAncestralSelectionGraph() {
    // cout << "deconstruction" << endl;
    // cout << name <<": selection_coefficient: " << selection_coefficient << endl;
    temp_branches.clear();
    for(Node* node : all_nodes) {
        delete node;
        node = NULL;
    }
    vector<Node*>().swap(all_nodes);
    for(Branch* branch : all_branches) {
        delete branch;
        branch = NULL;
    }
    vector<Branch*>().swap(all_branches);
}
