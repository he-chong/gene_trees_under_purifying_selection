/*********************************************************************************
*                                                                                *
*  This file is contains codes that implements the forward-in-time process of   *
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
#include "ancestral_selection_graph_composer.h"

using namespace std;


class bad_relation : public exception  
{  
public:  
    const char* what()  
    {  
        return "bad relation";  
    }   
}; 


class unknown_genotype: public exception
{
public:
    const char* what()
    {
        return "unknown genotype";
    }
};


class bad_branch_pattern : public exception  
{  
public:  
    const char* what()  
    {  
        return "bad branch pattern";  
    }   
}; 


class no_parent : public exception  
{  
public:  
    const char* what()  
    {  
        return "no parent";  
    }  
}; 


class no_child : public exception  
{  
public:  
    const char* what()  
    {  
        return "bad relation";  
    }  
}; 
 

class null_ua : public exception  
{  
public:  
    const char* what()  
    {  
        return "null UA";  
    }  
}; 


class unexpected_ua : public exception  
{  
public:  
    const char* what()  
    {  
        return "unexpected UA";  
    }  
};


class empty_component : public exception  
{  
public:  
    const char* what()  
    {  
        return "empty component";  
    }  
}; 


void AncestralSelectionGraphComposer::__assignGenotypes(Node* node)
{
    Branch * c_branch = node->getCBranch();
    Branch * i_branch = node->getIBranch();

    if(node->getGenotype() != UNKNOWN)
    {
        // cout << "UA" << endl;
        node->setRealPBranch(c_branch);
    }
    else
    {
        if (c_branch != NULL && i_branch == NULL)
        {
            if(c_branch->getAncestor()->getGenotype() == UNKNOWN) 
            {
                throw unknown_genotype();
            }
            else 
            {
                int c_genotype = c_branch->getAncestor()->getGenotype();
                node->setRealPBranch(c_branch);
                node->setGenotype(c_genotype);
            }
        }
        else if (c_branch != NULL && i_branch != NULL)
        {
            if(c_branch->getAncestor()->getGenotype() == UNKNOWN || i_branch->getAncestor()->getGenotype() == UNKNOWN)
            {
                throw unknown_genotype();
            }
            else
            {
                int c_genotype = c_branch->getAncestor()->getGenotype();
                int i_genotype = i_branch->getAncestor()->getGenotype();
                if (i_genotype == ADVANTAGEOUS)
                {
                    node->setRealPBranch(i_branch);
                    node->setGenotype(i_genotype);
                }
                else
                {
                    node->setRealPBranch(c_branch);
                    node->setGenotype(c_genotype);
                }
            }
        }
        else
        {
            throw bad_branch_pattern();
        }

        Branch* real_p_branch = node->getRealPBranch();
        Node* real_parent = real_p_branch->getAncestor();
        if(real_parent->getLBranch() == real_p_branch)
        {
            // cout << "set real l" << endl;
            real_parent->setRealLBranch(real_p_branch);
        }
        else if(real_parent->getRBranch() == real_p_branch)
        {
            // cout << "set real r" << endl;
            real_parent->setRealRBranch(real_p_branch);
        }
        else
        {
            throw bad_relation();
        }
    }

    for (int i=0; i<node->getRealPBranch()->getMutation(); i++)
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
}


void AncestralSelectionGraphComposer::__trim(Node * node)
{
    if (node->getName() == UNNAMED && node->getRealLBranch() == NULL && node->getRealRBranch() == NULL)
    {
        // cout << "trim" << endl;
        Branch* real_p_branch = node->getRealPBranch();
        Node *real_parent = node->getRealPBranch()->getAncestor();
        if (real_p_branch == real_parent->getRealLBranch())
        {
            real_parent->setRealLBranch(NULL);
        }
        else if (real_p_branch == real_parent->getRealRBranch())
        {
            real_parent->setRealRBranch(NULL);
        }
        else
        {
            throw no_parent();
        }
        __trim(real_parent);
    }
}


AncestralSelectionGraphComposer::__OutputResult* AncestralSelectionGraphComposer::__output(Node * node)
{
    // cout <<"output"<<endl;
    double height = node->getHeight();

    if (node->getName() != UNNAMED)
    {
        __OutputResult* returned = new __OutputResult(node->getName(), height);
        return returned;
    }
    else
    {
        Branch * real_l_branch = node->getRealLBranch();
        Branch * real_r_branch = node->getRealRBranch();
        if (real_l_branch != NULL && real_r_branch != NULL)
        { 
            stringstream ss;
            string last;
            Node* l_child = real_l_branch->getDescendant();
            Node* r_child = real_r_branch->getDescendant();
            __OutputResult* l_result = __output(l_child);
            __OutputResult* r_result = __output(r_child);

            ss << "(";
            ss << l_result->getLast();
            ss << ":";
            ss << l_result->getDistanceToLast() + node->getHeight() - l_child->getHeight();
            ss << ",";
            ss << r_result->getLast();
            ss << ":";
            ss << r_result->getDistanceToLast() + node->getHeight() - r_child->getHeight();
            ss << ")";
            ss >> last;
            __OutputResult* returned = new __OutputResult(last, 0);
            return returned;
        }
        if (real_l_branch != NULL && real_r_branch == NULL)
        {
            Node* l_child = real_l_branch->getDescendant();
            __OutputResult* l_result = __output(l_child);
            __OutputResult* returned = new __OutputResult(l_result->getLast(), l_result->getDistanceToLast() + node->getHeight() - l_child->getHeight());
            return returned;
        }
        else if (real_r_branch != NULL && real_l_branch == NULL)
        {
            Node* r_child = real_r_branch->getDescendant();
            __OutputResult* r_result = __output(r_child);
            __OutputResult* returned = new __OutputResult(r_result->getLast(), r_result->getDistanceToLast() + node->getHeight() - r_child->getHeight());
            return returned;
        }
        else
        {
            throw no_child();
        }
    }
}


string AncestralSelectionGraphComposer::pick_embedded_genealogy()
{
    __setComponents();
    __setUltimateAncestor();
    
    if(ultimate_ancestor == NULL)
    {
        throw null_ua();
    }

    if(components.empty())
    {
        throw empty_component();
    }

    if(ultimate_ancestor != components[0]->getOutNodes()[0])
    {
        throw unexpected_ua();
    }

    // Branch * ultimate_branch = ultimate_ancestor->getRealPBranch();
    uniform_real_distribution<double> u (0, 1);
    double x = u(composer_random_engine);
    double selection_coefficient = components[0]->getSelectionCoefficient();
    double mutation_rate = components[0]->getMutationRate();
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

    for(vector<AncestralSelectionGraphComponent*>::iterator cit=components.begin(); cit!=components.end(); cit++)
    {
        vector<Node*>& all_nodes = (*cit)->getAllNodes();
        for(vector<Node*>::reverse_iterator ait=all_nodes.rbegin(); ait!=all_nodes.rend(); ait++) 
        {
            __assignGenotypes(*ait);
        }
    }
    // cout << "assign ok" << endl;

    for(vector<AncestralSelectionGraphComponent*>::iterator cit=components.begin(); cit!=components.end(); cit++)
    {
        vector<Node*>& all_nodes = (*cit)->getAllNodes();
        for(vector<Node*>::reverse_iterator ait=all_nodes.rbegin(); ait!=all_nodes.rend(); ait++)
        {
            if((*ait)->getName() == UNNAMED && (*ait)->getRealLBranch()==NULL && (*ait)->getRealRBranch()==NULL)
            {
                __trim(*ait);
            }
        }
    }
    // cout << "trim ok" << endl;

    __OutputResult* ua_result = __output(ultimate_ancestor);
    string newick_tree = ua_result->getLast()+";";
    // cout << "output ok" << endl;

    return newick_tree;
}


AncestralSelectionGraphComposer::~AncestralSelectionGraphComposer()
{
    ultimate_ancestor = NULL;
    for(vector<AncestralSelectionGraphComponent*>::iterator it=components.begin(); it!=components.end(); it++)
    {
        (*it)=NULL;
    }
    components.clear();
    vector<AncestralSelectionGraphComponent*>().swap(components);
}