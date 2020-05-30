#include <iostream>
#include "for_ancestral_selection_graph.h"
#define UNKNOWN -1
#define UNNAMED ""

Node::Node()
{
    name = UNNAMED;
    genotype = UNKNOWN;
    l_branch = NULL;
    r_branch = NULL;
    c_branch = NULL;
    i_branch = NULL;
    real_branch = NULL;
}

Node::Node(string name)
{
    this->name = name;
    genotype = UNKNOWN;
    l_branch = NULL;
    r_branch = NULL;
    c_branch = NULL;
    i_branch = NULL;
    real_branch = NULL;
}

Node::~Node()
{
    l_branch = NULL;
    r_branch = NULL;
    c_branch = NULL;
    i_branch = NULL;
    real_branch = NULL;
}

string Node::getName()
{
    return name;
}

void Node::setName(string name)
{
    this->name = name;
}

int Node::getGenotype()
{
    return genotype;
}

void Node::setGenotype(int genotype)
{
    this->genotype = genotype;
}

Branch * Node::getLBranch()
{
    return l_branch;
}

Branch * Node::getRBranch()
{
    return r_branch;
}

Branch * Node::getCBranch()
{
    return c_branch;
}

Branch * Node::getIBranch()
{
    return i_branch;
}

Branch * Node::getRealBranch()
{
    return real_branch;
}

void Node::setLBranch(Branch * aBranch)
{
    l_branch = aBranch;
}

void Node::setRBranch(Branch * aBranch)
{
    r_branch = aBranch;
}

void Node::setCBranch(Branch * aBranch)
{
    c_branch = aBranch;
}

void Node::setIBranch(Branch * aBranch)
{
    i_branch = aBranch;
}

void Node::setRealBranch(Branch * aBranch)
{
    real_branch = aBranch;
}


Branch::Branch()
{
    mutation = 0;
    length = 0;
    ancestor = NULL;
    descendant = NULL;
}

Branch::~Branch()
{
    ancestor = NULL;
    descendant = NULL;
}

Node * Branch::getAncestor()
{
    return ancestor;
}

Node * Branch::getDescendant()
{
    return descendant;
}

void Branch::setAncestor(Node * ancestor)
{
    this->ancestor = ancestor;
}

void Branch::setDescendant(Node * descendant)
{
    this->descendant = descendant;
}

int Branch::getMutation()
{
    return mutation;
}

double Branch::getLength()
{
    return length;
}

void Branch::addMutation()
{
    mutation ++;
}

void Branch::addLength(double length)
{
    this->length += length;
}

void Branch::setLength(double length)
{
    this->length = length;
}