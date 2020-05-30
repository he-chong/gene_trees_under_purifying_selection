#include <iostream>
#include "for_ancestral_selection_graph.h"


Node::Node(string name, double height)
{
    this->name = name;
    this->height = height;
    genotype = UNKNOWN;
    l_branch = NULL;
    r_branch = NULL;
    c_branch = NULL;
    i_branch = NULL;
    real_p_branch = NULL;
    real_l_branch = NULL;
    real_r_branch = NULL;
}

Node::Node(double height)
{
    this->name = UNNAMED;
    this->height = height;
    genotype = UNKNOWN;
    l_branch = NULL;
    r_branch = NULL;
    c_branch = NULL;
    i_branch = NULL;
    real_p_branch = NULL;
    real_l_branch = NULL;
    real_r_branch = NULL;
}

Node::~Node()
{
    l_branch = NULL;
    r_branch = NULL;
    c_branch = NULL;
    i_branch = NULL;
    real_p_branch = NULL;
    real_l_branch = NULL;
    real_r_branch = NULL;
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

double Node::getHeight()
{
    return height;
}

void Node::setHeight(double height)
{
    this->height = height;
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

Branch * Node::getRealPBranch()
{
    return real_p_branch;
}

Branch * Node::getRealLBranch()
{
    return real_l_branch;
}

Branch * Node::getRealRBranch()
{
    return real_r_branch;
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

void Node::setRealPBranch(Branch * aBranch)
{
    real_p_branch = aBranch;
}

void Node::setRealLBranch(Branch * aBranch)
{
    real_l_branch = aBranch;
}

void Node::setRealRBranch(Branch * aBranch)
{
    real_r_branch = aBranch;
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