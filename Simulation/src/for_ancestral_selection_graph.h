#ifndef FOR_ANCESTRAL_SELECTION_GRAPH_H
#define FOR_ANCESTRAL_SELECTION_GRAPH_H
#include <vector>
#include <string>
using namespace std;

class Branch;

class Node
{
    private:
        string name;
        int genotype;
        double height;
        Branch *c_branch, *i_branch, *real_branch;
        Branch *l_branch, *r_branch;
        
    public:
        Node();
        ~Node();
        Node(string name);

        string getName();
        void setName(string name);
        string getNodeName();
        void setNodeName(string nodeName);        

        int getGenotype();
        void setGenotype(int genotype);
        
        Branch * getLBranch();
        Branch * getRBranch();
        Branch * getCBranch();
        Branch * getIBranch();
        Branch * getRealBranch();

        void setLBranch(Branch * aBranch);
        void setRBranch(Branch * aBranch);
        void setCBranch(Branch * aBranch);
        void setIBranch(Branch * aBranch);
        void setRealBranch(Branch * aBranch);

};

class Branch
{
    private:
        int mutation;
        double length;
        Node *ancestor, *descendant;

    public:
        Branch();
        ~Branch();
        Node * getAncestor();
        Node * getDescendant();
        void setAncestor(Node * ancestor);
        void setDescendant(Node * descendant);

        int getMutation();
        double getLength();
        void addMutation();
        void addLength(double length);
        void setLength(double length);
};

#endif