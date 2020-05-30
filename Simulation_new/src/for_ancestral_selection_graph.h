#ifndef FOR_ANCESTRAL_SELECTION_GRAPH_H
#define FOR_ANCESTRAL_SELECTION_GRAPH_H
#define UNKNOWN -1
#define UNNAMED ""
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
        Branch *c_branch, *i_branch, *real_p_branch;
        Branch *l_branch, *r_branch, *real_l_branch, *real_r_branch;
        
    public:
        Node(string name, double height);
        Node(double height);
        ~Node();

        string getName();
        void setName(string name);   

        int getGenotype();
        void setGenotype(int genotype);
        
        double getHeight();
        void setHeight(double);
        
        Branch * getLBranch();
        Branch * getRBranch();
        Branch * getCBranch();
        Branch * getIBranch();
        Branch * getRealPBranch();
        Branch * getRealLBranch();
        Branch * getRealRBranch();

        void setLBranch(Branch * aBranch);
        void setRBranch(Branch * aBranch);
        void setCBranch(Branch * aBranch);
        void setIBranch(Branch * aBranch);
        void setRealPBranch(Branch * aBranch);
        void setRealLBranch(Branch * aBranch);
        void setRealRBranch(Branch * aBranch);

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