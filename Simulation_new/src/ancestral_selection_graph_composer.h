#ifndef ANCESTRAL_SELECTION_GRAPH_COMPOSER_H
#define ANCESTRAL_SELECTION_GRAPH_COMPOSER_H
#define UNKNOWN -1
#define NONE -1
#define UNNAMED ""
#define ADVANTAGEOUS 1
#define DISADVANTAGEOUS 2
#include <vector>
#include <stack>
#include <string>
#include <random>
#include <chrono>
#include "ancestral_selection_graph_component.h"

using namespace std;


class AncestralSelectionGraphComposer
{
    protected:
        Node* ultimate_ancestor;
        vector<AncestralSelectionGraphComponent*> components;
        mt19937 composer_random_engine;

        class __OutputResult
        {
            private:
                double distance_to_last;
                string last;
            public:
                __OutputResult(string l, double dtl)
                {
                    last = l;
                    distance_to_last = dtl;
                }
                double getDistanceToLast() {return distance_to_last;}
                void setDistanceToLast(double dtl) {distance_to_last = dtl;}
                string getLast() {return last;}
                void setLast(string l) {last = l;}
        };
        
        void __assignGenotypes(Node* node);
        void __trim(Node* node);
        __OutputResult* __output(Node* node); 

        virtual void __setUltimateAncestor()=0;
        virtual void __setComponents()=0;

    public:
        AncestralSelectionGraphComposer():ultimate_ancestor(NULL), composer_random_engine(2280714801L^chrono::high_resolution_clock::now().time_since_epoch().count())
        {
        }
        virtual ~AncestralSelectionGraphComposer();

        string pick_embedded_genealogy();

        Node* getUltimateAncestor() {
            return ultimate_ancestor;
        }

        virtual void getExtantLineages(vector<Node*>& extantLineages)=0;

};


#endif