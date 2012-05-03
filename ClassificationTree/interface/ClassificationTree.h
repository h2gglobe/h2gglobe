#ifndef CLASSTREE_H
#define CLASSTREE_H

#define CLSDEBUGOUT
#define CLSASSERT

#ifdef CLSDEBUG
#include <iostream>
#define OUT(str) do { std::cout << __FUNCTION__ << ": " <<  str << std::endl; } while( false )
#else
#define OUT(str) do { } while ( false )
#endif

#ifdef CLSASSERT
#include <iostream>
#include <cstdlib>
#   define ASSERT(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            std::exit(-1);						\
        } \
    } while (false)
#else
#   define ASSERT(condition, message) do { } while (false)
#endif

#include <vector>

namespace Classification{

using namespace std;
typedef unsigned int Dimension;
typedef double Value;
typedef unsigned int Class;
typedef vector<Value> Point;

///////////////////////////////////////////
class DecisionNode{
public:
	DecisionNode(Dimension dimension_, Value value_)
	: dimension(dimension_), value(value_), left(0), right(0), parent(0){
		OUT("(" << dimension << "," << value << ")");
	};

	DecisionNode()
	: dimension(0), value(0), left(0), right(0), parent(0){
	};

	Dimension dimension;
	Value value;
	DecisionNode *left;   // Pointer to the left subtree.
	DecisionNode *right;  // Pointer to the right subtree.
	DecisionNode *parent;

	virtual unsigned int makeClasses();

	DecisionNode *addNodeLeft(DecisionNode *node){
		return addNode(node, true);
	}
	DecisionNode *addNodeRight(DecisionNode *node){
		return addNode(node, false);
	}

	DecisionNode *addNode(DecisionNode *nodeToAdd, bool addToLeft){
		// (Temporarily) assign
		// Makes the checking simpler as the parent will know
		// if the new node is to be added left or right.
		DecisionNode *&targetNode = addToLeft ? left : right;
		ASSERT( targetNode==0, "Cannot re-assign node.");
		targetNode = nodeToAdd;

		OUT("Check if can add node");
		// Check if can be added
		if( !canAddNode(nodeToAdd, nodeToAdd) ){
			OUT("Could not add a node.");
			targetNode = 0;
			return 0;
		}

		nodeToAdd->parent = this;
		return nodeToAdd;
	}

	bool canAddNode(DecisionNode *nodeToAdd, DecisionNode *daughter){
		OUT(" this: " <<"(" << dimension << "," << value << ")" << " "
				<< " toAdd: " "(" << nodeToAdd->dimension << "," << nodeToAdd->value << ")" << " "
				<< " daughter: " "(" << daughter->dimension << "," << daughter->value << ")"
				);
		// Condition in which we cannot add:
		// same dimension and smaller (larger) value on the right (left)
		if(dimension == nodeToAdd->dimension
				&& (
						(daughter == right && nodeToAdd->value < value)
						||
						(daughter == left  && nodeToAdd->value >= value)
				)
		){
			ASSERT( false, "Could not add a node."  );
			return false;
		}

		if(parent==0){
			OUT(" Reached the root." << endl );
			// if the root node is ok, then trickle back that it is ok.
			return true;
		}else{
			OUT(" Calling parent.");
			// otherwise pass it to the parent telling who we are
			return parent->canAddNode(nodeToAdd, this);
		}
	}

	virtual Class Classify(Point& point){
		OUT("(" << dimension << "," << value << ")");
		if (point[dimension] < value){
			OUT(" go left");
			return left->Classify(point);
		}else{
			OUT(" go right");
			return right->Classify(point);
		}
	}


};

///////////////////////////////////////////
class ClassNode: public DecisionNode{
public:
	ClassNode(){
		myid = id++;
		OUT(myid);
	};

	static Class id;
	Class myid;

	Class Classify(Point& point){
		OUT(myid);
		return myid;
	}
};

///////////////////////////////////////////
class Tree: public DecisionNode{
public:
	Tree()
	: nodes(0){
	};

	vector<DecisionNode *> nodes;

	int addNode(Dimension dimension, Value value, int parent, bool addToLeft=true){
		if(parent==-1){
			ASSERT(nodes.size() == 0, "Cannot redo the root.");
		}else{
			ASSERT(nodes.size() > parent, "Trying to add but parent not yet there");
		}

		nodes.push_back( new DecisionNode(dimension, value) );
		if(parent!=-1){
			nodes[parent]->addNode(nodes.back(), addToLeft);
		}
		return nodes.size();
	}

	unsigned int makeClasses(){
		ClassNode::id=0;
		return nodes[0]->makeClasses();
	};

	Class Classify(Point& point){
		return nodes[0]->Classify(point);
	}
};


//namespace
};
#endif
