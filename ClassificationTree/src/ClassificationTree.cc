#include "../interface/ClassificationTree.h"

namespace Classification{

Class ClassNode::id = 0;

unsigned int DecisionNode::makeClasses(){
	unsigned int classesMade = 0;
	if(!left){
		OUT(" Making Class to the left");
		left = new ClassNode();
		classesMade++;
	}else{
		classesMade += left->makeClasses();
	}

	if(!right){
		OUT(" Making Class to the right");
		right = new ClassNode();
		classesMade++;
	}else{
		classesMade += right->makeClasses();
	}
	return classesMade;
}


//namespace
}

