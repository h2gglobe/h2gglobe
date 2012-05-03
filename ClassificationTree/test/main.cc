#include <iostream>
#include "../interface/ClassificationTree.h"

using namespace Classification;
using namespace std;

void testAPoint(Tree &tree, double x, double y){
	Point test_point(2);

	test_point[0]=x;
	test_point[1]=y;
	cout << "Classifying (" << test_point[0] << "," << test_point[1] << ")" <<endl;
	cout << " belongs to class " << tree.Classify(test_point) << endl;
}

int main(void){

	cout << "Make tree" << endl;
	Tree tree;
	tree.addNode(0,1.4  ,-1);
	tree.addNode(1,0.08 ,0,true);
	tree.addNode(0,1.0  ,1,false);
	tree.addNode(1,0.04 ,0,false);

	cout << "Counting Classes" << endl;
	int leaves = tree.makeClasses();
	cout << " found " << leaves << " classes" << endl;
	cout<<endl;

	testAPoint(tree, 0.9, 0.02);
	testAPoint(tree, 1.3, 0.20);
	testAPoint(tree, 1.7, 0.20);

	return 0;
}
