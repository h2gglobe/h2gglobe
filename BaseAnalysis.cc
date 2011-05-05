#include "BaseAnalysis.h"

BaseAnalysis::BaseAnalysis()
{}

BaseAnalysis::~BaseAnalysis() 
{}

bool operator == (BaseAnalysis * a, const std::string & b) { return a->name() == b; };
