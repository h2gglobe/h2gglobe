#include "BaseSmearer.h"

BaseSmearer::BaseSmearer()
{}

BaseSmearer::~BaseSmearer() 
{}

bool operator == (BaseSmearer * a, const std::string & b) { return a->name() == b; };
