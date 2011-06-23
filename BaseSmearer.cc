#include "BaseSmearer.h"

BaseSmearer::BaseSmearer()
{}

BaseSmearer::~BaseSmearer() 
{}

bool operator == (BaseSmearer * a, const std::string & b) { return a->name() == b; };

BaseDiPhotonSmearer::BaseDiPhotonSmearer()
{}

BaseDiPhotonSmearer::~BaseDiPhotonSmearer() 
{}

bool operator == (BaseDiPhotonSmearer * a, const std::string & b) { return a->name() == b; };

BaseGenLevelSmearer::BaseGenLevelSmearer()
{}

BaseGenLevelSmearer::~BaseGenLevelSmearer() 
{}

bool operator == (BaseGenLevelSmearer * a, const std::string & b) { return a->name() == b; };

