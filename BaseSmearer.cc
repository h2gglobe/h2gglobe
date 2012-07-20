#include "BaseSmearer.h"

int BaseSmearer::nRegisteredSmerers_ = 0;

BaseSmearer::BaseSmearer():
	smearerId_(-1)
{}

BaseSmearer::~BaseSmearer() 
{}

bool operator == (BaseSmearer * a, const std::string & b) { return a->name() == b; };

void BaseSmearer::registerMe()
{
	smearerId_ = nRegisteredSmerers_;
	++nRegisteredSmerers_;
}


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

