#ifndef __BASESMEARER__
#define __BASESMEARER__

#include "LoopAll.h"

// ------------------------------------------------------------------------------------
/**
 * \class BaseSmearer
 *
 * Defines the interface to be implemented by any smearing function on photon informations
 *
 *
 */

class PhotonReducedInfo;

class BaseSmearer 
{
public:
	// ! C-TOR
	BaseSmearer();
	
	// ! D-TOR
	virtual ~BaseSmearer();
	
	// ! the class name
	virtual const std::string & name() const = 0; 
	
	// !  
	operator const std::string & () const { return this->name(); };
	
	// ! return smeared photon informations
	virtual bool smearPhoton( PhotonReducedInfo & pho, float & weight, int run, float syst_shift=0. ) const = 0;

	/// virtual bool smearDiPhoton( TLorentzVector & 4p, const TVector3 & selVtx, const TVector3 & trueVtx);
};

// ! Used to search analyzers by name 
bool operator == (BaseSmearer * a, const std::string & b);


class BaseDiPhotonSmearer 
{
public:
	// ! C-TOR
	BaseDiPhotonSmearer();
	
	// ! D-TOR
	virtual ~BaseDiPhotonSmearer();
	
	// ! the class name
	virtual const std::string & name() const = 0; 
	
	// !  
	operator const std::string & () const { return this->name(); };
	
	virtual bool smearDiPhoton( TLorentzVector & p4, TVector3 & selVtx, float & weight, const int & category, 
				    const int & genMassPoint, const TVector3 & trueVtx, float & idMVA1, float & idMVA2, float syst_shift=0.) const = 0 ;
};

// ! Used to search analyzers by name 
bool operator == (BaseDiPhotonSmearer * a, const std::string & b);


class BaseGenLevelSmearer 
{
public:
	// ! C-TOR
	BaseGenLevelSmearer();
	
	// ! D-TOR
	virtual ~BaseGenLevelSmearer();
	
	// ! the class name
	virtual const std::string & name() const = 0; 
	
	// !  
	operator const std::string & () const { return this->name(); };
	
	virtual bool smearEvent(  float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift=0.) const = 0 ;
};

// ! Used to search analyzers by name 
bool operator == (BaseGenLevelSmearer * a, const std::string & b);

#endif
