#include "InterferenceSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include <assert.h>

#include "Macros/Normalization_8TeV.h"

InterferenceSmearer::InterferenceSmearer(Normalization_8TeV* norm, float correction, float error) : norm_(norm), correction_(correction), error_(error)
{
  name_="InterferenceSmearer";
}

InterferenceSmearer::~InterferenceSmearer()
{
}

bool InterferenceSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  int genMassPoint;

  if( sample_type >= 0 ) { return true; }
  genMassPoint = std::round(norm_->GetMass(sample_type));

  if( norm_->GetProcess(sample_type) != "ggh" ) {
	  return true;
  }
  if( genMassPoint > 150 ) { genMassPoint=150; } // Warning: missing k-factor
  if( genMassPoint == 100 ) { genMassPoint=105; }  // Warning: missing k-factor
  
  assert( genMassPoint % 5 == 0 );
  
  weight = 1.-(correction_+syst_shift*error_);

  return true;
}


