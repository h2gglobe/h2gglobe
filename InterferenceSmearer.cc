#include "InterferenceSmearer.h"
#include "PhotonReducedInfo.h"
#include "TRandom3.h"
#include <assert.h>

InterferenceSmearer::InterferenceSmearer(float correction, float error) : correction_(correction), error_(error)
{
  name_="InterferenceSmearer";
}

InterferenceSmearer::~InterferenceSmearer()
{
}

bool InterferenceSmearer::smearEvent( float & weight, const TLorentzVector & p4, const int nPu, const int sample_type, float syst_shift ) const 
{
  int genMassPoint;

  if      (sample_type == -1 ) genMassPoint=90;
  else if (sample_type == -5 ) genMassPoint=95;
  else if (sample_type == -9 ) genMassPoint=100;
  else if (sample_type == -13) genMassPoint=105;
  else if (sample_type == -17) genMassPoint=110;
  else if (sample_type == -21) genMassPoint=115;
  else if (sample_type == -25) genMassPoint=120;
  else if (sample_type == -29) genMassPoint=130;
  else if (sample_type == -33) genMassPoint=140;
  else if (sample_type == -37) genMassPoint=125;
  else if (sample_type == -41) genMassPoint=135;
  else if (sample_type == -45) genMassPoint=145;
  else if (sample_type == -49) genMassPoint=150;
  else if (sample_type == -53) genMassPoint=121;  
  else if (sample_type == -57) genMassPoint=123;  
  else if (sample_type == -61) genMassPoint=150;  // FIXME  this is mass=155; remapped into 150: no kfactors currently available for 155;
  else if (sample_type == -65) genMassPoint=150;  // FIXME  this is mass=160; remapped into 150: no kfactors currently available for 160;
  else if (sample_type == -69) genMassPoint=105;  // FIXME  this is mass=100; remapped into 105: no kfactors currently available for 100;

  else if (sample_type <=-62) assert(0);   // this is the case of non-existing sample
  else    return true;                     // this is the case of backgrounds

  weight = 1.-(correction_+syst_shift*error_);

  return true;
}


