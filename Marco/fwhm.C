
//effsigma function from Chris
Double_t fwhm(TF1 * f)
{

  float y=f->GetMaximum();
  float x=f->GetMaximumX();

  float low1=100;
  float low2=120;
  float high1=120;
  float high2=140;
  float low=0;
  float high=0;

  for (int i=1; i<300; i++) {
    low=(low1+low2)/2.;
    high=(high1+high2)/2.;


    //cout<<fabs(f->Eval(low)-y/2.) <<" "<<fabs( f->Eval(high)-y/2. )<<endl;
    if(fabs(f->Eval(low)-y/2. )<0.000001 && fabs( f->Eval(high)-y/2.  )<0.000001)
      break;

    //std::cout<<"low, high, "<<low<<" "<<high<<" "<<f->Eval(low)<<" "<<f->Eval(high)<<" "<<y/2.<<endl;

    if(f->Eval(low)>y/2.) 
      low2=low;
    else
      low1=low;

    if(f->Eval(high)>y/2.) 
      high1=high;
    else
      high2=high;
  }

  return (high-low)/2.35/120*100;

}
