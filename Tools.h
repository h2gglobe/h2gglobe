#ifndef Tools_h
#define Tools_h

#include "TRandom.h"

#define PI 3.141592654
#define TWOPI 6.283185308

//ADDED MJ
Float_t pho_Et[MAX_PHOTONS];

struct PhotonCandidate{
    TLorentzVector *p4;
    TVector3 *calopos;
    bool pixSeed;
    double trkIso;
    double ecalIso;
    double hcalIso;
    double sieie;
    double hoe;
    double r9;
};

//Added NCKW
bool PhoP4greater(PhotonCandidate p1, PhotonCandidate p2){
	
     return p1.p4->Pt() > p2.p4->Pt();

}

float DeltaPhi(float phi1, float phi2) {
  float deltaphi;
  if(phi1<0) phi1+=TWOPI;
  if(phi2<0) phi2+=TWOPI;
  deltaphi=fabs(phi1-phi2);
  if(deltaphi>TWOPI) deltaphi-=TWOPI;
  if(deltaphi>PI) deltaphi=TWOPI-deltaphi;
  return deltaphi;
}

float DeltaPhiWithSign(float phi1, float phi2) {
  float deltaphi;
  if(phi1<0) phi1+=TWOPI;
  if(phi2<0) phi2+=TWOPI;
  deltaphi=(phi1-phi2);
  if(deltaphi>TWOPI) deltaphi-=TWOPI;
  if(deltaphi<-TWOPI) deltaphi+=TWOPI;
  if(deltaphi>PI) deltaphi=-(TWOPI-deltaphi);
  if(deltaphi<-PI) deltaphi=(TWOPI+deltaphi);
  return deltaphi;
}

float DeltaEtaWithSign(float eta1, float eta2) {
  return eta1-eta2;
}

float DeltaTheta(float the1, float phi1, float the2, float phi2) {
  float a;
  float x[2][3];
  x[0][0]=sin(the1)*cos(phi1);
  x[0][1]=sin(the1)*sin(phi1);
  x[0][2]=cos(the1);
  x[1][0]=sin(the2)*cos(phi2);
  x[1][1]=sin(the2)*sin(phi2);
  x[1][2]=cos(the2);
  a=acos(x[0][0]*x[1][0]+x[0][1]*x[1][1]+x[0][2]*x[1][2]);
  return a;
}

float DeltaR(float phi1, float phi2, float eta1, float eta2) {
  float deltaphi=DeltaPhi(phi1,phi2);
  float deltaeta=fabs(eta1-eta2);
  float delta = sqrt(deltaphi*deltaphi+ deltaeta*deltaeta);
  //  cout<<delta<<" "<<deltaeta<<" "<<deltaphi<<endl;
  return delta;
}

float DeltaR(float phi1, float phi2, float eta1, float eta2, float r1, float z1orig, float z1new) {

  float x=0.;
  float y=0.;
  float z0=0.;
  float the1=2*atan(exp(-eta1));

  x = r1*sin(the1)*cos(phi1);
  y = r1*sin(the1)*sin(phi1);
  z0 = r1*cos(the1);		     
  z0=z0-z1new;
  the1=acos(z0/sqrt(x*x+y*y+z0*z0));
  float eta1new=-log(tan(the1/2.));

  float deltaphi=DeltaPhi(phi1,phi2);
	       //float deltaeta=fabs(eta1-eta2);
  float deltaeta=fabs(eta1new-eta2);
  float delta = sqrt(deltaphi*deltaphi+ deltaeta*deltaeta);
  //  cout<<delta<<" "<<deltaeta<<" "<<deltaphi<<endl;
  return delta;
}

float DeltaPhiSign(float phi1, float phi2) {
  float deltaphi;
  if(phi1<0) phi1+=TWOPI;
  if(phi2<0) phi2+=TWOPI;
  deltaphi=phi1-phi2;
  if(deltaphi>PI) deltaphi=deltaphi-TWOPI;
  if(deltaphi<-PI) deltaphi=deltaphi+TWOPI;
  return deltaphi;
}


float CorrEtaVertex(float eta, float r, float zvtx) {
  float mythe=0.;
  float x=0.;
  float z0=0.;
  mythe=2*atan(exp(-eta));
  x = r*sin(mythe);
  z0 = r*cos(mythe);		     
  z0=z0-zvtx;
  mythe=acos(z0/sqrt(x*x+z0*z0));
  return -log(tan(mythe/2.));
}

float InvMass(float p1, float the1, float phi1, float p2, float the2, float phi2) {

  float px1 = 0.;
  float py1 = 0.;
  float pz1 = 0.;
  float px2 = 0.;
  float py2 = 0.;
  float pz2 = 0.;


  px1=p1*sin(the1)*cos(phi1);
  py1=p1*sin(the1)*sin(phi1);
  pz1=p1*cos(the1);
    
  px2=p2*sin(the2)*cos(phi2);
  py2=p2*sin(the2)*sin(phi2);
  pz2=p2*cos(the2);

  float invMass=sqrt( (p1+p2)*(p1+p2) - (px1+px2)*(px1+px2) - (py1+py2)*(py1+py2) - (pz1+pz2)*(pz1+pz2));
  return invMass;
}

float PtPair(float et1, float eta1, float phi1, float et2, float eta2, float phi2) {

  //float p1=0.;
  //float p2=0.;
  //float the1=2*atan(exp(-eta1));
  //float the2=2*atan(exp(-eta2));
  float px1 = 0.;
  float py1 = 0.;
  //float pz1 = 0.;
  float px2 = 0.;
  float py2 = 0.;
  //float pz2 = 0.;

  //p1=et1/sin(the1);
  //p2=et2/sin(the2);

  //px1=p1*sin(the1)*cos(phi1);
  //py1=p1*sin(the1)*sin(phi1);
  //pz1=p1*cos(the1);
    
  //px2=p2*sin(the2)*cos(phi2);
  //py2=p2*sin(the2)*sin(phi2);
  //pz2=p2*cos(the2);

  px1=et1*cos(phi1);
  py1=et1*sin(phi1);

  px2=et2*cos(phi2);
  py2=et2*sin(phi2);

  float ptPair = sqrt((px1+px2)*(px1+px2) + (py1+py2)*(py1+py2));
  return ptPair;
}

float InvMass(float et1, float eta1, float phi1, float r1, float et2, float eta2, float phi2, float r2, float zVtxFound)
{

  float px1 = 0.;
  float py1 = 0.;
  float pz1 = 0.;
  float px2 = 0.;
  float py2 = 0.;
  float pz2 = 0.;
  float p1=0.;
  float p2=0.;

  //cout<<"InvMass eta1="<<eta1<<" r1 "<<r1<<endl;
  //cout<<"InvMass eta2="<<eta2<<" r2 "<<r2<<endl;

  float the1=2*atan(exp(-eta1));
  float the2=2*atan(exp(-eta2));

  p1=et1/sin(the1);
  p2=et2/sin(the2);

  float x=0.;
  float y=0.;
  float z0=0.;

  x = r1*sin(the1)*cos(phi1);
  y = r1*sin(the1)*sin(phi1);
  z0 = r1*cos(the1);		     
  z0=z0-zVtxFound;
  the1=acos(z0/sqrt(x*x+y*y+z0*z0));
  //float eta1n=-log(tan(the1/2.));
  //if(fabs(eta1)>2.4&&fabs(zVtxFound)>12.) cout<<"*** eta1 "<<eta1<<" "<<eta1n<<" "<<zVtxFound<<" "<<fabs(eta1-eta1n)<<" the1 "<<the1<<"x,y,z0 "<<x<<" "<<y<<" "<<z0<<endl;

  x = r2*sin(the2)*cos(phi2);
  y = r2*sin(the2)*sin(phi2);
  z0 = r2*cos(the2);		     
  z0=z0-zVtxFound;
  the2=acos(z0/sqrt(x*x+y*y+z0*z0));

  //float eta2n=-log(tan(the2/2.));
  //if(fabs(eta2)>2.4&&fabs(zVtxFound)>12.) cout<<"*** eta2 "<<eta2<<" "<<eta2n<<" "<<zVtxFound<<" "<<fabs(eta2-eta2n)<<" the2 "<<the2<<"x,y,z0 "<<x<<" "<<y<<" "<<z0<<endl;

  px1=p1*sin(the1)*cos(phi1);
  py1=p1*sin(the1)*sin(phi1);
  pz1=p1*cos(the1);
    
  px2=p2*sin(the2)*cos(phi2);
  py2=p2*sin(the2)*sin(phi2);
  pz2=p2*cos(the2);

  float invMass2=( (p1+p2)*(p1+p2) - (px1+px2)*(px1+px2) - (py1+py2)*(py1+py2) - (pz1+pz2)*(pz1+pz2));
  //float invMass=sqrt( (p1+p2)*(p1+p2) - (px1+px2)*(px1+px2) - (py1+py2)*(py1+py2) - (pz1+pz2)*(pz1+pz2));
  if(invMass2<=0.) return 0.;
  float invMass=sqrt(invMass2);
  //cout<<invMass<<" "<<invMass2<<" "<<p1<<" "<<p2<<" "<<px1<<" "<<px2<<" "<<py1<<" "<<py2<<" "<<pz1<<" "<<pz2<<endl;
  return invMass;
}

int EventPid(int pid1, int pid2) {

  if(pid1==22) {
    if(pid2==22) return 0;
    if(pid2==111)return 1;
    if(pid2==221)return 2;
    if(pid2==331)return 3;
    if(pid2==113)return 4;
    if(pid2==223)return 5;
    if((int) fabs( (float) pid2)==11) return 8;
    if(pid2==130)return 9;
  }
  else if(pid1==111) {
    if(pid2==22) return 1;
    if(pid2==111)return 11;
    if(pid2==221)return 12;
    if(pid2==331)return 13;
    if(pid2==113)return 14;
    if(pid2==223)return 15;
    if((int) fabs( (float) pid2)==11) return 18;
    if(pid2==130)return 19;
  }
  else if(pid1==221) {
    if(pid2==22) return 2;
    if(pid2==111)return 12;
    if(pid2==221)return 22;
    if(pid2==331)return 23;
    if(pid2==113)return 24;
    if(pid2==223)return 25;
    if((int) fabs( (float) pid2)==11) return 28;
    if(pid2==130)return 29;
  }
  else if(pid1==331) {
    if(pid2==22) return 3;
    if(pid2==111)return 13;
    if(pid2==221)return 23;
    if(pid2==331)return 33;
    if(pid2==113)return 34;
    if(pid2==223)return 35;
    if((int) fabs( (float) pid2)==11) return 38;
    if(pid2==130)return 39;
  }
  else if(pid1==113) {
    if(pid2==22) return 4;
    if(pid2==111)return 14;
    if(pid2==221)return 24;
    if(pid2==331)return 34;
    if(pid2==113)return 44;
    if(pid2==223)return 45;
    if((int) fabs( (float) pid2)==11) return 48;
    if(pid2==130)return 49;
  }
  else if(pid1==223) {
    if(pid2==22) return 5;
    if(pid2==111)return 15;
    if(pid2==221)return 25;
    if(pid2==331)return 35;
    if(pid2==113)return 45;
    if(pid2==223)return 55;
    if((int) fabs( (float) pid2)==11) return 58;
    if(pid2==130)return 59;
  }
  else if((int) fabs((float) pid1)==11) {
    if(pid2==22) return 8;
    if(pid2==111)return 18;
    if(pid2==221)return 28;
    if(pid2==331)return 38;
    if(pid2==113)return 48;
    if(pid2==223)return 58;
    if((int) fabs( (float) pid2)==11) return 88;
    if(pid2==130)return 89;
  }
  else if(pid1==130) {
    if(pid2==22) return 9;
    if(pid2==111)return 19;
    if(pid2==221)return 29;
    if(pid2==331)return 39;
    if(pid2==113)return 49;
    if(pid2==223)return 59;
    if((int) fabs( (float) pid2)==11) return 89;
    if(pid2==130)return 99;
  }

  //cout<<"PID1 PID2 ************* "<<pid1<<" "<<pid2<<endl;
  return -1;
}

float CorrectSciso(float scisonew, int iflag)
{
  // iflag 0 for endcaps and 1 for barrel (hybrid and island) 
  // barrel eta from 1.2 to 1.5
  Double_t Fraction_endcaps=0.;
  Double_t Low=0.;
  Double_t High=0.;
  Double_t par[3];
  if (iflag==0) {
    Fraction_endcaps=0.42;
    Low=0.15;
    High=1.;
    par[0]=  8.40449e-02;
    par[1]=  3.85324e-01;                          
    par[2]=  2.61439e-01;                         
  }
  else if (iflag==1) {
    Fraction_endcaps=0.195;
    Low=0.15;
    High=1.;
    par[0]=  2.38174e-01;
    par[1]=  -3.41506e+00;
  }
  //      TRandom r;  frac=r.Rndm();
  Double_t frac=gRandom->Rndm();
  //cout<<" frac " <<frac<<endl;
  if (frac<Fraction_endcaps) {
    bool done=false;
    //int iii=1;
    while(!done) {
      Double_t x1,x2;
      x1=gRandom->Rndm();
      x2=gRandom->Rndm();
      x1=Low+x1*(High-Low);
      Double_t arg1=0.;
      Double_t f=0.;
      if(iflag==0) {
	arg1 = (x1 - par[1])/(par[2]);
	f = par[0]*exp(-0.5*pow(arg1,2));
      }
      else if(iflag==1) {
	f=par[0]*exp(par[1]*x1);
      }	  
      //cout<< x1<< " " <<x1<< " " <<x2<< " " <<f<<endl;
      if (x2<f) {
	scisonew=x1;
	done=true;
      }
    }
  }
  return scisonew;
}
	
#endif
