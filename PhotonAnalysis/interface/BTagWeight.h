#ifndef __BTAGWEIGHT_H__
#define __BTAGWEIGHT_H__

#include <math.h>
#include <iostream>
 
class BTagWeight 
{
 public:
  BTagWeight(int jmin, int jmax,float bEff, float cEff, float  udsgEff) : 
    maxTags(jmax), minTags(jmin),
    meanBeff(bEff), meanCeff(cEff), 
    meanUDSGeff(udsgEff) {}

    bool filter(int t);
    float weight(int b,int c,int l,float sfb,float sfc,float sfl, int tags);
    int comb(int n, int k);
    int fact(int n);

 private:
    int maxTags;
    int minTags;
    float meanBeff;
    float meanCeff;
    float meanUDSGeff;
};

float BTagWeight::weight(int b,int c,int l,float sfb,float sfc,float sfl, int tags)
{
  if(!filter(tags))
    {
      //   std::cout << "This event should not pass the selection, what is it doing here?" << std::endl;
      return 0;
    }
  int njets=b+c+l;

  /* int comb= 1 << njets;
 for(int i=0;i < comb; i++)
  {
    

  }
  */
  float pMC=0;
  float pData=0;
  for(int ib=0;ib<=b;ib++)
    for(int ic=0;ic<=c;ic++)
      for(int il=0;il<=l;il++)
	{
	  int t=ib+ic+il;
	  if(!filter(t)) continue;
       
	  // how many equivalent ways 
	  int totalCombinations=comb(b,ib)*comb(c,ic)*comb(l,il);
      
	  //       std::cout <<  totalCombinations << std::endl;
	  pMC+=1.*totalCombinations * pow(meanBeff,ib)*pow(1.-meanBeff,b-ib)*pow(meanCeff,ic)*pow(1.-meanCeff,c-ic)*pow(meanUDSGeff,il)*pow(1.-meanUDSGeff,l-il);
	  pData+=1.*totalCombinations *  pow(meanBeff*sfb,ib)*pow(1.-meanBeff*sfb,b-ib)*pow(meanCeff*sfc,ic)*pow(1.-meanCeff*sfc,c-ic)*pow(meanUDSGeff*sfl,il)*pow(1.-meanUDSGeff*sfl,l-il);
	  //	  std::cout << totalCombinations << " " << pow(meanBeff,ib)*pow(1.-meanBeff,b-ib)*pow(meanCeff,ic)*pow(1.-meanCeff,c-ic)*pow(meanUDSGeff,il)*pow(1.-meanUDSGeff,l-il) << " ";   
	  //	  std::cout <<  pow(meanBeff*sfb,ib)*pow(1.-meanBeff*sfb,b-ib)*pow(meanCeff*sfc,ic)*pow(1.-meanCeff*sfc,c-ic)*pow(meanUDSGeff*sfl,il)*pow(1.-meanUDSGeff*sfl,l-il) << std::endl ;
	  /* std::cout << pow(meanBeff*sfb,ib)*pow(1.-meanBeff*sfb,b-ib)*pow(meanCeff*sfc,ic)*pow(1.-meanCeff*sfc,c-ic)*pow(meanUDSGeff*sfl,il)*pow(1.-meanUDSGeff*sfl,l-il)/(pow(meanBeff,ib)*pow(1.-meanBeff,b-ib)*pow(meanCeff,ic)*pow(1.-meanCeff,c-ic)*pow(meanUDSGeff,il)*pow(1.-meanUDSGeff,l-il)) << std::endl; 
      std::cout << pMC << " " << pData <<  " " <<  pData/pMC << std::endl; 
	  */   
	} 
  if(pMC==0) return 0; 
  return pData/pMC;
}

int BTagWeight::fact(int n)
{
  if(n < 1) return 1;
  int r=1;
  for(int i=n;i > 1; i--)
    r*=i;
  
  return r;
}
int BTagWeight::comb(int n, int k)
{
  return fact(n)/fact(k)/fact(n-k);
}

bool BTagWeight::filter(int t)
{
  return (t >= minTags && t <= maxTags);
}



#endif

