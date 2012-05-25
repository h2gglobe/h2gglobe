#include <vector>
#include <map>
#include <iostream>

#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TString.h"

double GetBR(double mass) {

  map <double, double> BranchingRatioMap;
  BranchingRatioMap[90]=0.00123;
  BranchingRatioMap[95]=0.0014;
  BranchingRatioMap[100]=0.00159;
  BranchingRatioMap[105]=0.00178;
  BranchingRatioMap[110]=0.00197;
  BranchingRatioMap[110.5]=0.00199;
  BranchingRatioMap[111]=0.002;
  BranchingRatioMap[111.5]=0.00202;
  BranchingRatioMap[112]=0.00204;
  BranchingRatioMap[112.5]=0.00205;
  BranchingRatioMap[113]=0.00207;
  BranchingRatioMap[113.5]=0.00209;
  BranchingRatioMap[114]=0.0021;
  BranchingRatioMap[114.5]=0.00212;
  BranchingRatioMap[115]=0.00213;
  BranchingRatioMap[115.5]=0.00215;
  BranchingRatioMap[116]=0.00216;
  BranchingRatioMap[116.5]=0.00217;
  BranchingRatioMap[117]=0.00218;
  BranchingRatioMap[117.5]=0.0022;
  BranchingRatioMap[118]=0.00221;
  BranchingRatioMap[118.5]=0.00222;
  BranchingRatioMap[119]=0.00223;
  BranchingRatioMap[119.5]=0.00224;
  BranchingRatioMap[120]=0.00225;
  BranchingRatioMap[120.5]=0.00226;
  BranchingRatioMap[121]=0.00226;
  BranchingRatioMap[121.5]=0.00227;
  BranchingRatioMap[122]=0.00228;
  BranchingRatioMap[122.5]=0.00228;
  BranchingRatioMap[123]=0.00228;
  BranchingRatioMap[123.5]=0.00229;
  BranchingRatioMap[124]=0.00229;
  BranchingRatioMap[124.5]=0.00229;
  BranchingRatioMap[125]=0.00229;
  BranchingRatioMap[125.5]=0.00229;
  BranchingRatioMap[126]=0.00229;
  BranchingRatioMap[126.5]=0.00229;
  BranchingRatioMap[127]=0.00229;
  BranchingRatioMap[127.5]=0.00229;
  BranchingRatioMap[128]=0.00228;
  BranchingRatioMap[128.5]=0.00228;
  BranchingRatioMap[129]=0.00227;
  BranchingRatioMap[129.5]=0.00227;
  BranchingRatioMap[130]=0.00226;
  BranchingRatioMap[130.5]=0.00225;
  BranchingRatioMap[131]=0.00224;
  BranchingRatioMap[131.5]=0.00223;
  BranchingRatioMap[132]=0.00222;
  BranchingRatioMap[132.5]=0.00221;
  BranchingRatioMap[133]=0.00219;
  BranchingRatioMap[133.5]=0.00218;
  BranchingRatioMap[134]=0.00217;
  BranchingRatioMap[134.5]=0.00215;
  BranchingRatioMap[135]=0.00213;
  BranchingRatioMap[135.5]=0.00212;
  BranchingRatioMap[136]=0.0021;
  BranchingRatioMap[136.5]=0.00208;
  BranchingRatioMap[137]=0.00206;
  BranchingRatioMap[137.5]=0.00204;
  BranchingRatioMap[138]=0.00202;
  BranchingRatioMap[138.5]=0.002;
  BranchingRatioMap[139]=0.00198;
  BranchingRatioMap[139.5]=0.00196;
  BranchingRatioMap[140]=0.00193;
  BranchingRatioMap[141]=0.00189;
  BranchingRatioMap[142]=0.00184;
  BranchingRatioMap[143]=0.00178;
  BranchingRatioMap[144]=0.00173;
  BranchingRatioMap[145]=0.00167;
  BranchingRatioMap[146]=0.00162;
  BranchingRatioMap[147]=0.00156;
  BranchingRatioMap[148]=0.00149;
  BranchingRatioMap[149]=0.00143;
  BranchingRatioMap[150]=0.00136;
  BranchingRatioMap[151]=0.0013;
  BranchingRatioMap[152]=0.00123;
  BranchingRatioMap[153]=0.00115;
  BranchingRatioMap[154]=0.00108;
  BranchingRatioMap[155]=0.000999;
  BranchingRatioMap[156]=0.000914;
  BranchingRatioMap[157]=0.000825;
  BranchingRatioMap[158]=0.000729;
  BranchingRatioMap[159]=0.000628;
  BranchingRatioMap[160]=0.000532;
  BranchingRatioMap[162]=0.00037;
  BranchingRatioMap[164]=0.000259;
  BranchingRatioMap[166]=0.000208;
  BranchingRatioMap[168]=0.000178;
  BranchingRatioMap[170]=0.000158;
  BranchingRatioMap[172]=0.000143;
  BranchingRatioMap[174]=0.000132;
  BranchingRatioMap[176]=0.000122;
  BranchingRatioMap[178]=0.000113;
  BranchingRatioMap[180]=0.000105;
  BranchingRatioMap[182]=0.0000968;
  BranchingRatioMap[184]=0.0000881;
  BranchingRatioMap[186]=0.0000809;
  BranchingRatioMap[188]=0.0000752;
  BranchingRatioMap[190]=0.0000705;
  BranchingRatioMap[192]=0.0000666;
  BranchingRatioMap[194]=0.0000632;
  BranchingRatioMap[196]=0.0000602;
  BranchingRatioMap[198]=0.0000575;
  BranchingRatioMap[200]=0.0000551;
  BranchingRatioMap[202]=0.0000529;
  BranchingRatioMap[204]=0.0000508;
  BranchingRatioMap[206]=0.0000489;
  BranchingRatioMap[208]=0.0000471;
  BranchingRatioMap[210]=0.0000454;
  BranchingRatioMap[212]=0.0000439;
  BranchingRatioMap[214]=0.0000424;
  BranchingRatioMap[216]=0.000041;
  BranchingRatioMap[218]=0.0000396;
  BranchingRatioMap[220]=0.0000384;
  BranchingRatioMap[222]=0.0000372;
  BranchingRatioMap[224]=0.000036;
  BranchingRatioMap[226]=0.0000349;
  BranchingRatioMap[228]=0.0000339;
  BranchingRatioMap[230]=0.0000328;
  BranchingRatioMap[232]=0.0000319;
  BranchingRatioMap[234]=0.0000309;
  BranchingRatioMap[236]=0.0000301;
  BranchingRatioMap[238]=0.0000292;
  BranchingRatioMap[240]=0.0000284;
  BranchingRatioMap[242]=0.0000276;
  BranchingRatioMap[244]=0.0000268;
  BranchingRatioMap[246]=0.0000261;
  BranchingRatioMap[248]=0.0000254;
  BranchingRatioMap[250]=0.0000247;

  for (std::map<double, double>::const_iterator iter = BranchingRatioMap.begin();  iter != BranchingRatioMap.end(); ++iter) {
    if (mass==iter->first) return iter->second;
    if (mass>iter->first) {
      double lowmass = iter->first;
      double lowbr = iter->second;
      ++iter;
      if (mass<iter->first) {
        double highmass = iter->first;
        double highbr = iter->second;
        double br = (highbr-lowbr)/(highmass-lowmass)*(mass-lowmass)+lowbr;
        return br;
      }
      --iter;
    }
  }
  
  std::cout << "Warning branching ratio outside range of 90-250GeV!!!!" << std::endl;
  exit(1);
  
}

double GetXsection(double mass) {

  map <double, double> XSectionMap;
  XSectionMap[80]=46.0593;
  XSectionMap[81]=44.9777;
  XSectionMap[82]=43.9344;
  XSectionMap[83]=42.928; 
  XSectionMap[84]=41.9566;
  XSectionMap[85]=41.0179;
  XSectionMap[86]=40.1112;
  XSectionMap[87]=39.2347;
  XSectionMap[88]=38.3868;
  XSectionMap[89]=37.5664;
  XSectionMap[90]=36.7719;
  XSectionMap[91]=36.0024;
  XSectionMap[92]=35.2565;
  XSectionMap[93]=34.5342;
  XSectionMap[94]=33.833; 
  XSectionMap[95]=33.1525;
  XSectionMap[96]=32.4924;
  XSectionMap[97]=31.8511;
  XSectionMap[98]=31.2283;
  XSectionMap[99]=30.6236;
  XSectionMap[100]=30.0891;
  XSectionMap[101]=29.5154;
  XSectionMap[102]=28.9585;
  XSectionMap[103]=28.4107;
  XSectionMap[104]=27.8856;
  XSectionMap[105]=27.3602;
  XSectionMap[106]=26.8628;
  XSectionMap[107]=26.3911;
  XSectionMap[108]=25.9199;
  XSectionMap[109]=25.4603;
  XSectionMap[110]=25.0119;
  XSectionMap[110.5]=24.7918;
  XSectionMap[111]=24.5743;
  XSectionMap[111.5]=24.3597;
  XSectionMap[112]=24.1486;
  XSectionMap[112.5]=23.9401;
  XSectionMap[113]=23.7343;
  XSectionMap[113.5]=23.5312;
  XSectionMap[114]=23.3306;
  XSectionMap[114.5]=23.1325;
  XSectionMap[115]=22.9369;
  XSectionMap[115.5]=22.7442;
  XSectionMap[116]=22.5534;
  XSectionMap[116.5]=22.3651;
  XSectionMap[117]=22.179; 
  XSectionMap[117.5]=21.9949;
  XSectionMap[118]=21.8134;
  XSectionMap[118.5]=21.6341;
  XSectionMap[119]=21.457; 
  XSectionMap[119.5]=21.282; 
  XSectionMap[120]=21.109; 
  XSectionMap[120.5]=20.9381;
  XSectionMap[121]=20.7692;
  XSectionMap[121.5]=20.6023;
  XSectionMap[122]=20.4373;
  XSectionMap[122.5]=20.2742;
  XSectionMap[123]=20.1131;
  XSectionMap[123.5]=19.9538;
  XSectionMap[124]=19.7963;
  XSectionMap[124.5]=19.6407;
  XSectionMap[125]=19.4868;
  XSectionMap[125.5]=19.3347;
  XSectionMap[126]=19.1843;
  XSectionMap[126.5]=19.0357;
  XSectionMap[127]=18.8887;
  XSectionMap[127.5]=18.7433;
  XSectionMap[128]=18.5996;
  XSectionMap[128.5]=18.4574;
  XSectionMap[129]=18.3169;
  XSectionMap[129.5]=18.1778;
  XSectionMap[130]=18.0403;
  XSectionMap[130.5]=17.9043;
  XSectionMap[131]=17.7698;
  XSectionMap[131.5]=17.6368;
  XSectionMap[132]=17.5053;
  XSectionMap[132.5]=17.3752;
  XSectionMap[133]=17.2464;
  XSectionMap[133.5]=17.1191;
  XSectionMap[134]=16.993; 
  XSectionMap[134.5]=16.8683;
  XSectionMap[135]=16.7448;
  XSectionMap[135.5]=16.6226;
  XSectionMap[136]=16.5017;
  XSectionMap[136.5]=16.382; 
  XSectionMap[137]=16.2635;
  XSectionMap[137.5]=16.1463;
  XSectionMap[138]=16.0303;
  XSectionMap[138.5]=15.9155;
  XSectionMap[139]=15.8019;
  XSectionMap[139.5]=15.6894;
  XSectionMap[140]=15.5782;
  XSectionMap[141]=15.3592;
  XSectionMap[142]=15.1446;
  XSectionMap[143]=14.9342;
  XSectionMap[144]=14.7278;
  XSectionMap[145]=14.5253;
  XSectionMap[146]=14.3267;
  XSectionMap[147]=14.1318;
  XSectionMap[148]=13.9404;
  XSectionMap[149]=13.7522;
  XSectionMap[150]=13.5672;
  XSectionMap[151]=13.3852;
  XSectionMap[152]=13.2057;
  XSectionMap[153]=13.0284;
  XSectionMap[154]=12.8527;
  XSectionMap[155]=12.6779;
  XSectionMap[156]=12.5028;
  XSectionMap[157]=12.3256;
  XSectionMap[158]=12.1439;
  XSectionMap[159]=11.9554;
  XSectionMap[160]=11.7624;
  XSectionMap[162]=11.3963;
  XSectionMap[164]=11.0387;
  XSectionMap[166]=10.7032;
  XSectionMap[168]=10.3946;
  XSectionMap[170]=10.1077;
  XSectionMap[172]=9.8376;
  XSectionMap[174]=9.5807;
  XSectionMap[176]=9.3337;
  XSectionMap[178]=9.0926;
  XSectionMap[180]=8.8521;
  XSectionMap[182]=8.6143;
  XSectionMap[184]=8.3939;
  XSectionMap[186]=8.179; 
  XSectionMap[188]=7.9739;
  XSectionMap[190]=7.7803;
  XSectionMap[192]=7.5975;
  XSectionMap[194]=7.4242;
  XSectionMap[196]=7.259; 
  XSectionMap[198]=7.101; 
  XSectionMap[200]=6.9497;
  XSectionMap[202]=6.8042;
  XSectionMap[204]=6.6642;
  XSectionMap[206]=6.5296;
  XSectionMap[208]=6.3999;
  XSectionMap[210]=6.2744;
  XSectionMap[212]=6.1533;
  XSectionMap[214]=6.0364;
  XSectionMap[216]=5.9234;
  XSectionMap[218]=5.8142;
  XSectionMap[220]=5.7085;
  XSectionMap[222]=5.6061;
  XSectionMap[224]=5.5071;
  XSectionMap[226]=5.4112;
  XSectionMap[228]=5.3183;
  XSectionMap[230]=5.2283;
  XSectionMap[232]=5.141; 
  XSectionMap[234]=5.0563;
  XSectionMap[236]=4.974; 
  XSectionMap[238]=4.8944;
  XSectionMap[240]=4.817; 
  XSectionMap[242]=4.7418;
  XSectionMap[244]=4.6689;
  XSectionMap[246]=4.5981;
  XSectionMap[248]=4.5294;
  XSectionMap[250]=4.4627;
  XSectionMap[252]=4.398; 
  XSectionMap[254]=4.3351;
  XSectionMap[256]=4.2741;
  XSectionMap[258]=4.2149;
  XSectionMap[260]=4.1574;
  XSectionMap[262]=4.1016;
  XSectionMap[264]=4.0709;
  XSectionMap[266]=3.9951;
  XSectionMap[268]=3.9442;
  XSectionMap[270]=3.8949;
  XSectionMap[272]=3.8471;
  XSectionMap[274]=3.8009;
  XSectionMap[276]=3.7561;
  XSectionMap[278]=3.7128;
  XSectionMap[280]=3.6709;
  XSectionMap[282]=3.6302;
  XSectionMap[284]=3.591; 
  XSectionMap[286]=3.553; 
  XSectionMap[288]=3.5164;
  XSectionMap[290]=3.4811;
  XSectionMap[295]=3.3986;
  XSectionMap[300]=3.324; 
  
  for (std::map<double, double>::const_iterator iter = XSectionMap.begin();  iter != XSectionMap.end(); ++iter) {
    if (mass==iter->first) return iter->second;
    if (mass>iter->first) {
      double lowmass = iter->first;
      double lowxsec = iter->second;
      ++iter;
      if (mass<iter->first) {
        double highmass = iter->first;
        double highxsec = iter->second;
        double xsec = (highxsec-lowxsec)/(highmass-lowmass)*(mass-lowmass)+lowxsec;
        return xsec;
      }
      --iter;
    }
  }

  std::cout << "Warning cross section outside range of 80-300GeV!!!!" << std::endl;
  exit(1);

}

double GetXsection(double mass, TString HistName) {

  map <double, double> XSectionMap;
  if (HistName.Contains("ggh")) {
    XSectionMap[80]=46.0593;
    XSectionMap[81]=44.9777;
    XSectionMap[82]=43.9344;
    XSectionMap[83]=42.928; 
    XSectionMap[84]=41.9566;
    XSectionMap[85]=41.0179;
    XSectionMap[86]=40.1112;
    XSectionMap[87]=39.2347;
    XSectionMap[88]=38.3868;
    XSectionMap[89]=37.5664;
    XSectionMap[90]=36.7719;
    XSectionMap[91]=36.0024;
    XSectionMap[92]=35.2565;
    XSectionMap[93]=34.5342;
    XSectionMap[94]=33.833; 
    XSectionMap[95]=33.1525;
    XSectionMap[96]=32.4924;
    XSectionMap[97]=31.8511;
    XSectionMap[98]=31.2283;
    XSectionMap[99]=30.6236;
    XSectionMap[100]=30.0891;
    XSectionMap[101]=29.5154;
    XSectionMap[102]=28.9585;
    XSectionMap[103]=28.4107;
    XSectionMap[104]=27.8856;
    XSectionMap[105]=27.3602;
    XSectionMap[106]=26.8628;
    XSectionMap[107]=26.3911;
    XSectionMap[108]=25.9199;
    XSectionMap[109]=25.4603;
    XSectionMap[110]=25.0119;
    XSectionMap[110.5]=24.7918;
    XSectionMap[111]=24.5743;
    XSectionMap[111.5]=24.3597;
    XSectionMap[112]=24.1486;
    XSectionMap[112.5]=23.9401;
    XSectionMap[113]=23.7343;
    XSectionMap[113.5]=23.5312;
    XSectionMap[114]=23.3306;
    XSectionMap[114.5]=23.1325;
    XSectionMap[115]=22.9369;
    XSectionMap[115.5]=22.7442;
    XSectionMap[116]=22.5534;
    XSectionMap[116.5]=22.3651;
    XSectionMap[117]=22.179; 
    XSectionMap[117.5]=21.9949;
    XSectionMap[118]=21.8134;
    XSectionMap[118.5]=21.6341;
    XSectionMap[119]=21.457; 
    XSectionMap[119.5]=21.282; 
    XSectionMap[120]=21.109; 
    XSectionMap[120.5]=20.9381;
    XSectionMap[121]=20.7692;
    XSectionMap[121.5]=20.6023;
    XSectionMap[122]=20.4373;
    XSectionMap[122.5]=20.2742;
    XSectionMap[123]=20.1131;
    XSectionMap[123.5]=19.9538;
    XSectionMap[124]=19.7963;
    XSectionMap[124.5]=19.6407;
    XSectionMap[125]=19.4868;
    XSectionMap[125.5]=19.3347;
    XSectionMap[126]=19.1843;
    XSectionMap[126.5]=19.0357;
    XSectionMap[127]=18.8887;
    XSectionMap[127.5]=18.7433;
    XSectionMap[128]=18.5996;
    XSectionMap[128.5]=18.4574;
    XSectionMap[129]=18.3169;
    XSectionMap[129.5]=18.1778;
    XSectionMap[130]=18.0403;
    XSectionMap[130.5]=17.9043;
    XSectionMap[131]=17.7698;
    XSectionMap[131.5]=17.6368;
    XSectionMap[132]=17.5053;
    XSectionMap[132.5]=17.3752;
    XSectionMap[133]=17.2464;
    XSectionMap[133.5]=17.1191;
    XSectionMap[134]=16.993; 
    XSectionMap[134.5]=16.8683;
    XSectionMap[135]=16.7448;
    XSectionMap[135.5]=16.6226;
    XSectionMap[136]=16.5017;
    XSectionMap[136.5]=16.382; 
    XSectionMap[137]=16.2635;
    XSectionMap[137.5]=16.1463;
    XSectionMap[138]=16.0303;
    XSectionMap[138.5]=15.9155;
    XSectionMap[139]=15.8019;
    XSectionMap[139.5]=15.6894;
    XSectionMap[140]=15.5782;
    XSectionMap[141]=15.3592;
    XSectionMap[142]=15.1446;
    XSectionMap[143]=14.9342;
    XSectionMap[144]=14.7278;
    XSectionMap[145]=14.5253;
    XSectionMap[146]=14.3267;
    XSectionMap[147]=14.1318;
    XSectionMap[148]=13.9404;
    XSectionMap[149]=13.7522;
    XSectionMap[150]=13.5672;
    XSectionMap[151]=13.3852;
    XSectionMap[152]=13.2057;
    XSectionMap[153]=13.0284;
    XSectionMap[154]=12.8527;
    XSectionMap[155]=12.6779;
    XSectionMap[156]=12.5028;
    XSectionMap[157]=12.3256;
    XSectionMap[158]=12.1439;
    XSectionMap[159]=11.9554;
    XSectionMap[160]=11.7624;
    XSectionMap[162]=11.3963;
    XSectionMap[164]=11.0387;
    XSectionMap[166]=10.7032;
    XSectionMap[168]=10.3946;
    XSectionMap[170]=10.1077;
    XSectionMap[172]=9.8376;
    XSectionMap[174]=9.5807;
    XSectionMap[176]=9.3337;
    XSectionMap[178]=9.0926;
    XSectionMap[180]=8.8521;
    XSectionMap[182]=8.6143;
    XSectionMap[184]=8.3939;
    XSectionMap[186]=8.179; 
    XSectionMap[188]=7.9739;
    XSectionMap[190]=7.7803;
    XSectionMap[192]=7.5975;
    XSectionMap[194]=7.4242;
    XSectionMap[196]=7.259; 
    XSectionMap[198]=7.101; 
    XSectionMap[200]=6.9497;
    XSectionMap[202]=6.8042;
    XSectionMap[204]=6.6642;
    XSectionMap[206]=6.5296;
    XSectionMap[208]=6.3999;
    XSectionMap[210]=6.2744;
    XSectionMap[212]=6.1533;
    XSectionMap[214]=6.0364;
    XSectionMap[216]=5.9234;
    XSectionMap[218]=5.8142;
    XSectionMap[220]=5.7085;
    XSectionMap[222]=5.6061;
    XSectionMap[224]=5.5071;
    XSectionMap[226]=5.4112;
    XSectionMap[228]=5.3183;
    XSectionMap[230]=5.2283;
    XSectionMap[232]=5.141; 
    XSectionMap[234]=5.0563;
    XSectionMap[236]=4.974; 
    XSectionMap[238]=4.8944;
    XSectionMap[240]=4.817; 
    XSectionMap[242]=4.7418;
    XSectionMap[244]=4.6689;
    XSectionMap[246]=4.5981;
    XSectionMap[248]=4.5294;
    XSectionMap[250]=4.4627;
    XSectionMap[252]=4.398; 
    XSectionMap[254]=4.3351;
    XSectionMap[256]=4.2741;
    XSectionMap[258]=4.2149;
    XSectionMap[260]=4.1574;
    XSectionMap[262]=4.1016;
    XSectionMap[264]=4.0709;
    XSectionMap[266]=3.9951;
    XSectionMap[268]=3.9442;
    XSectionMap[270]=3.8949;
    XSectionMap[272]=3.8471;
    XSectionMap[274]=3.8009;
    XSectionMap[276]=3.7561;
    XSectionMap[278]=3.7128;
    XSectionMap[280]=3.6709;
    XSectionMap[282]=3.6302;
    XSectionMap[284]=3.591; 
    XSectionMap[286]=3.553; 
    XSectionMap[288]=3.5164;
    XSectionMap[290]=3.4811;
    XSectionMap[295]=3.3986;
    XSectionMap[300]=3.324;
  } else {
    std::cout << "Warning ggh not found in histname!!!!" << std::endl;
    exit(1);
  }

  for (std::map<double, double>::const_iterator iter = XSectionMap.begin();  iter != XSectionMap.end(); ++iter) {
    if (mass==iter->first) return iter->second;
    if (mass>iter->first) {
      double lowmass = iter->first;
      double lowxsec = iter->second;
      ++iter;
      if (mass<iter->first) {
        double highmass = iter->first;
        double highxsec = iter->second;
        double xsec = (highxsec-lowxsec)/(highmass-lowmass)*(mass-lowmass)+lowxsec;
        return xsec;
      }
      --iter;
    }
  }

  std::cout << "Warning cross section outside range of 80-300GeV!!!!" << std::endl;
  exit(1);

}

double GetNorm(double mass1, TH1F* hist1, double mass2, TH1F* hist2, double mass) {

  double br = GetBR(mass);
  double br1 = GetBR(mass1);
  double br2 = GetBR(mass2);
  
  double xsec = GetXsection(mass, hist1->GetName());
  double xsec1 = GetXsection(mass1, hist1->GetName());
  double xsec2 = GetXsection(mass2, hist2->GetName());
  
  double alpha = 1.0*(mass-mass1)/(mass2-mass1);
  double effAcc1 = hist1->Integral()/(xsec1*br1);
  double effAcc2 = hist2->Integral()/(xsec2*br2);

  double Normalization = (xsec*br)*(effAcc1 + alpha * (effAcc2 - effAcc1));
  return Normalization;
  
}

void CheckNorm(double Min, double Max, double Step, TString histname="") {

  vector <double> Mass;
  vector <double> BranchingRatio;
  vector <double> XSection;
  for (double i=Min; i<Max; i+=Step) {
    Mass.push_back(i);
    BranchingRatio.push_back(GetBR(i));
    if (histname=="") XSection.push_back(GetXsection(i));
    else XSection.push_back(GetXsection(i,histname));
  }

  TGraph* BranchGraph = new TGraph(Mass.size(),&Mass[0],&BranchingRatio[0]);
  TGraph* XSectionGraph = new TGraph(Mass.size(),&Mass[0],&XSection[0]);
  BranchGraph->SetTitle("Interpolated Branching Ratios");
  XSectionGraph->SetTitle("Interpolated Cross Sections");
  BranchGraph->SetMarkerStyle(20);
  XSectionGraph->SetMarkerStyle(20);
  BranchGraph->SetMarkerSize(1);
  XSectionGraph->SetMarkerSize(1);

  TCanvas* c1 = new TCanvas("c1","c1",800,650);
  c1->cd();
  BranchGraph->Draw("AP");
  c1->SaveAs("BranchingRatios.png");
  c1->Clear();
  XSectionGraph->Draw("AP");
  c1->SaveAs("XSections.png");

  delete BranchGraph;
  delete XSectionGraph;
  delete c1;

}
