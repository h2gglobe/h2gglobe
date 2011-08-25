#include <vector>
#include <map>
#include <iostream>

#include "TH1F.h"
#include "TGraph.h"
#include "TCanvas.h"

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
  XSectionMap[90]=33.8959;
  XSectionMap[95]=30.5228;
  XSectionMap[100]=27.5471;
  XSectionMap[105]=24.9582;
  XSectionMap[110]=22.7112;
  XSectionMap[110.5]=22.4929;
  XSectionMap[111]=22.2949;
  XSectionMap[111.5]=22.0882;
  XSectionMap[112]=21.891;
  XSectionMap[112.5]=21.6841;
  XSectionMap[113]=21.4987;
  XSectionMap[113.5]=21.3024;
  XSectionMap[114]=21.1075;
  XSectionMap[114.5]=20.9231;
  XSectionMap[115]=20.7379;
  XSectionMap[115.5]=20.5543;
  XSectionMap[116]=20.3799;
  XSectionMap[116.5]=20.197;
  XSectionMap[117]=20.0243;
  XSectionMap[117.5]=19.851;
  XSectionMap[118]=19.6789;
  XSectionMap[118.5]=19.5072;
  XSectionMap[119]=19.3457;
  XSectionMap[119.5]=19.18438;
  XSectionMap[120]=19.01246;
  XSectionMap[120.5]=18.85166;
  XSectionMap[121]=18.70108;
  XSectionMap[121.5]=18.54062;
  XSectionMap[122]=18.38147;
  XSectionMap[122.5]=18.22164;
  XSectionMap[123]=18.08183;
  XSectionMap[123.5]=17.93234;
  XSectionMap[124]=17.78306;
  XSectionMap[124.5]=17.63389;
  XSectionMap[125]=17.49604;
  XSectionMap[125.5]=17.3472;
  XSectionMap[126]=17.20858;
  XSectionMap[126.5]=17.07017;
  XSectionMap[127]=16.93297;
  XSectionMap[127.5]=16.79489;
  XSectionMap[128]=16.65702;
  XSectionMap[128.5]=16.52027;
  XSectionMap[129]=16.39272;
  XSectionMap[129.5]=16.2553;
  XSectionMap[130]=16.12918;
  XSectionMap[130.5]=16.00218;
  XSectionMap[131]=15.87639;
  XSectionMap[131.5]=15.75971;
  XSectionMap[132]=15.63414;
  XSectionMap[132.5]=15.50779;
  XSectionMap[133]=15.39265;
  XSectionMap[133.5]=15.26662;
  XSectionMap[134]=15.1617;
  XSectionMap[134.5]=15.047;
  XSectionMap[135]=14.9324;
  XSectionMap[135.5]=14.81791;
  XSectionMap[136]=14.70354;
  XSectionMap[136.5]=14.58927;
  XSectionMap[137]=14.47532;
  XSectionMap[137.5]=14.38228;
  XSectionMap[138]=14.26845;
  XSectionMap[138.5]=14.15482;
  XSectionMap[139]=14.05231;
  XSectionMap[139.5]=13.94891;
  XSectionMap[140]=13.84562;
  XSectionMap[141]=13.64057;
  XSectionMap[142]=13.44507;
  XSectionMap[143]=13.24;
  XSectionMap[144]=13.05526;
  XSectionMap[145]=12.86195;
  XSectionMap[146]=12.68876;
  XSectionMap[147]=12.5064;
  XSectionMap[148]=12.32437;
  XSectionMap[149]=12.15286;
  XSectionMap[150]=11.98179;
  XSectionMap[151]=11.82075;
  XSectionMap[152]=11.64993;
  XSectionMap[153]=11.48934;
  XSectionMap[154]=11.33298;
  XSectionMap[155]=11.17384;
  XSectionMap[156]=11.01263;
  XSectionMap[157]=10.85104;
  XSectionMap[158]=10.68418;
  XSectionMap[159]=10.52334;
  XSectionMap[160]=10.36062;
  XSectionMap[162]=10.02523;
  XSectionMap[164]=9.69762;
  XSectionMap[166]=9.39728;
  XSectionMap[168]=9.1301;
  XSectionMap[170]=8.87739;
  XSectionMap[172]=8.63504;
  XSectionMap[174]=8.40245;
  XSectionMap[176]=8.17764;
  XSectionMap[178]=7.96305;
  XSectionMap[180]=7.75479;
  XSectionMap[182]=7.54702;
  XSectionMap[184]=7.35255;
  XSectionMap[186]=7.16541;
  XSectionMap[188]=6.98324;
  XSectionMap[190]=6.80952;
  XSectionMap[192]=6.64487;
  XSectionMap[194]=6.49068;
  XSectionMap[196]=6.34284;
  XSectionMap[198]=6.20153;
  XSectionMap[200]=6.06875;
  XSectionMap[202]=5.93861;
  XSectionMap[204]=5.81288;
  XSectionMap[206]=5.69375;
  XSectionMap[208]=5.57598;
  XSectionMap[210]=5.46177;
  XSectionMap[212]=5.35397;
  XSectionMap[214]=5.24869;
  XSectionMap[216]=5.14975;
  XSectionMap[218]=5.05121;
  XSectionMap[220]=4.95707;
  XSectionMap[222]=4.86323;
  XSectionMap[224]=4.77178;
  XSectionMap[226]=4.68459;
  XSectionMap[228]=4.59966;
  XSectionMap[230]=4.51619;
  XSectionMap[232]=4.43594;
  XSectionMap[234]=4.35701;
  XSectionMap[236]=4.28231;
  XSectionMap[238]=4.20786;
  XSectionMap[240]=4.135713;
  XSectionMap[242]=4.067887;
  XSectionMap[244]=4.003181;
  XSectionMap[246]=3.939776;
  XSectionMap[248]=3.87848;
  XSectionMap[250]=3.819473;

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

  std::cout << "Warning cross section outside range of 90-250GeV!!!!" << std::endl;
  exit(1);

}

double GetNorm(double mass1, TH1F* hist1, double mass2, TH1F* hist2, double mass) {

  double br = GetBR(mass);
  double br1 = GetBR(mass1);
  double br2 = GetBR(mass2);
  
  double xsec = GetXsection(mass);
  double xsec1 = GetXsection(mass1);
  double xsec2 = GetXsection(mass2);
  
  double alpha = 1.0*(mass-mass1)/(mass2-mass1);
  double effAcc1 = hist1->Integral()/(xsec1*br1);
  double effAcc2 = hist2->Integral()/(xsec2*br2);

  double Normalization = (xsec*br)*(effAcc1 + alpha * (effAcc2 - effAcc1));
  return Normalization;
  
}

void CheckNorm(double Min, double Max, double Step) {

  vector <double> Mass;
  vector <double> BranchingRatio;
  vector <double> XSection;
  for (double i=Min; i<Max; i+=Step) {
    Mass.push_back(i);
    BranchingRatio.push_back(GetBR(i));
    XSection.push_back(GetXsection(i));
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
