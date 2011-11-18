~capalmer/CMS/Higgs2gg/Limits/nov14_forAN/428p7/src/HiggsAnalysis/HiggsTo2photons/h2gglobe_new/PhotonAnalysis_scripts/histograms_CMS-HGG_4686_nov17_tryExcl.root


# FOR LOOPINTERACTIVE:

root -l
.L Marco/tdrstyle.C
setTDRStyle();
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  //gSystem->Load("VertexAnalysis/lib/libh2gglobeVertexAnalysis.so");
gSystem->Load("libLoopAll.so");
LoopAll* m=new LoopAll();
m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_CMS-HGG_4686_nov17_tryExcl.root","Hgg");
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_CMS-HGG_4686_nov17_tryExcl.root");



m->myPlotInteractiveSetup("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root","Hgg");
m->myPlotInteractive("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root");



m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_CMS-HGG_4686_nov17_tryExcl.root","Hgg");
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_CMS-HGG_4686_nov17_tryExcl.root");


m->myPlotInteractiveSetup("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root","Hgg");
m->myPlotInteractive("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root");

m->myPlotInteractiveSetup("PhotonAnalysis_scripts/histograms_CMS-HGG_4686_nov17_tryExcl.root","Hgg");
m->myPlotInteractive("PhotonAnalysis_scripts/histograms_CMS-HGG_4686_nov17_tryExcl.root");


m->myPlotInteractiveSetup("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root","Hgg");
m->myPlotInteractive("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root");




Instructions 16/11/2011
-----------------------

cd CMSSW_4_2_8/src
cvs co -r regression_Sept30 HiggsAnalysis/HiggsToGammaGamma
cvs co -r branch_for_42X -d HiggsAnalysis/HiggsTo2photons UserCode/HiggsAnalysis/HiggsTo2photons 

rm -r  HiggsAnalysis/HiggsTo2photons/h2gglobe

#check out the head of h2gglobe (we should make a tag

cvs co -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe

cp Marco/HistoContainer.cc ./HistoContainer.cc
cp Marco/HistoContainer.h ./HistoContainer.h

cp Marco/CounterContainer.cc ./CounterContainer.cc
cp Marco/CounterContainer.h ./CounterContainer.h

cp Marco/minimal_statanalysis_input.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
cp Marco/cuts.dat PhotonAnalysis_scripts/cuts.dat
cp Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
cp Marco/inputfiles.dat PhotonAnalysis_scripts/inputfiles.dat
# cp Marco/looper.py PhotonAnalysis_scripts/looper.py
# cp Marco/looper_input.dat PhotonAnalysis_scripts/looper_input.dat
cp Marco/photonanalysis.dat PhotonAnalysis_scripts/photonanalysis.dat
cp Marco/reduction_output.dat PhotonAnalysis_scripts/reduction_output.dat

cp Marco/configProducer.py PhotonAnalysis_scripts/python/configProducer.py
# cp Marco/PhotonAnalysis.cc PhotonAnalysis/src/PhotonAnalysis.cc

cp Marco/datafiles_5fb.dat PhotonAnalysis_scripts/datafiles_5fb.dat
cp Marco/statanalysis.dat PhotonAnalysis_scripts/statanalysis.dat

cp Marco/StatAnalysisExclusive.h PhotonAnalysis/interface/StatAnalysisExclusive.h
cp Marco/StatAnalysisExclusive.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
cp Marco/statanalysisexclusive.dat PhotonAnalysis_scripts/statanalysisexclusive.dat
cp Marco/fitter.py PhotonAnalysis_scripts/fitter.py

cp Marco/Makefile ./Makefile




diff Marco/HistoContainer.cc ./HistoContainer.cc
diff Marco/HistoContainer.h ./HistoContainer.h

diff Marco/CounterContainer.cc ./CounterContainer.cc
diff Marco/CounterContainer.h ./CounterContainer.h

diff Marco/minimal_statanalysis_input.dat PhotonAnalysis_scripts/minimal_statanalysis_input.dat
diff Marco/cuts.dat PhotonAnalysis_scripts/cuts.dat
diff Marco/plotvariables.dat PhotonAnalysis_scripts/plotvariables.dat
diff Marco/inputfiles.dat PhotonAnalysis_scripts/inputfiles.dat
# diff Marco/looper.py PhotonAnalysis_scripts/looper.py
# diff Marco/looper_input.dat PhotonAnalysis_scripts/looper_input.dat
diff Marco/photonanalysis.dat PhotonAnalysis_scripts/photonanalysis.dat
diff Marco/reduction_output.dat PhotonAnalysis_scripts/reduction_output.dat

diff Marco/configProducer.py PhotonAnalysis_scripts/python/configProducer.py
# diff Marco/PhotonAnalysis.cc PhotonAnalysis/src/PhotonAnalysis.cc

diff Marco/datafiles_5fb.dat PhotonAnalysis_scripts/datafiles_5fb.dat
diff Marco/statanalysis.dat PhotonAnalysis_scripts/statanalysis.dat

diff Marco/StatAnalysisExclusive.h PhotonAnalysis/interface/StatAnalysisExclusive.h
diff Marco/StatAnalysisExclusive.cc PhotonAnalysis/src/StatAnalysisExclusive.cc
diff Marco/statanalysisexclusive.dat PhotonAnalysis_scripts/statanalysisexclusive.dat
diff Marco/fitter.py PhotonAnalysis_scripts/fitter.py

diff Marco/Makefile ./Makefile






uncomment the second to last line in CommonParameters.h

Modify in PhotonAnalysis/interface/PhotonAnalysis.h the line similar into:
#include "../../../../HiggsToGammaGamma/interface/GBRForest.h"

Then make

make clean; make -j 30

cd PhotonAnalysis_scripts/.
python fitter.py -i datafiles_5fb.dat --dryRun
python fitter.py -i datafiles_5fb.dat 

Plotinteractive and myprintcountersnew are still to be fixed

When finished if you want co commit everything:

rm CommonParameters.h 
rm HistoContainer.cc HistoContainer.h PhotonAnalysis_scripts/cuts.dat PhotonAnalysis_scripts/plotvariables.dat PhotonAnalysis_scripts/inputfiles.dat PhotonAnalysis_scripts/python/configProducer.py  PhotonAnalysis/src/PhotonAnalysis.cc PhotonAnalysis_scripts/looper.py PhotonAnalysis_scripts/looper_input.dat PhotonAnalysis_scripts/photonanalysis.dat PhotonAnalysis_scripts/reduction_output.dat PhotonAnalysis_scripts/datafiles_5fb.dat PhotonAnalysis_scripts/statanalysis.dat
rm PhotonAnalysis/interface/StatAnalysisExclusive.h PhotonAnalysis/src/StatAnalysisExclusive.cc Makefile PhotonAnalysis_scripts/statanalysisexclusive.dat PhotonAnalysis_scripts/fitter.py
rm CounterContainer.cc CounterContainer.h PhotonAnalysis_scripts/minimal_statanalysis_input.dat

cvs update -A

=====================================
=====================================
=====================================
=====================================
END FOR NOW
=====================================
=====================================
=====================================
=====================================




UCSD:
-----

cd ~/Aug2011/mpieri

source /home/users/mpieri/cmsset_default.sh
export SCRAM_ARCH="slc5_amd64_gcc434"

cmsenv
export CVS_RSH=ssh
export CVSROOT=:ext:mpieri@cmscvs.cern.ch:/cvs_server/repositories/CMSSW

#setenv CVS_RSH ssh
#setenv CVSROOT :ext:mpieri@cmscvs.cern.ch:/cvs_server/repositories/CMSSW

kinit mpieri@CERN.CH

scramv1 project CMSSW CMSSW_4_2_8

cd CMSSW_4_2_8/src
cvs co -r regression_Sept30 HiggsAnalysis/HiggsToGammaGamma
cvs co -r branch_for_42X -d HiggsAnalysis/HiggsTo2photons UserCode/HiggsAnalysis/HiggsTo2photons 

rm -r  HiggsAnalysis/HiggsTo2photons/h2gglobe

#cvs co -r nw_16_11_11 -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe

cvs co -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe

cmsenv
scramv1 b

cd HiggsAnalysis/HiggsTo2photons/h2gglobe
make clean; make



Hi All,

I made a tag for the h2gglobe package I used for producing the Workspace from which we obtained the limit plot yesterday (nw_16_11_11 ).
Below is a recipe for making the workspace at LXPLUS using the 4.69/fb and the mix of Summer11 and Fall11 MC since some people have asked
how to produce it.

cmsrel CMSSW_4_2_8
cd CMSSW_4_2_8/src

cvs co -r nw_16_11_11 -dHiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe
cvs co HiggsAnalysis/HiggsToGammaGamma (HEAD is ok, but maybe Chris Palmer can comment on what he used for Ntuple production)

scramv1 b -j8
cd HiggsAnalysis/HiggsTo2photons/h2gglobe

cmsenv
make clean ; make -j8

cd PhotonAnalysis_scripts
python fitter.py -i datafiles_5fb.dat --dryRun

the last command will configure to make sure all the files are there for reweighting etc and the ntuples themselves exist.
Note this also create a file called datafiles_5fb.dat.pevents which you should not alter. If you change the file datafiles_5fb.dat, then you should
delete the datafiles_5fb.dat.pevents and rerun the above command.

next, you can actually run the jobs by the following:

python fitter.py -i datafiles_5fb.dat -n 30 -j 0

which runs the first of 30 jobs (its up to you how many you submit but i find 30 completes in a reasonable time)
You can submit each of the jobs individually to the batch machine at LXPLUS if you know how. I put a little script here
~nckw/pubic/makesub.py
which will create submission scripts for each job, just do

python makesub.py 30 $PWD datafiles_5fb.dat

and then you can submit each one to batch with the usual
qsub -q 1nh sub0.sh  etc...
1nh should be long enough but check your logs in case not.

You will now have the 30 subjobs (CMS-HGG_4686pb_(0-29).root ) which you should move into a folder called cms-hgg-workspaces,
then finally

python combiner.py -i filestocombine.dat

will merge those files into one called CMS-HGG_4686pb.root

This is the workspace which goes to the binned limit setting for Generated MC points.
There are Macros in h2gglobe/Macros which do the interpolation and diagnostics but I haven't the time to explain them all here. Please contact me if you want to take the binned limit setting through. Otherwise, from that workspace, Josh has a script to overlay a signal model to do the unbinned one.

Let me know if anything is unclear or you run into any problems

Cheers,
Nick


  hfile->cd();

//int dummy = fscanf(file,"%d plot=%d idummy=%d",&Nvar, &typplotall, &idummystart);
//h2d[i], typplot[i], histoncat[i], nbinsx[i], nbinsy[i], lowlim[i],highlim[i],lowlim2[i],highlim2[i],varnamescread[i]
  plotvartree = new TTree("plotvariables","globe plotvariables provenance information");
  plotvartree->Branch("Nvar", &Nvar, "Nvar/I");
  plotvartree->Branch("typplotall", &typplotall, "typplotall/I");
  plotvartree->Branch("doplot", &plothistoplotitPI, "plothistoplotitPI[Nvar]/I");
  plotvartree->Branch("h2d", &h2d, "h2d[Nvar]/I");
  plotvartree->Branch("typplot", &typplot, "typplot[Nvar]/I");
  plotvartree->Branch("histoncat", &histoncat, "histoncat[Nvar]/I");
  plotvartree->Branch("histoncatindtonames", &histoncatindtonames, "histoncatindtonames[Nvar]/I");
  plotvartree->Branch("nbinsx", &nbinsx, "nbinsx[Nvar]/I");
  plotvartree->Branch("nbinsy", &nbinsy, "nbinsy[Nvar]/I");
  plotvartree->Branch("lowlim", &lowlim, "lowlim[Nvar]/F");
  plotvartree->Branch("highlim", &highlim, "highlim[Nvar]/F");
  plotvartree->Branch("lowlim2", &lowlim2, "lowlim2[Nvar]/F");
  plotvartree->Branch("highlim2", &highlim2, "highlim2[Nvar]/F");

  tca_xaxislabels = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("xaxislabels", "TClonesArray", &tca_xaxislabels, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_xaxislabels)[iplot]) TObjString(); 
    ((TObjString *)tca_xaxislabels->At(iplot))->SetString(xaxislabel[iplot]);
  }

  tca_yaxislabels = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("yaxislabels", "TClonesArray", &tca_yaxislabels, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_yaxislabels)[iplot]) TObjString(); 
    ((TObjString *)tca_yaxislabels->At(iplot))->SetString(yaxislabel[iplot]);
  }

  tca_plotvarnames = new TClonesArray("TObjString",Nvar);
  plotvartree->Branch("plotvarnames", "TClonesArray", &tca_plotvarnames, 32000, 0);
  for(int iplot=0;iplot!=Nvar;++iplot) { 
    new ((*tca_plotvarnames)[iplot]) TObjString(); 
    ((TObjString *)tca_plotvarnames->At(iplot))->SetString(varnamescread[iplot]);
  }

  plotvartree->Branch("Nvarcats", &Nvarcats, "Nvarcats/I");
  plotvartree->Branch("catid", &catid, "catid[Nvarcats]/I");
  plotvartree->Branch("ncats", &ncats, "ncats[Nvarcats]/I");
  tca_plotvarcatnames = new TClonesArray("TObjString",Ncatvar);
  plotvartree->Branch("plotvarcatnames", "TClonesArray", &tca_plotvarcatnames, 32000, 0);
  int catvartemp=0;
  for(int i=0; i<Nvarcats; i++) {
    for(int j=0; j<ncats[i]; j++) {
      new ((*tca_plotvarcatnames)[catvartemp]) TObjString(); 
      ((TObjString *)tca_plotvarcatnames->At(catvartemp))->SetString(catnames[i][j]);
      catvartemp++;
    } 
  } 
  std::cout << "Ncatvar: " << Ncatvar << std::endl;
  std::cout << "catvartemp: " << catvartemp << std::endl;

  plotvartree->Fill();
  plotvartree->Write(0,TObject::kWriteDelete);

//int dummy = fscanf(file,"%d intL=%f rtree=%d %s %s\n",&nfiles, &intlumi, &makeOutputTree, outFilNam, histFilNam);
//int dummy = fscanf(file,"typ=%d ind=%d draw=%d Nam=%s Fil=%s tot=%d red=%d lum=%f xsec=%f kfac=%f scal=%f\n", &itype[i], &histoindfromfiles[i], &histoplotit[i], filesshortnam[i], files[i],&ntot[i],&nred[i], &lumi[i], &xsec[i], &kfactor[i], &scale[i]);
  inputfiletree = new TTree("inputfiles","globe inputfiles provenance information");
  inputfiletree->Branch("nfiles", &mp->nfiles, "nfiles/I");
  inputfiletree->Branch("nindfiles", &mp->nindfiles, "nindfiles/I");
  inputfiletree->Branch("intlumi", &mp->intlumi, "intlumi/F");
  inputfiletree->Branch("makeOutputTree", &mp->makeOutputTree, "makeOutputTree/I");
  tca_histfilename = new TClonesArray("TObjString",1);
  inputfiletree->Branch("histfilename", "TClonesArray", &tca_histfilename, 32000, 0);
  new ((*tca_histfilename)[0]) TObjString(); 
  ((TObjString *)tca_histfilename->At(0))->SetString(mp->histFilNam);
  inputfiletree->Branch("itype", &mp->itype, "itype[nfiles]/I");
  inputfiletree->Branch("histoind", &mp->histoindfromfiles, "histoindfromfiles[nfiles]/I");
  inputfiletree->Branch("infoind", &mp->infoind, "infoind[nindfiles]/I");
  inputfiletree->Branch("histoplotit", &mp->histoplotit, "histoplotit[nfiles]/I");
  inputfiletree->Branch("ntot", &mp->ntot, "ntot[nfiles]/I");
  inputfiletree->Branch("nred", &mp->nred, "nred[nfiles]/I");
  inputfiletree->Branch("lumi", &mp->lumi, "lumi[nfiles]/F");
  inputfiletree->Branch("xsec", &mp->xsec, "xsec[nfiles]/F");
  inputfiletree->Branch("kfactor", &mp->kfactor, "kfactor[nfiles]/F");
  inputfiletree->Branch("scale", &mp->scale, "scale[nfiles]/F");
  tca_inshortnames = new TClonesArray("TObjString",1);
  tca_infilenames = new TClonesArray("TObjString",mp->nfiles);
  inputfiletree->Branch("inshortnames", "TClonesArray", &tca_inshortnames, 32000, 0);
  inputfiletree->Branch("infilenames", "TClonesArray", &tca_infilenames, 32000, 0);
  for(int ifile=0;ifile!=mp->nfiles;++ifile) { 
    new ((*tca_inshortnames)[ifile]) TObjString(); 
    new ((*tca_infilenames)[ifile]) TObjString(); 
    ((TObjString *)tca_inshortnames->At(ifile))->SetString(mp->filesshortnam[ifile]);
    ((TObjString *)tca_infilenames->At(ifile))->SetString(mp->files[ifile]);
  }
  inputfiletree->Fill();
  inputfiletree->Write(0,TObject::kWriteDelete);




In mpUtil:

void mpUtil::PlotInteractive(TString tag, TString inputfile) { 

      loops->myPlotInteractiveSetup(this,inputfile,tag); 
      if(tag == "Hgg") {
        //loops->myPlotInteractiveHggSetup(this,inputfile);
      } else if(tag == "Hww") {
	//loops->myPlotInteractiveHwwSetup(this,inputfile); 
      } else if(tag == "Hzz") {
      } else if(tag == "Mwl") {
      } else if(tag == "Elizabeth") {
	//loops->myPlotInteractiveElizabethSetup(this,inputfile); 
      }
      loops->myPlotInteractive(this,inputfile); 
}

Script:




root -l
//.L matteo/tdrstyle.C
//setTDRStyle();
  gSystem->Load("libPhysics.so");
  gSystem->Load("libCore.so");
  //gSystem->Load("VertexAnalysis/lib/libh2gglobeVertexAnalysis.so");
gSystem->Load("libLoopAll.so");
LoopAll* m=new LoopAll();
m->myPlotInteractiveSetup("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root","Hgg");
m->myPlotInteractive("Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root");



 Marco/histograms_CMS-HGG_4686_nov17_tryExcl.root



mpUtil* m=new mpUtil();
m->PlotInteractive("Hgg","oct13_ucsdcuts_smnewnew2.root");


varname[0]: all_mass
varname[1]: pt
varname[2]: eta
varname[3]: mass
varname[4]: helicityAngle
varname[5]: decayAngle
varname[6]: pho_pt
varname[7]: pho1_pt
varname[8]: pho2_pt
varname[9]: pho_eta
varname[10]: pho1_eta
varname[11]: pho2_eta
varname[12]: pho_r9
varname[13]: pho_n
varname[14]: mass_pf
varname[15]: mass_ff
varname[16]: pho1_pt_presel
varname[17]: pho2_pt_presel
varname[18]: pho1_pt_sel
varname[19]: pho2_pt_sel
varname[20]: pho1_eta_presel
varname[21]: pho2_eta_presel
varname[22]: pho1_eta_sel
varname[23]: pho2_eta_sel
varname[24]: cut_Mgg_nminus1
varname[25]: cut_Mgg_sequential
varname[26]: cut_VBFLeadJPt_nminus1
varname[27]: cut_VBFLeadJPt_sequential
varname[28]: cut_VBFSubJPt_nminus1
varname[29]: cut_VBFSubJPt_sequential
varname[30]: cut_VBF_Mjj_nminus1
varname[31]: cut_VBF_Mjj_sequential
varname[32]: cut_VBF_dEta_nminus1
varname[33]: cut_VBF_dEta_sequential
varname[34]: cut_VBF_Zep_nminus1
varname[35]: cut_VBF_Zep_sequential
varname[36]: cut_VBF_dPhi_nminus1
varname[37]: cut_VBF_dPhi_sequential
varname[38]: cut_VBF_Mgg_nminus1
varname[39]: cut_VBF_Mgg_sequential
varname[40]: cut_VBF_Mgg4_nminus1
varname[41]: cut_VBF_Mgg4_sequential

