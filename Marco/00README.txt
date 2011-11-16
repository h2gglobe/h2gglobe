UCSD:
-----

cd ~/Aug2011/

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
cvs co -r branches_four_42X -d HiggsAnalysis/HiggsTo2photons UserCode/HiggsAnalysis/HiggsTo2photons 

#cvs co -r nw_11_11_11_regr_otf -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe

cvs co -d HiggsAnalysis/HiggsTo2photons/h2gglobe UserCode/HiggsAnalysis/HiggsTo2photons/h2gglobe

scramv1 b

cd HiggsAnalysis/HiggsTo2photons/h2gglobe
make clean; make