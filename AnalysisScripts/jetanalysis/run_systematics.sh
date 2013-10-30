#!/bin/bash

export DISPLAY=""

dir="yr3_systematics_v1"
mkdir $dir

tunes="TuneD6T TuneD6TUEOFF TuneP0 TuneP0UEOFF TuneProPT0 TuneProPT0UEOFF TuneProQ20 TuneProQ20UEOFF TuneUEOFF TuneZ2Star"
## tunes="TuneP0 TuneProPT0 TuneUEOFF TuneZ2Star TuneD6T TuneProQ20"

for tune in $tunes; do
    rm jetanalysis/tmp_datafiles_mvavbftag_tunes.dat{,.pevents}
    ./mk_fitter.py -i jetanalysis/datafiles_mvavbftag_tunes.dat -n 5 -l ${tune} -o ${dir}/mavtag_${tune}/sub &&  \
     	./submit_fitter.py -q 8nh -d ${dir}/mavtag_${tune}
    ### rm jetanalysis/tmp_datafiles_cutbasedvbftag_tunes.dat{,.pevents}
    ### ./mk_fitter.py -i jetanalysis/datafiles_cutbasedvbftag_tunes.dat -n 5 -l ${tune} -o ${dir}/cutbtag_${tune}/sub && \
    ###	./submit_fitter.py -q 8nh -d ${dir}/cutbtag_${tune}
done

## dir="yr3_systematics_v1"
## mkdir $dir

jecUncs=""
## jecUncs="jerCentral"
## jecUncs="nominal jecUp jecDown jerCentral jerUp jerDown"

for syst in $jecUncs; do
    ./mk_fitter.py -i jetanalysis/datafiles_mvavbftag_jec.dat -n 5 -l ${syst} -o ${dir}/mvatag_${syst}/sub && \
	./submit_fitter.py -q 8nh -d ${dir}/mvatag_${syst}
    ./mk_fitter.py -i jetanalysis/datafiles_cutbasedvbftag_jec.dat -n 5 -l ${syst} -o ${dir}/cutbtag_${syst}/sub && \
	./submit_fitter.py -q 8nh -d ${dir}/cutbtag_${syst}
done

btagEff=" noBtagSF"
#btagEff=" noBtagSF nominalBtagSF shiftBtagEffUp_bc shiftBtagEffDown_bc shiftBtagEffUp_l shiftBtagEffDown_l"


for systb in $btagEff; do
    ./mk_fitter.py -i jetanalysis/datafiles_mva_btag.dat -n 5 -l ${systb} -o ${dir}/mvabtag_${systb}/sub && \                                               
        ./submit_fitter.py -q 8nh -d ${dir}/btag_${systb}
    ./mk_fitter.py -i jetanalysis/datafiles_cutbased_btag.dat -n 15 -l ${systb} -o ${dir}/cutbtag_${systb}/sub && \
	./submit_fitter.py -q 8nh -d ${dir}/btag_${systb}
done