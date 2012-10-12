#!/bin/bash

export DISPLAY=""

tunes=""
## tunes="TuneP0 TuneProPT0 TuneUEOFF TuneZ2Star TuneD6T TuneProQ20"

for tune in $tunes; do
    rm jetanalysis/tmp_datafiles_mvavbftag_tunes.dat{,.pevents}
    ./mk_fitter.py -i jetanalysis/datafiles_mvavbftag_tunes.dat -n 5 -l ${tune} -o tunes_v2/mavtag_${tune}/sub &&  \
     	./submit_fitter.py -q 8nh -d tunes_v2/mavtag_${tune}
    rm jetanalysis/tmp_datafiles_cutbasedvbftag_tunes.dat{,.pevents}
    ./mk_fitter.py -i jetanalysis/datafiles_cutbasedvbftag_tunes.dat -n 5 -l ${tune} -o tunes_v2/cutbtag_${tune}/sub && \
	./submit_fitter.py -q 8nh -d tunes_v2/cutbtag_${tune}
done


## jecUncs="jerCentral"
jecUncs="nominal jecUp jecDown jerCentral jerUp jerDown"

for syst in $jecUncs; do
    ./mk_fitter.py -i jetanalysis/datafiles_mvavbftag_jec.dat -n 5 -l ${syst} -o jec_v3/mvatag_${syst}/sub && \
	./submit_fitter.py -q 8nh -d jec_v3/mvatag_${syst}
    ./mk_fitter.py -i jetanalysis/datafiles_cutbasedvbftag_jec.dat -n 5 -l ${syst} -o jec_v3/cutbtag_${syst}/sub && \
	./submit_fitter.py -q 8nh -d jec_v3/cutbtag_${syst}
done
