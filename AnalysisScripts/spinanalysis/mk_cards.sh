#!/bin/bash
dfile=$1
sfile=$2
./mk_spin_card.py -n datacard_spinanalysis.txt --nCosThetaCats 5 --nKinCats 4 -B ${dfile} -S ${sfile}
./mk_spin_card.py -n datacard_spinanalysis_qqbar.txt --nCosThetaCats 5 --nKinCats 4 -B ${dfile} -S ${sfile} -Q 
./mk_spin_card.py -n datacard_spinanalysis_justSM.txt --nCosThetaCats 5 --nKinCats 4 -B ${dfile} -S ${sfile} -s
./mk_spin_card.py -n datacard_spinanalysis_justGravGG.txt --nCosThetaCats 5 --nKinCats 4 -B ${dfile} -S ${sfile} -g
./mk_spin_card.py -n datacard_spinanalysis_justGravQQ.txt --nCosThetaCats 5 --nKinCats 4 -B ${dfile} -S ${sfile} -q

#./mk_spin_card.py -n datacard_spinanalysis_noCTcats.txt --nCosThetaCats 1 --nKinCats 4 -B ${dfile} -S ${sfile}
#./mk_spin_card.py -n datacard_spinanalysis_qqbar_noCTcats.txt --nCosThetaCats 1 --nKinCats 4 -B ${dfile} -S ${sfile} -Q 
#./mk_spin_card.py -n datacard_spinanalysis_justSM_noCTcats.txt --nCosThetaCats 1 --nKinCats 4 -B ${dfile} -S ${sfile} -s
#./mk_spin_card.py -n datacard_spinanalysis_justGravGG_noCTcats.txt --nCosThetaCats 1 --nKinCats 4 -B ${dfile} -S ${sfile} -g
#./mk_spin_card.py -n datacard_spinanalysis_justGravQQ_noCTcats.txt --nCosThetaCats 1 --nKinCats 4 -B ${dfile} -S ${sfile} -q

