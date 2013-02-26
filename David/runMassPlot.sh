#!/bin/csh

setenv OUTDIR "/afs/cern.ch/user/f/futyand/www/hgg/moriond_preapproval/dataMC/cic"
#mkdir $OUTDIR/mass $OUTDIR/pho1_pt $OUTDIR/pho2_pt $OUTDIR/pho1_r9 $OUTDIR/pho2_r9 $OUTDIR/pho1_eta $OUTDIR/pho2_eta $OUTDIR/pt $OUTDIR/eta
mkdir $OUTDIR/mass $OUTDIR/pho1_pt $OUTDIR/pho2_pt $OUTDIR/pho_r9 $OUTDIR/pho1_eta $OUTDIR/pho2_eta $OUTDIR/pt $OUTDIR/eta

root -b -l -q 'massPlot.C("0","all_mass","cic")'
root -b -l -q 'massPlot.C("1","all_mass","cic")'
root -b -l -q 'massPlot.C("2","all_mass","cic")'
root -b -l -q 'massPlot.C("3","all_mass","cic")'
root -b -l -q 'massPlot.C("4","all_mass","cic")'
root -b -l -q 'massPlot.C("5","all_mass","cic")'
root -b -l -q 'massPlot.C("6","all_mass","cic")'

root -b -l -q 'massPlot.C("0","pho1_pt","cic")'
root -b -l -q 'massPlot.C("1","pho1_pt","cic")'
root -b -l -q 'massPlot.C("2","pho1_pt","cic")'
root -b -l -q 'massPlot.C("3","pho1_pt","cic")'
root -b -l -q 'massPlot.C("4","pho1_pt","cic")'

#root -b -l -q 'massPlot.C("0","pho1_r9","cic")'
#root -b -l -q 'massPlot.C("1","pho1_r9","cic")'
#root -b -l -q 'massPlot.C("2","pho1_r9","cic")'
#root -b -l -q 'massPlot.C("3","pho1_r9","cic")'
#root -b -l -q 'massPlot.C("4","pho1_r9","cic")'

root -b -l -q 'massPlot.C("0","pho1_eta","cic")'
root -b -l -q 'massPlot.C("1","pho1_eta","cic")'
root -b -l -q 'massPlot.C("2","pho1_eta","cic")'
root -b -l -q 'massPlot.C("3","pho1_eta","cic")'
root -b -l -q 'massPlot.C("4","pho1_eta","cic")'

root -b -l -q 'massPlot.C("0","pho2_pt","cic")'
root -b -l -q 'massPlot.C("1","pho2_pt","cic")'
root -b -l -q 'massPlot.C("2","pho2_pt","cic")'
root -b -l -q 'massPlot.C("3","pho2_pt","cic")'
root -b -l -q 'massPlot.C("4","pho2_pt","cic")'

#root -b -l -q 'massPlot.C("0","pho2_r9","cic")'
#root -b -l -q 'massPlot.C("1","pho2_r9","cic")'
#root -b -l -q 'massPlot.C("2","pho2_r9","cic")'
#root -b -l -q 'massPlot.C("3","pho2_r9","cic")'
#root -b -l -q 'massPlot.C("4","pho2_r9","cic")'

root -b -l -q 'massPlot.C("0","pho2_eta","cic")'
root -b -l -q 'massPlot.C("1","pho2_eta","cic")'
root -b -l -q 'massPlot.C("2","pho2_eta","cic")'
root -b -l -q 'massPlot.C("3","pho2_eta","cic")'
root -b -l -q 'massPlot.C("4","pho2_eta","cic")'

root -b -l -q 'massPlot.C("0","pt","cic")'
root -b -l -q 'massPlot.C("1","pt","cic")'
root -b -l -q 'massPlot.C("2","pt","cic")'
root -b -l -q 'massPlot.C("3","pt","cic")'
root -b -l -q 'massPlot.C("4","pt","cic")'

root -b -l -q 'massPlot.C("0","eta","cic")'
root -b -l -q 'massPlot.C("1","eta","cic")'
root -b -l -q 'massPlot.C("2","eta","cic")'
root -b -l -q 'massPlot.C("3","eta","cic")'
root -b -l -q 'massPlot.C("4","eta","cic")'

root -b -l -q 'massPlot.C("0","pho_r9","cic")'
root -b -l -q 'massPlot.C("1","pho_r9","cic")'
root -b -l -q 'massPlot.C("2","pho_r9","cic")'
root -b -l -q 'massPlot.C("3","pho_r9","cic")'
root -b -l -q 'massPlot.C("4","pho_r9","cic")'
