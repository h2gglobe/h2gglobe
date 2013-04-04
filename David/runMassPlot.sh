#!/bin/csh

# replace massfit with cic to make plots for cut based analysis

setenv OUTDIR "/afs/cern.ch/user/f/futyand/www/hgg/moriond_preapproval/dataMC/massfit"
mkdir $OUTDIR/mass $OUTDIR/pho1_pt $OUTDIR/pho2_pt $OUTDIR/pho_r9 $OUTDIR/pho1_eta $OUTDIR/pho2_eta $OUTDIR/pt $OUTDIR/eta

root -b -l -q 'massPlot.C("0","all_mass","massfit")'
root -b -l -q 'massPlot.C("1","all_mass","massfit")'
root -b -l -q 'massPlot.C("2","all_mass","massfit")'
root -b -l -q 'massPlot.C("3","all_mass","massfit")'
root -b -l -q 'massPlot.C("4","all_mass","massfit")'
root -b -l -q 'massPlot.C("5","all_mass","massfit")'
root -b -l -q 'massPlot.C("6","all_mass","massfit")'

root -b -l -q 'massPlot.C("0","pho1_pt","massfit")'
root -b -l -q 'massPlot.C("1","pho1_pt","massfit")'
root -b -l -q 'massPlot.C("2","pho1_pt","massfit")'
root -b -l -q 'massPlot.C("3","pho1_pt","massfit")'
root -b -l -q 'massPlot.C("4","pho1_pt","massfit")'

root -b -l -q 'massPlot.C("0","pho1_eta","massfit")'
root -b -l -q 'massPlot.C("1","pho1_eta","massfit")'
root -b -l -q 'massPlot.C("2","pho1_eta","massfit")'
root -b -l -q 'massPlot.C("3","pho1_eta","massfit")'
root -b -l -q 'massPlot.C("4","pho1_eta","massfit")'

root -b -l -q 'massPlot.C("0","pho2_pt","massfit")'
root -b -l -q 'massPlot.C("1","pho2_pt","massfit")'
root -b -l -q 'massPlot.C("2","pho2_pt","massfit")'
root -b -l -q 'massPlot.C("3","pho2_pt","massfit")'
root -b -l -q 'massPlot.C("4","pho2_pt","massfit")'

root -b -l -q 'massPlot.C("0","pho2_eta","massfit")'
root -b -l -q 'massPlot.C("1","pho2_eta","massfit")'
root -b -l -q 'massPlot.C("2","pho2_eta","massfit")'
root -b -l -q 'massPlot.C("3","pho2_eta","massfit")'
root -b -l -q 'massPlot.C("4","pho2_eta","massfit")'

root -b -l -q 'massPlot.C("0","pt","massfit")'
root -b -l -q 'massPlot.C("1","pt","massfit")'
root -b -l -q 'massPlot.C("2","pt","massfit")'
root -b -l -q 'massPlot.C("3","pt","massfit")'
root -b -l -q 'massPlot.C("4","pt","massfit")'

root -b -l -q 'massPlot.C("0","eta","massfit")'
root -b -l -q 'massPlot.C("1","eta","massfit")'
root -b -l -q 'massPlot.C("2","eta","massfit")'
root -b -l -q 'massPlot.C("3","eta","massfit")'
root -b -l -q 'massPlot.C("4","eta","massfit")'

root -b -l -q 'massPlot.C("0","pho_r9","massfit")'
root -b -l -q 'massPlot.C("1","pho_r9","massfit")'
root -b -l -q 'massPlot.C("2","pho_r9","massfit")'
root -b -l -q 'massPlot.C("3","pho_r9","massfit")'
root -b -l -q 'massPlot.C("4","pho_r9","massfit")'
