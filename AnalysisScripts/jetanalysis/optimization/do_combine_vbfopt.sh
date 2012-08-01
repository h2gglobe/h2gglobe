#!/bin/bash

export DISPLAY=""

dir=$(echo $1 | sed 's%/$%%')

cd ../Macros/

export PATH=$PWD/ResultScripts:$PATH

[[ ! -e $dir ]] && ln -s ../AnalysisScripts/$dir

datacards="datacard_full_x6_fit.txt datacard_vbf.txt datacard_vbf_plus_boosted.txt datacard_vbf_plus_boosted_x6.txt datacard_vbf_x6.txt"

for dc in $datacards; do
    [[ ! -f $dir/$dc ]] && cp -p ../AnalysisScripts/jetanalysis/optimization/$dc $dir
done

[[ ! -f $dir/CMS-HGG_test_interpolated.root ]] && ./massInterpolator.py --outputMasses=125 --inputMasses=125 --doSmoothing -i $dir/CMS-HGG_test.root
[[ ! -f $dir/CMS-HGG_test_x6_interpolated.root ]] && ./massInterpolator.py --outputMasses=125 --inputMasses=125 --doSmoothing -i $dir/CMS-HGG_test.root -o $dir/CMS-HGG_test_x6 -k 6.

cd $dir

parallel -j 5 --eta 'combine --verbose=2 -n expected_pval_vbf_x6 -M ProfileLikelihood --signif -t -1 --expectSignal=1 --pvalue -m {} -S 0 datacard_vbf_x6.txt 2>&1 | tee combine_{}_pvalue_vbf_exp_x6.log' ::: 125

parallel -j 5 --eta 'combine --verbose=2 -n expected_pval_vbf -M ProfileLikelihood --signif -t -1 --expectSignal=1 --pvalue -m {} -S 0 datacard_vbf.txt 2>&1 | tee combine_{}_pvalue_vbf_exp.log' ::: 125

parallel -j 5 --eta 'combine --verbose=2 -n expected_pval -M ProfileLikelihood --signif -t -1 --expectSignal=1 --pvalue -m {} -S 0 datacard_vbf_plus_boosted.txt 2>&1 | tee combine_{}_pvalue_exp.log' ::: 125

parallel -j 5 --eta 'combine --verbose=2 -n expected_pval_x6 -M ProfileLikelihood --signif -t -1 --expectSignal=1 --pvalue -m {} -S 0 datacard_vbf_plus_boosted_x6.txt 2>&1 | tee combine_{}_pvalue_exp_x6.log' ::: 125
  
cd ..

multiDimFit.py -w $dir -d datacard_vbf_plus_boosted.txt -m 125. -M rV --npoints 100 --expected 1 -f
mkdir $dir/rV125 
mv $dir/*rV125*.* $dir/rV125

multiDimFit.py -w $dir -d datacard_vbf_plus_boosted_x6.txt -m 125. -M rV --npoints 100 --expected 1 -f 
mkdir $dir/rV125_x6
mv $dir/*rV125*.* $dir/rV125_x6

multiDimFit.py -w $dir -d datacard_full_x6_fit.txt -m 125. -M cVcF --npoints 1000 --expected 1 -f 
mkdir $dir/cVcF125_x6
mv $dir/*cVcF125*.* $dir/cVcF125_x6

multiDimFit.py -w $dir -d datacard_full_x6_fit.txt -m 125. -M rVrF --npoints 1000 --expected 1 -f 
mkdir $dir/rVrF125_x6
mv $dir/*rVrF125*.* $dir/rVrF125_x6

