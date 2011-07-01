#!/bin/tcsh
foreach i ( `/bin/ls -d * | grep $1` )
	cd ${i}/res/
	set mass = `echo ${i} | sed 's|\([0-9.]*\).*|\1|'`
	hadd -f ${mass}FrequentistGrid.root *.root
	cd ../../
end
hadd FrequentistGrid.root *GeVFrequentistGrid/res/*FrequentistGrid.root
