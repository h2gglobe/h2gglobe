#!/bin/tcsh
eval `scramv1 runtime -csh`
foreach i ( `seq 110 0.5 140` )
	echo "Making ${i}GeVmodel.root"
	text2workspace.py -m ${i} -D data_mass cms-hgg-datacard_CMS-HGG.root_parBKG.txt -b -o RooModel/${i}GeVmodel.root
end
