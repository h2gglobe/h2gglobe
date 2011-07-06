#!/bin/tcsh
eval `scramv1 runtime -csh`
foreach i ( `seq 110 0.5 140` )
	echo "Making ${i}GeVmodel.root"
	text2workspace.py -m ${i} -D data_mass cms-hgg-datacard_CMS-HGG_204pb.root_parBKG.txt -b -o RooModel204/${i}GeVmodel.root
	text2workspace.py -m ${i} -D data_mass cms-hgg-datacard_CMS-HGG_951pb.root_parBKG.txt -b -o RooModel951/${i}GeVmodel.root
end
