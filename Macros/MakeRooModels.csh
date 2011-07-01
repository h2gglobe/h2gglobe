#!/bin/tcsh
eval `scramv1 runtime -csh`
if (! -d RooModel) then
	mkdir RooModel
endif
set DATACARD = $1
if ($?DATACARD) then
	set DATACARD = cms-hgg-datacard_parBKG.txt
endif
foreach i ( `seq 105 0.5 140` )
	echo "Making ${i}GeVmodel.root"
	text2workspace.py -m ${i} -D data_mass $DATACARD -b -o RooModel/${i}GeVmodel.root
end
