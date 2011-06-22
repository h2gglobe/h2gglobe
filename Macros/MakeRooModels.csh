#!/bin/tcsh
foreach i ( `seq 105 5 140` )
	text2workspace.py -m ${i}.000000 -D data_mass cms-hgg-datacard_parBKG_3sigma.txt -b -o ${i}GeVmodel.root
end
