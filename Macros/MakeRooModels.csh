#!/bin/tcsh
eval `scramv1 runtime -csh`
if (! -d $2) then
	mkdir $2
endif
foreach i ( `seq 110 0.5 150` )
	echo "Making ${i}GeVmodel.root"
	text2workspace.py -m ${i} $1 -b -o ${2}/${i}GeVmodel.root
end
