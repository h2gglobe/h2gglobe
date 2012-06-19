if (! -d $2) then
        mkdir $2
endif
cd $2
if ( "$3" == "" ) then
	set maxthreads = "4"
else
	set maxthreads = $3
endif
echo "Max Threads: $maxthreads"
foreach i (`seq 110 0.5 150`)
	set threads = `ps | grep combine | wc | awk '{print $1}'`
	while ($threads >= $maxthreads)
		echo "The Computer is Tired: Resting for 10 seconds"
		sleep 10
		set threads = `ps | grep combine | wc | awk '{print $1}'`
	end
	echo "Subitting Mass ${i}GeV"
	if ( -d ../$1) then
		combine ../$1/${i}GeVmodel.root -m $i -M Asymptotic --minimizerStrategy=1 --minosAlgo=stepping -s -1 --picky >&! higgsCombineTest.Asymptotic.mH${i}.log &
	else if ( -f ../$1) then
		combine ../$1 -m $i -M Asymptotic --minimizerStrategy=1 --minosAlgo=stepping -s -1 --picky >&! higgsCombineTest.Asymptotic.mH${i}.log &
	endif
	sleep 2
end
sleep 120
foreach i (`/bin/ls higgsCombineTest.Asymptotic.*.root`)
	set targetfile = `echo "$i" | sed 's|higgsCombineTest.Asymptotic.mH\(1[1-5][0-9][\.5]*\)\.[0-9-][0-9-]*.root|higgsCombineTest.Asymptotic.mH\1.root|'`
	mv $i $targetfile
end
cd ..
