#!/bin/bash
#1 Datacard
#2 DOOBSERVED
#3 DOEXPECTED
#4 DOCATEGORIES
#5 THREADS

DEBUG=0
CUSTUMCAT=""

SCRAM_ARCH=slc5_amd64_gcc462
export SCRAM_ARCH
eval `scramv1 runtime -sh`

if [ ! -e "$1" ]; then
	echo "Error, $1 not found! Please give valid datacrd."
	exit 1
fi

if [ "$5" == "" ]; then
	MAXTHREADS=`cat /proc/cpuinfo | grep processor | awk '{print $3}' | tail -1 | xargs -i echo '{}+1' | bc`
else
	MAXTHREADS=$5
fi
echo "Max Threads: $MAXTHREADS"

if [ "$2" == "" ]; then DOOBSERVED=1; else DOOBSERVED=$2; fi
if [ "$3" == "" ]; then DOEXPECTED=1; else DOEXPECTED=$3; fi
if [ "$4" == "" ]; then DOCATEGORIES=1; else DOCATEGORIES=$4; fi

MASSES=`seq 110 0.5 150`
THREADS=0
DIRNAME=`echo $1 | sed 's|.txt||'`
JOBNAME=$DIRNAME

if [ ! -d $DIRNAME -a $DOOBSERVED -eq 1 ]; then mkdir $DIRNAME; fi
if [ ! -d ${DIRNAME}_Expected -a $DOEXPECTED -eq 1 ]; then mkdir ${DIRNAME}_Expected; fi

echo "Calculating Observed and Expected PValues"
if [ $DOOBSERVED -eq 1 -o $DOEXPECTED -eq 1 ]; then
	for MASS in $MASSES; do
		if [ $DEBUG -gt 1 ]; then echo "THREADS: $THREADS and MAXTHREADS: $MAXTHREADS"; fi
		while [ ! $THREADS -lt $MAXTHREADS ]; do
			sleep 2
			THREADS=`ps | grep combine | wc -l`
		done
		if [ $DOOBSERVED -eq 1 ]; then combine $1 -m $MASS -M ProfileLikelihood -t 0 -s -1 -n PValue --signif --pvalue >& $DIRNAME/higgsCombinePValue.ProfileLikelihood.mH${MASS}.log & fi
		if [ $DOEXPECTED -eq 1 ]; then combine $1 -m $MASS -M ProfileLikelihood -t -1 -s -1 -n PValueExpected --signif --pvalue --expectSignal=1 >& ${DIRNAME}_Expected/higgsCombinePValueExpected.ProfileLikelihood.mH${MASS}.log & fi
		sleep 0.25
		THREADS=`ps | grep combine | wc -l`
	done
	while [ ! $THREADS -eq 0 ]; do
		sleep 2
		THREADS=`ps | grep combine | wc -l`
	done
	if [ $DOOBSERVED -eq 1 ]; then hadd higgsCombinePValue.${JOBNAME}.ProfileLikelihood.root higgsCombinePValue.ProfileLikelihood.mH[0-9][0-9]*.[0-9-][0-9]*.root >& $DIRNAME/higgsCombinePValue.ProfileLikelihood.log; fi
	if [ $DOEXPECTED -eq 1 ]; then hadd higgsCombinePValueExpected.${JOBNAME}.ProfileLikelihood.root higgsCombinePValueExpected.ProfileLikelihood.mH[0-9][0-9]*.[0-9-][0-9]*.root >& ${DIRNAME}_Expected/higgsCombinePValueExpected.ProfileLikelihood.log; fi
	if [ $DOOBSERVED -eq 1 ]; then mv higgsCombinePValue.ProfileLikelihood.mH[0-9][0-9]*.[0-9-][0-9]*.root $DIRNAME; fi
	if [ $DOEXPECTED -eq 1 ]; then mv higgsCombinePValueExpected.ProfileLikelihood.mH[0-9][0-9]*.[0-9-][0-9]*.root ${DIRNAME}_Expected; fi
fi

if [ $DOCATEGORIES -eq 1 ]; then
	CATAGORIES=`grep bin $1 | grep cat | grep -v combine | head -n 1 | sed 's|bin[ \t][ \t]*||'`
	if [ $DEBUG -gt 1 ]; then echo $CATAGORIES; fi
	for CAT in $CATAGORIES
	  do
	  echo "Making DataCard for $CAT"
	  VETO=`echo $CATAGORIES | sed -e "s|$CAT||" | sed 's|^[ ]||' | sed 's|[ ]$||' | sed 's/[ ][ ]*/|ch1_/g' | sed 's|^|ch1_|'`
	  OUTFILE=`echo $1 | sed "s|.txt|_$CAT.txt|"`
	  if [ $DEBUG -gt 1 ]; then echo $VETO; fi
	  combineCards.py --xc="$VETO" $1 >& $OUTFILE
	  THREADS=0
	  DIRNAME=`echo $OUTFILE | sed 's|.txt||'`
	  if [ ! -d $DIRNAME ]; then
		  mkdir $DIRNAME
	  fi
	  echo "Calculating PValues for $CAT"
	  for MASS in $MASSES; do
		  if [ $DEBUG -gt 1 ]; then echo "THREADS: $THREADS and MAXTHREADS: $MAXTHREADS"; fi
		  while [ ! $THREADS -lt $MAXTHREADS ]; do
			  sleep 2
			  THREADS=`ps | grep combine | wc -l`
		  done
		  combine $OUTFILE -m $MASS -M ProfileLikelihood -t 0 -s -1 -n PValue$CAT --signif --pvalue >& $DIRNAME/higgsCombinePValuePerCat.ProfileLikelihood.$CAT.mH${MASS}.log &
		  sleep 0.25
		  THREADS=`ps | grep combine | wc -l`
	  done
	  while [ ! $THREADS -eq 0 ]; do
		  sleep 2
		  THREADS=`ps | grep combine | wc -l`
	  done
	  hadd higgsCombine.${JOBNAME}.ProfileLikelihood.$CAT.root higgsCombinePValue$CAT.ProfileLikelihood.mH[0-9][0-9]*.[0-9-][0-9]*.root >& $DIRNAME/higgsCombine.ProfileLikelihood.$CAT.log
	  mv higgsCombinePValue$CAT.ProfileLikelihood.mH[0-9][0-9]*.[0-9-][0-9]*.root $DIRNAME
	done
fi

if [ "$CUSTUMCAT" != "" ]; then
	echo "Making DataCard for Custum Catagory $CUSTUMCAT"
	for CAT in $CATAGORIES
	  do
	  echo "Making DataCard for $CAT"
	  OUTFILE=`echo $1 | sed "s|.txt|_$CAT.txt|"`
	  if [ $DEBUG -gt 1 ]; then echo $CUSTUMCAT; fi
	  combineCards.py --xc="$CUSTUMCAT" $1 >& $OUTFILE
	  THREADS=0
	  DIRNAME=`echo $OUTFILE | sed 's|.txt||'`
	  if [ ! -d $DIRNAME ]; then
		  mkdir $DIRNAME
	  fi
	  echo "Calculating PValues for $CAT"
	  for MASS in $MASSES; do
		  if [ $DEBUG -gt 1 ]; then echo "THREADS: $THREADS and MAXTHREADS: $MAXTHREADS"; fi
		  while [ ! $THREADS -lt $MAXTHREADS ]; do
			  sleep 2
			  THREADS=`ps | grep combine | wc -l`
		  done
		  combine $OUTFILE -m $MASS -M ProfileLikelihood -t 0 -s -1 -n PValue$CAT --signif --pvalue >& $DIRNAME/higgsCombinePValuePerCat.ProfileLikelihood.$CAT.mH${MASS}.log &
		  sleep 0.25
		  THREADS=`ps | grep combine | wc -l`
	  done
	  while [ ! $THREADS -eq 0 ]; do
		  sleep 2
		  THREADS=`ps | grep combine | wc -l`
	  done
	  hadd higgsCombine.${JOBNAME}.ProfileLikelihood.$CAT.root higgsCombinePValue$CAT.ProfileLikelihood.mH[0-9][0-9]*.[0-9-][0-9]*.root >& $DIRNAME/higgsCombine.ProfileLikelihood.$CAT.log
	  mv higgsCombinePValue$CAT.ProfileLikelihood.mH[0-9][0-9]*.[0-9-][0-9]*.root $DIRNAME
	done
fi

echo -e '\nDone'
