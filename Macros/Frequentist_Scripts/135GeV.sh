
#!/bin/bash
#############################################################
#
# Driver script for creating Hybrid or Frequentist grids
#
# author: Giovanni Petrucciani, UCSD                       
#         from a similar script by Luca Lista, INFN        
#
##############################################################

i="$1"
if [ "$i" = "" ]; then
  echo "Error: missing job index"
  exit 1;
fi
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
if [ "$n" = "" ]; then
  n="$2"
fi
if [ "$n" = "" ]; then
  echo "Error: missing number of experiments"
  exit 2;
fi

echo "## Starting at $(date)"
(( ($i + 0) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 0.0
(( ($i + 1) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 0.333333333333
(( ($i + 2) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 0.666666666667
(( ($i + 3) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 1.0
(( ($i + 4) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 1.33333333333
(( ($i + 5) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 1.66666666667
(( ($i + 6) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 2.0
(( ($i + 7) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 2.33333333333
(( ($i + 8) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 2.66666666667
(( ($i + 9) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 3.0
(( ($i + 10) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 3.33333333333
(( ($i + 11) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 3.66666666667
(( ($i + 12) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 4.0
(( ($i + 13) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 4.33333333333
(( ($i + 14) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 4.66666666667
(( ($i + 15) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 5.0
(( ($i + 16) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 5.33333333333
(( ($i + 17) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 5.66666666667
(( ($i + 18) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 6.0
(( ($i + 19) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 6.33333333333
(( ($i + 20) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 6.66666666667
(( ($i + 21) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 7.0
(( ($i + 22) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 7.33333333333
(( ($i + 23) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 7.66666666667
(( ($i + 24) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 8.0
(( ($i + 25) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 8.33333333333
(( ($i + 26) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 8.66666666667
(( ($i + 27) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 9.0
(( ($i + 28) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 9.33333333333
(( ($i + 29) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 9.66666666667
(( ($i + 30) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 10.0
(( ($i + 31) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 10.3333333333
(( ($i + 32) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 10.6666666667
(( ($i + 33) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 11.0
(( ($i + 34) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 11.3333333333
(( ($i + 35) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 11.6666666667
(( ($i + 36) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 12.0
(( ($i + 37) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 12.3333333333
(( ($i + 38) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 12.6666666667
(( ($i + 39) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 13.0
(( ($i + 40) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 13.3333333333
(( ($i + 41) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 13.6666666667
(( ($i + 42) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 14.0
(( ($i + 43) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 14.3333333333
(( ($i + 44) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 14.6666666667
(( ($i + 45) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 15.0
(( ($i + 46) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 15.3333333333
(( ($i + 47) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 15.6666666667
(( ($i + 48) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 16.0
(( ($i + 49) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 16.3333333333
(( ($i + 50) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 16.6666666667
(( ($i + 51) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 17.0
(( ($i + 52) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 17.3333333333
(( ($i + 53) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 17.6666666667
(( ($i + 54) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 18.0
(( ($i + 55) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 18.3333333333
(( ($i + 56) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 18.6666666667
(( ($i + 57) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 19.0
(( ($i + 58) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 19.3333333333
(( ($i + 59) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 19.6666666667
(( ($i + 60) % 1 == 0 )) &&  ./combine 135GeVmodel.root -M HybridNew -D data_mass -m 135 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 135GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 20.0
if [ -d /tmp/drberry/ ]; then
	rm /tmp/drberry/*.txt
fi
/bin/ls
hadd 135GeVFrequentist.root higgsCombine*.root
echo "## Done at $(date)"
