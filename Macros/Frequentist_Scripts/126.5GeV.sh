
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

/bin/ls
echo "## Starting at $(date)"
(( ($i + 0) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 0.0
(( ($i + 1) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 1.01694915254
(( ($i + 2) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 2.03389830508
(( ($i + 3) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 3.05084745763
(( ($i + 4) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 4.06779661017
(( ($i + 5) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 5.08474576271
(( ($i + 6) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 6.10169491525
(( ($i + 7) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 7.1186440678
(( ($i + 8) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 8.13559322034
(( ($i + 9) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 9.15254237288
(( ($i + 10) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 10.1694915254
(( ($i + 11) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 11.186440678
(( ($i + 12) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 12.2033898305
(( ($i + 13) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 13.2203389831
(( ($i + 14) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 14.2372881356
(( ($i + 15) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 15.2542372881
(( ($i + 16) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 16.2711864407
(( ($i + 17) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 17.2881355932
(( ($i + 18) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 18.3050847458
(( ($i + 19) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 19.3220338983
(( ($i + 20) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 20.3389830508
(( ($i + 21) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 21.3559322034
(( ($i + 22) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 22.3728813559
(( ($i + 23) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 23.3898305085
(( ($i + 24) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 24.406779661
(( ($i + 25) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 25.4237288136
(( ($i + 26) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 26.4406779661
(( ($i + 27) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 27.4576271186
(( ($i + 28) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 28.4745762712
(( ($i + 29) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 29.4915254237
(( ($i + 30) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 30.5084745763
(( ($i + 31) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 31.5254237288
(( ($i + 32) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 32.5423728814
(( ($i + 33) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 33.5593220339
(( ($i + 34) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 34.5762711864
(( ($i + 35) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 35.593220339
(( ($i + 36) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 36.6101694915
(( ($i + 37) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 37.6271186441
(( ($i + 38) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 38.6440677966
(( ($i + 39) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 39.6610169492
(( ($i + 40) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 40.6779661017
(( ($i + 41) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 41.6949152542
(( ($i + 42) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 42.7118644068
(( ($i + 43) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 43.7288135593
(( ($i + 44) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 44.7457627119
(( ($i + 45) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 45.7627118644
(( ($i + 46) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 46.7796610169
(( ($i + 47) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 47.7966101695
(( ($i + 48) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 48.813559322
(( ($i + 49) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 49.8305084746
(( ($i + 50) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 50.8474576271
(( ($i + 51) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 51.8644067797
(( ($i + 52) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 52.8813559322
(( ($i + 53) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 53.8983050847
(( ($i + 54) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 54.9152542373
(( ($i + 55) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 55.9322033898
(( ($i + 56) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 56.9491525424
(( ($i + 57) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 57.9661016949
(( ($i + 58) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 58.9830508475
(( ($i + 59) % 1 == 0 )) &&  ./combine 126.5GeVmodel.root -M HybridNew -D data_mass -m 126.5 --freq --fork 1 -T 50 --clsAcc 0 -v 0 -n 126.5GeV --saveHybridResult --saveToys -s -1 -i $n --singlePoint 60.0

if [ -d /tmp/drberry/ ]; then
        rm /tmp/drberry/*.txt
fi
/bin/ls
hadd 126.5GeV.root higgsCombine*.root
echo "## Done at $(date)"
