source version.sh

./count_all_events.sh data_2012 > ${version}_summary.txt
###./count_all_events.sh mc_sig_summer12_s7 >> ${version}_summary.txt
./count_all_events.sh mc_sig_summer12_s7_03_06 >> ${version}_summary.txt
### ./count_all_events.sh mc_sig_summer12_s7_03_05 >> ${version}_summary.txt
### ./count_all_events.sh mc_sig_summer12_s7_03_04 >> ${version}_summary.txt
./count_all_events.sh mc_bkg_summer12_s7 >> ${version}_summary.txt
### ./count_all_events.sh mc_bkg_summer12_s7_b >> ${version}_summary.txt
