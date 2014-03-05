#!/bin/bash
if [ -e outputToy ]; then
	rm -rf outputToy
fi
mkdir outputToy
echo "max events from CRAB: $MaxEvents"
n="$MaxEvents"
if [ "$n" = "" ]; then
	n="$1"
fi
if [ "$n" = "" ]; then
	echo "Error: missing number of experiments"
	exit 2;
fi
if ( ./combine hgg_card_spin_comb_bernsteins.root -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T $n -i 1 --clsAcc 0 --fullBToys --redefineSignalPOIs x --setPhysicsModelParameters fqq=0.0 --freezeNuisances fqq -s 0 -n job0fqq0.00 && ./combine hgg_card_spin_comb_bernsteins.root -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T $n -i 1 --clsAcc 0 --fullBToys --redefineSignalPOIs x --setPhysicsModelParameters fqq=0.25 --freezeNuisances fqq -s 0 -n job0fqq0.25 && ./combine hgg_card_spin_comb_bernsteins.root -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T $n -i 1 --clsAcc 0 --fullBToys --redefineSignalPOIs x --setPhysicsModelParameters fqq=0.50 --freezeNuisances fqq -s 0 -n job0fqq0.50 && ./combine hgg_card_spin_comb_bernsteins.root -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T $n -i 1 --clsAcc 0 --fullBToys --redefineSignalPOIs x --setPhysicsModelParameters fqq=0.75 --freezeNuisances fqq -s 0 -n job0fqq0.75 && ./combine hgg_card_spin_comb_bernsteins.root -M HybridNew --testStat=TEV --generateExt=1 --generateNuis=0 --fitNuis=1 --singlePoint 1 --saveHybridResult -T $n -i 1 --clsAcc 0 --fullBToys --redefineSignalPOIs x --setPhysicsModelParameters fqq=1.00 --freezeNuisances fqq -s 0 -n job0fqq1.00 >& outputToy/subqqbar.sh.log )
	then touch outputToy/subqqbar.sh.done
	else touch outputToy/subqqbar.sh.fail
fi
mv higgsCombine*.root outputToy/ 
echo "pack the results"
tar cvfz outputToy.tgz outputToy
