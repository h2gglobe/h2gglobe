#!/bin/bash

./mk_spin_card.py -n datacard_spinanalysis.txt
./mk_spin_card.py -n datacard_spinanalysis_qqbar.txt -q
./mk_spin_card.py -n datacard_spinanalysis_justSM.txt -s
./mk_spin_card.py -n datacard_spinanalysis_justGrav.txt -g

