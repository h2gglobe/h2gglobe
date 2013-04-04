#!/bin/csh

root -b -l -q 'diphomva_inputs.C(0,0)'
root -b -l -q diphomva_inputs.C
root -b -l -q bdtout_syst.C
root -b -l -q idmva.C
root -b -l -q idmvaInputs_EB.C
root -b -l -q idmvaInputs_EE.C
root -b -l -q sigmaE.C
root -b -l -q 'diphomva_inputs_cat.C(0)'
root -b -l -q diphomva_inputs_cat.C
root -b -l -q ptRescale.C
root -b -l -q rhoVsNvtx.C
root -b -l -q 'zee_mass.C(0)'
root -b -l -q zee_mass.C
