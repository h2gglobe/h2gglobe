## mkdir logs

### ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 3vars    -T pervtx.json               2>&1 | tee logs/3vars.log &
### ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 5vars    -T pervtx.json,5vars.json    2>&1 | tee logs/5vars.log & 
### ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 3vars_pt -T pervtx.json,pt.json       2>&1 | tee logs/3vars_pt.log &
### ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 5vars_pt -T pervtx.json,5vars_pt.json 2>&1 | tee logs/5vars_pt.log &
### ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -O -l 5varsOpt -T pervtx.json,5varsOpt.json 2>&1 | tee logs/5varsOpt.log &
### ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 5varsNoConv -T pervtx.json,5varsOpt.json 2>&1 | tee logs/5varsNoConv.log &
./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 3varsNoConv -T pervtx.json,3varsOpt.json 2>&1 | tee logs/3varsNoConv.log &
## ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 6varsNoConv -T pervtx.json,6varsOpt.json 2>&1 | tee logs/6varsNoConv.log &

### ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 5varsNoConvKS -T pervtx.json,5varsOptKS.json 2>&1 | tee logs/5varsNoConvKS.log & 

./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 3varsConvNlegs -T pervtx.json,3varsConvNlegs.json  2>&1 | tee logs/3varsConvNlegs.log &
#### 
#### ./classify.py -i ../../AnalysisScripts/ -f vtxOptReduction_v2.root -l 3varsConv -T pervtx.json,3varsConv.json  2>&1 | tee logs/3varsConv.log &

wait

