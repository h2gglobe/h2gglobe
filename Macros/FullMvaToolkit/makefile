CC = g++
ROOFITBASE=$(ROOFITSYS)
ROOFITLIBS=$(ROOFITBASE)/lib
ROOFITINCS=$(ROOFITBASE)/include
#BOOSTBASE=$(CMSSW_DATA_PATH)/../external/boost/1.47.0
BOOSTBASE=$(shell scram tool tag boost BOOST_BASE)
BOOSTLIBS=$(BOOSTBASE)/lib
BOOSTINCS=$(BOOSTBASE)/include

DEPENDENTS = `root-config --cflags --libs` -L$(ROOFITLIBS) -I$(ROOFITINCS) -lMinuit -lFoam -lRooFit -lRooFitCore -I$(BOOSTINCS) -L$(BOOSTLIBS) -lboost_program_options 
FILES = RunFullMvaAnalysis.cc src/FMTBase.cc src/FMTFit.cc src/FMTRebin.cc src/FMTSetup.cc src/FMTSigInterp.cc src/FMTPlots.cc ../Normalization_8TeV.cc
OUT_EXE = runIt.exe

all: RunFullMvaAnalysis

RunFullMvaAnalysis: $(FILES)
	$(CC) -g -O2 -o $(OUT_EXE) $(FILES) $(DEPENDENTS)

clean:
	rm $(OUT_EXE)
