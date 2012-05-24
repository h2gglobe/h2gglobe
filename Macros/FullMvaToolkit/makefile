CC = g++
ROOFITBASE=$(ROOFITSYS)
ROOFITLIBS=$(ROOFITBASE)/lib
ROOFITINCS=$(ROOFITBASE)/include
DEPENDENTS = `root-config --cflags --libs` -L$(ROOFITLIBS) -I$(ROOFITINCS) -lMinuit -lFoam -lRooFit -lRooFitCore -lboost_program_options
FILES = RunFullMvaAnalysis.cc src/FMTBase.cc src/FMTFit.cc src/FMTRebin.cc src/FMTSetup.cc src/FMTSigInterp.cc
OUT_EXE = runIt.exe

all: RunFullMvaAnalysis

RunFullMvaAnalysis: $(FILES)
	$(CC) -g -O2 -o $(OUT_EXE) $(FILES) $(DEPENDENTS)

clean:
	rm $(OUT_EXE)
