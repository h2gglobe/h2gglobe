#------------------------------------------------------------------------------
# setup for BOOST (take the CMSSW boost if defined)
#------------------------------------------------------------------------------

CXXFLAGS += -I$(BOOST_DIR)/include
LDFLAGS += -L$(BOOST_DIR)/lib

#------------------------------------------------------------------------------

ALL_EXECUTABLES = fitSignal

all:		$(ALL_EXECUTABLES)


fitSignal:	fitSignal.o XMLUtils.o ParametricSignalModelConfig.o
		g++ -L$$ROOFITSYS/lib -L$(BOOST_DIR)/lib -lRooFit -lRooFitCore -lFoam -lMinuit -lMathMore -lboost_regex -g  $^ -o $@ $$($$ROOTSYS/bin/root-config --glibs --cflags) -I$$ROOFITSYS/include -I$(BOOST_DIR)/include -lXMLParser 


%.o: 	%.C
		g++ -L$$ROOFITSYS/lib -lRooFit -lRooFitCore -lFoam -lMinuit -lMathMore -g -c $^ $$($$ROOTSYS/bin/root-config --glibs --cflags) -I$$ROOFITSYS/include -I$$BOOST_DIR/include -Wall

%.o: 	%.cc
		g++ -L$$ROOFITSYS/lib -lRooFit -lRooFitCore -lFoam -lMinuit -lMathMore -g -c $^ $$($$ROOTSYS/bin/root-config --glibs --cflags) -I$$ROOFITSYS/include -I$$BOOST_DIR/include -Wall


clean:
		-rm -f *.o
		-rm -f *.so
		-rm -f *~
		-rm -f $(ALL_EXECUTABLES)
