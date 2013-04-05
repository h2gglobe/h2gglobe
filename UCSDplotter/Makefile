ROOFIT_BASE=$(ROOFITSYS)
CXXFLAGS += $(shell $$ROOTSYS/bin/root-config --cflags) -I$(ROOFIT_BASE)/include
LDLIBS += $(shell $$ROOTSYS/bin/root-config --glibs)
LDFLAGS += -L$(ROOFIT_BASE)/lib -L$(ROOTSYS)/lib -lRooFitCore -lRooFit -lTMVA 

# for building the shared libraries
CXXFLAGS += -fPIC

# include debugging information
CXXFLAGS += -g

#----------------------------------------
# boost
#----------------------------------------
CXXFLAGS += -I$(BOOST_BASE)/include
LDFLAGS += -L$(BOOST_BASE)/lib

LDLIBS +=  -lboost_regex

#----------------------------------------
CC = $(CXX)

ALL_TARGETS = plotter JimAnalysis.so 

all:			$(ALL_TARGETS)

plotter:		plotter.o parser.o HistoContainer.o ../RooContainer.o Plotting.o ../BaseSmearer.o ../EnergySmearer.o \
			../EfficiencySmearer.o	../DiPhoEfficiencySmearer.o ../KFactorSmearer.o ../InterferenceSmearer.o


JimAnalysis.so:		JimAnalysis.o 
			$(LINK.cc) $^ -shared -o $@ 


clean:
			-rm -f *.o
			-rm -f $(ALL_TARGETS)
			-rm -f *~
