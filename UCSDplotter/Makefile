CXXFLAGS += $(shell $$ROOTSYS/bin/root-config --cflags)
LDLIBS += $(shell $$ROOTSYS/bin/root-config --glibs)

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

plotter:		plotter.o HistoContainer.o 


JimAnalysis.so:		JimAnalysis.o 
			$(LINK.cc) $^ -shared -o $@ 


clean:
			-rm -f *.o
			-rm -f $(ALL_TARGETS)
			-rm -f *~
