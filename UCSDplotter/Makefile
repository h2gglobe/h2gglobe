#BOOST_BASE=$(BOOST_BASE)
ROOFIT_BASE=$(ROOFITSYS)

CXXFLAGS += $(shell $$ROOTSYS/bin/root-config --cflags) -I$(ROOFIT_BASE)/include
LDLIBS += $(shell $$ROOTSYS/bin/root-config --glibs)
LDFLAGS += -L $(ROOFIT_BASE)/lib -lRooFitCore -lRooFit -lTMVA -lCore

# for building the shared libraries
CXXFLAGS += -Wall -fPIC -g

#----------------------------------------
# boost
#----------------------------------------
CXXFLAGS += -I$(BOOST_BASE)/include
LDFLAGS += -L$(BOOST_BASE)/lib $(shell $$ROOTSYS/bin/root-config --glibs)

#LDLIBS += -lboost_regex-mt
LDLIBS += -lboost_regex

#----------------------------------------
CC = $(CXX)

.SUFFIXES: .o

ALL_TARGETS = workspacer MyAnalysis.so

all:			$(ALL_TARGETS)

workspacer:		workspacer.o parser.o ../RooContainer.o

MyAnalysis.so:		MyAnalysis.o
			$(LINK.cc) $^ -shared -o $@ 

clean:
			-rm -f *.o
			-rm -f $(ALL_TARGETS)
			-rm -f *~
