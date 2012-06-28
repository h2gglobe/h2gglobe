CXXFLAGS += $(shell $$ROOTSYS/bin/root-config --cflags)
LDLIBS += $(shell $$ROOTSYS/bin/root-config --glibs)

#----------------------------------------
# boost
#----------------------------------------
CXXFLAGS += -I$(BOOST_BASE)/include
LDFLAGS += -L$(BOOST_BASE)/lib

LDLIBS +=  -lboost_regex

#----------------------------------------
CC = $(CXX)



plotter:		plotter.o HistoContainer.o 



clean:
			-rm -f *.o
			-rm -f plotter
			-rm -f *~
