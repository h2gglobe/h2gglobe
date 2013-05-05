# name of the library
LIBNAME = SpinAnalysis

#Necessary to use shell built-in commands
SHELL=bash

$(shell mkdir -p bin)
$(shell mkdir -p obj)
$(shell mkdir -p lib)

# figure out os
UNAME := $(shell uname)
BASEDIR = $(shell pwd)

USERINCLUDES += -I$(ROOTSYS)/include
ifeq ($(UNAME),Darwin)
$(info Found OS Darwin - using mac osx build)
BOOST_PATH=/Applications/boost_1_49_0
USERINCLUDES += -I$(BOOST_PATH)/include/
USERLIBS += -L$(BOOST_PATH)/libs/ -lboost_regex -lboost_program_options
endif
ifeq ($(UNAME),Linux)
$(info Found OS Linux - using linux build)
BOOST_PATH=$(CMSSW_DATA_PATH)/../external/boost/1.47.0
USERINCLUDES += -I$(BOOST_PATH)/include/
USERINCLUDES += -I $(CMSSW_BASE)/src/ -I $(CMSSW_RELEASE_BASE)/src/
USERLIBS += -L$(BOOST_PATH)/lib/ -lboost_regex -lboost_program_options
USERINCLUDES += -I$(ROOFITSYS)/include
USERLIBS += -L$(ROOFITSYS)/lib -lRooFit -lRooFitCore
endif

USERLIBS += $(shell root-config --glibs) -lTMVA -lMinuit -lFoam -lRooFit -lRooFitCore -lMathMore
CXXFLAGS = -Wall -g -O2 -Wl,-rpath=$(BASEDIR)/lib 
LDFLAGS = -shared -W 

CXX=clang++
LD=clang++
CLANGPATH := $(shell type -p clang++)

ifeq ($(CLANGPATH),)
#$(warning clang++ not found, reverting to g++!)
CXX=g++
LD=g++
endif

INCLUDE +=  -I $(BASEDIR)/interface
INCLUDE += $(USERINCLUDES)
CXXFLAGS += $(INCLUDE)

LIBS += $(USERLIBS)

# this is where you  would normaly have the lib, bin and source directories
LIBDIR = $(BASEDIR)/lib
EXEDIR = $(BASEDIR)/bin
MACRODIR = $(BASEDIR)/src
SRCDIR = $(BASEDIR)/src
OBJDIR = $(BASEDIR)/obj
TESTDIR = $(BASEDIR)/test
DOCDIR= $(BASEDIR)/docs
OBJ_EXT=o
TEST_EXT=cpp

SRCS=$(wildcard $(BASEDIR)/src/*.cc)
EXES=$(wildcard $(BASEDIR)/test/*.cpp)
HEADERS=$(wildcard $(BASEDIR)/interface/*.h)
OBJS=$(subst $(SRCDIR), $(OBJDIR),$(subst cc,$(OBJ_EXT),$(SRCS)))
BINS=$(subst $(TESTDIR), $(EXEDIR),$(subst .$(TEST_EXT),,$(EXES)))

all:  $(BINS)

$(EXEDIR)/%: $(TESTDIR)/%.cpp $(LIBDIR)/lib$(LIBNAME).so $(HEADERS)
	@echo Compiling executable $@
	@$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME) $(USERLIBSX)


$(OBJDIR)/%.$(OBJ_EXT): $(SRCDIR)/%.cc
	@echo Making object $@
	@$(CXX) $(CXXFLAGS) -fPIC -c $<  -o $@

$(LIBDIR)/lib$(LIBNAME).so:$(OBJS)
	@echo Building shared library $@
	@$(LD) $(LDFLAGS) -o $(LIBDIR)/lib$(LIBNAME).so  $^ $(LIBS)
	

lib: $(LIBDIR)/lib$(LIBNAME).so

vars:
	@echo "LIBS: " $(LIBS)
	@echo "CXXFLAGS: " $(CXXFLAGS)
	@echo "Source files: " $(SRCS) 
	@echo "Object files: " $(OBJS)
	@echo "Executables:  " $(TARGETS)

clean:
	rm -rf $(OBJS) $(LIBDIR)/lib$(LIBNAME).so $(BINS)


