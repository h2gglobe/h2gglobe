# name of the library
LIBNAME = FullMvaToolkit

#Necessary to use shell built-in commands
SHELL=bash

$(shell mkdir -p bin)
$(shell mkdir -p lib)
$(shell mkdir -p obj)
$(shell cp ../../libLoopAll.so lib/)

# figure out os
UNAME := $(shell uname)
BASEDIR = $(shell pwd)

USERINCLUDES += -I$(ROOTSYS)/include
BOOST_PATH=$(CMSSW_DATA_PATH)/../external/boost/1.47.0
USERINCLUDES += -I$(BOOST_PATH)/include/
USERINCLUDES += -I$(CMS_PATH)/$(SCRAM_ARCH)/lcg/roofit/5.32.00-cms5/include/
USERINCLUDES += -I $(CMSSW_BASE)/src/ -I $(CMSSW_RELEASE_BASE)/src/
USERLIBS += -L$(BOOST_PATH)/lib/ -lboost_regex -lboost_program_options
USERLIBS += -L$(ROOFITSYS)/lib/ -lRooFit -lRooFitCore
USERLIBS += -L$(BASEDIR)/../../ -lLoopAll

USERLIBS += $(shell root-config --glibs) -lTMVA -lMinuit -lFoam -lRooFit -lRooFitCore 
CXXFLAGS = -g -O2 -Wl,-rpath=$(BASEDIR)/lib 
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
	$(shell LD_LIBRARY_PATH+=:$(BASEDIR)/lib)
	$(shell LD_LIBRARY_PATH+=:$(BASEDIR)/../../)

$(EXEDIR)/%: $(TESTDIR)/%.cpp $(LIBDIR)/lib$(LIBNAME).so $(BASEDIR)/interface/*.h
	@echo Compiling executable $@
	@$(CXX) -o $@ $(CXXFLAGS) $< $(LIBS) -L$(LIBDIR) -l$(LIBNAME) $(USERLIBSX)


$(OBJDIR)/%.$(OBJ_EXT): $(SRCDIR)/%.cc
	@echo Making object $@
	@$(CXX) $(CXXFLAGS) -fPIC -c $<  -o $@ 

$(LIBDIR)/lib$(LIBNAME).so: $(OBJS)
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
	@rm -rf $(OBJS) $(LIBDIR)/lib$(LIBNAME).so $(BINS)
	@echo Cleaned


