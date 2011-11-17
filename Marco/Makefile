# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch

SrcSuf        = cc

#------------------------------------------------------------------------------


#------------------------------------------------------------------------------

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    

LOOPALL = LoopAll

LOOPALLSO = libLoopAll.$(DllSuf)

# integrate VertexAnalysis sub-package
VTX=VertexAnalysis
VTXSRC=$(wildcard $(VTX)/src/*.$(SrcSuf))
VTXOBS=$(patsubst %$(SrcSuf), %$(ObjSuf), $(VTXSRC))

PHO=PhotonAnalysis
PHOSRC=$(wildcard $(PHO)/src/*.$(SrcSuf))
PHOOBS=$(patsubst %$(SrcSuf), %$(ObjSuf), $(PHOSRC))

LOOPALLO = LoopAll.$(ObjSuf) \
	   LoopAllDict.$(ObjSuf) \
	   dict.$(ObjSuf) \
	   HistoContainer.o \
	   BaseAnalysis.o \
	   BaseSmearer.o \
	   EnergySmearer.o \
	   EfficiencySmearer.o \
	   DiPhoEfficiencySmearer.o \
	   KFactorSmearer.o \
	   CounterContainer.o \
	   SampleContainer.o \
	   PhotonReducedInfo.o \
	   RooContainer.o \
	   Cut.o \
	   TriggerSelection.o \
	   PhotonFix.o \
           $(VTXOBS) $(PHOOBS)

DICTS = LoopAll.h BaseAnalysis.h BaseSmearer.h EnergySmearer.h EfficiencySmearer.h DiPhoEfficiencySmearer.h EnergySmearer.h KFactorSmearer.h SampleContainer.h PhotonFix.h \
	VertexAnalysis/interface/VertexAlgoParameters.h\
	PhotonAnalysis/interface/PhotonAnalysis.h\
	PhotonAnalysis/interface/StatAnalysis.h\
	PhotonAnalysis/interface/StatAnalysisExclusive.h\
	RooContainer.h


ROOFIT_BASE=$(ROOFITSYS)
LDFLAGS+=-L$(ROOFIT_BASE)/lib $(ROOTLIBS) -lRooFitCore -lRooFit -lTMVA
CXXFLAGS+=-I$(ROOFIT_BASE)/include  -I$(CMSSW_BASE)/src

all: $(LOOPALL)

clean:
	@rm -f $(LOOPALLO) core *Dict.* *.so dict.cpp

.SUFFIXES: .$(SrcSuf)

CXXFLAGS+=-I$(shell pwd)

$(LOOPALL):  $(LOOPALLO)
	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTLIBS)  $(LOOPALLO) $(OutPutOpt) $(LOOPALLSO)
	@echo "$(LOOPALLSO) done"

LoopAll.$(ObjSuf): CommonParameters.h LoopAll.h Tools.h \
	branchdef/Limits.h branchdef/treedef.h branchdef/newclonesarray.h \
	branchdef/treebranch.h branchdef/setbranchaddress.h branchdef/getentry.h branchdef/getbranch.h branchdef/branchdef.h \
	GeneralFunctions_cc.h GeneralFunctions_h.h \
	BaseAnalysis.cc BaseAnalysis.h \
	BaseSmearer.cc BaseSmearer.h \
	EnergySmearer.cc EnergySmearer.h \
	KFactorSmearer.cc KFactorSmearer.h \
	PhotonReducedInfo.cc PhotonReducedInfo.h \
	HistoContainer.cc HistoContainer.h \
	CounterContainer.cc CounterContainer.h \
	SampleContainer.cc SampleContainer.h \
	RooContainer.cc RooContainer.h \
	TriggerSelection.h TriggerSelection.cc \
	PhotonFix.h PhotonFix.cc \
	Cut.cc Cut.h $(VTXSRC) $(PHOSRC)

LoopAllDict.$(SrcSuf): CommonParameters.h LoopAll.h \
	branchdef/Limits.h branchdef/treedef.h \
	GeneralFunctions_h.h \
	BaseAnalysis.h \
	BaseSmearer.h \
	EnergySmearer.h \
	HistoContainer.h \
	CounterContainer.h \
	SampleContainer.h \
	RooContainer.h \
	PhotonFix.h \
	Cut.h \
	VertexAnalysis/interface/VertexAlgoParameters.h PhotonAnalysis/interface/PhotonAnalysis.h PhotonAnalysis/interface/StatAnalysis.h PhotonAnalysis/interface/StatAnalysisExclusive.h

	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c -I$(ROOFIT_BASE)/include -I$(CMSSW_BASE)/src  $(DICTS)

dict.cpp:
	@rootcint -f dict.cpp -c -p LinkDef.h 

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -g -c $< -o $@

