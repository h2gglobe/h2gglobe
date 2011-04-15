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

LOOPALLO = LoopAll.$(ObjSuf) \
           Util.$(ObjSuf) \
	   LoopAllDict.$(ObjSuf) \
	   dict.$(ObjSuf) \
	   HistoContainer.o \
	   CounterContainer.o \
	   SampleContainer.o \
	   RooContainer.o \
	   Cut.o \
           $(VTXOBS)

all: $(LOOPALL)

clean:
	@rm -f $(LOOPALLO) core *Dict.* *.so

.SUFFIXES: .$(SrcSuf)


$(LOOPALL):  $(LOOPALLO)
	$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTLIBS)  $(LOOPALLO) $(OutPutOpt) $(LOOPALLSO)
	@echo "$(LOOPALLSO) done"

LoopAll.$(ObjSuf): CommonParameters.h LoopAll.h Util.h Tools.h \
	branchdef/Limits.h branchdef/treedef.h branchdef/newclonesarray.h \
	branchdef/treebranch.h branchdef/setbranchaddress.h branchdef/getentry.h branchdef/getbranch.h branchdef/branchdef.h \
	PhotonAnalysis/PhotonAnalysisFunctions_h.h PhotonAnalysis/PhotonAnalysisFunctions_cc.h \
	GeneralFunctions_cc.h GeneralFunctions_h.h \
	HistoContainer.cc HistoContainer.h \
	CounterContainer.cc CounterContainer.h \
	SampleContainer.cc SampleContainer.h \
	RooContainer.cc RooContainer.h \
	Cut.cc Cut.h $(VTXSRC)

mpUtil.$(ObjSuf): CommonParameters.h LoopAll.h Util.h \
	branchdef/Limits.h branchdef/treedef.h branchdef/newclonesarray.h \
	branchdef/treebranch.h branchdef/setbranchaddress.h branchdef/getentry.h branchdef/getbranch.h branchdef/branchdef.h \
	PhotonAnalysis/PhotonAnalysisFunctions_h.h \
	GeneralFunctions_h.h \
	HistoContainer.h \
	CounterContainer.h \
	SampleContainer.h \
	RooContainer.h \
	Cut.h

LoopAllDict.$(SrcSuf): CommonParameters.h LoopAll.h Util.h \
	branchdef/Limits.h branchdef/treedef.h \
	PhotonAnalysis/PhotonAnalysisFunctions_h.h \
	GeneralFunctions_h.h \
	HistoContainer.h \
	CounterContainer.h \
	SampleContainer.h \
	RooContainer.h \
	Cut.h

	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c LoopAll.h Util.h 
	@rootcint -f dict.cpp -c -p LinkDef.h 

.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -g -c $< -o $@

