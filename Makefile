# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

include Makefile.arch

SrcSuf        = cc
HeadSuf       = h
ObjSuf        = o
DepSuf        = d

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) .$(DllSuf)
.PHONY:    

#------------------------------------------------------------------------------

LOOPALLSO = libLoopAll.$(DllSuf)

##
## Files in main directory
## 

## Headers
MainHead=$(wildcard *.$(HeadSuf)) 
MainHead+=$(wildcard branchdef/*.$(HeadSuf)) 

## Sources
MainSrc=$(filter-out dict.cc, $(filter-out LoopAllDict.cc,$(wildcard *.$(SrcSuf))) )
MainSrc+=LoopAllDict.cc dict.cc Macros/Normalization_8TeV.cc Macros/MassInterpolator.cc
MainObjs=$(patsubst %$(SrcSuf), %$(ObjSuf), $(MainSrc))

## ROOT dictionary
LinkDef=LinkDef.h 
MainDicts=LoopAll.h
MainDicts+=$(wildcard Base*.$(HeadSuf))
MainDicts+=$(wildcard *Smearer.$(HeadSuf))
MainDicts+=$(wildcard *Container.$(HeadSuf))
MainDicts+=PhotonFix.h MassResolution.h HtmlHelper.h Macros/Normalization_8TeV.h Macros/MassInterpolator.h


##
## Subdirectories
##
SubPkgs=PhotonAnalysis VertexAnalysis JetAnalysis PhotonJetAnalysis 
SubPkgsDict=VertexAnalysis/interface/VertexAlgoParameters.h

##
## Flags and external dependecies
## 
ROOFIT_BASE=$(ROOFITSYS)
LDFLAGS+=-L$(ROOFIT_BASE)/lib $(ROOTLIBS) -lRooFitCore -lRooFit -lTMVA 
LDFLAGS+= $(patsubst %, -L%, $(shell echo ${LD_LIBRARY_PATH} | tr ':' '\n')) -lFWCorePythonParameterSet -lFWCoreParameterSet -lCMGToolsExternal -lCondFormatsJetMETObjects
CXXFLAGS+=-I$(ROOFIT_BASE)/include -I$(CMSSW_BASE)/src  -I$(CMSSW_RELEASE_BASE)/src 
CXXFLAGS+= $(patsubst %, -I%, $(shell echo ${CMSSW_FWLITE_INCLUDE_PATH} | tr ':' '\n'))
CXXFLAGS+=-I$(shell pwd) -g

##
## Code from users
##
-include Makefile.user
SubPkgs     += $(UserPkgs)
SubPkgsDict += $(UserDict)

############################################################################################################################
###
### Do not modify below this point, unless you have a good reason
###
############################################################################################################################ 
ifneq ($(SubPkgs),)
SubPkgsHead=$(foreach Pack,$(SubPkgs),$(wildcard $(Pack)/interface/*.$(HeadSuf)))
SubPkgsSrc=$(foreach Pack,$(SubPkgs),$(wildcard $(Pack)/src/*.$(SrcSuf)))
SubPkgsObjs=$(patsubst %$(SrcSuf), %$(ObjSuf), $(SubPkgsSrc))
_SubPkgsDict=$(foreach Pack,$(SubPkgs),$(wildcard $(Pack)/interface/*Analysis.$(HeadSuf))) $(SubPkgsDict)
else
SubPkgsHead=
SubPkgsSrc=
SubPkgsDict=
_SubPkgsDict=
endif

## 
Objs = $(MainObjs) $(SubPkgsObjs)
Dicts = $(MainDicts) $(_SubPkgsDict) 
Deps = $(patsubst %$(ObjSuf), %$(DepSuf), $(Objs))
ExtPacks=.extraTags

## Targets
all: $(ExtPacks)
	@$(MAKE)  $(LOOPALLSO)

print:
	@echo "Subpackages:"
	@echo "------------"
	@echo $(SubPkgs)
	@echo

	@echo "Sources:"
	@echo "--------"	
	@echo $(MainSrc) | tr ' ' '\n'
	@echo $(SubPkgsSrc)  | tr ' ' '\n'
	@echo 

	@echo "Dictionary sources:"
	@echo "-------------------"
	@echo $(MainDicts)  | tr ' ' '\n'
	@echo $(_SubPkgsDict)  | tr ' ' '\n'
	@echo 

	@echo "CXXFLAGS: "
	@echo "-------------------"	
	@echo "$(CXXFLAGS) " | tr ' ' '\n'
	@echo 

	@echo "LDFLAGS: "
	@echo "-------------------"
	@echo "$(LDFLAGS)" | tr ' ' '\n'
	@echo

clean:
	@rm -fv $(Objs) $(Deps) $(LOOPALL) *[dD]ict.* .extraTags

.extraTags: extraTags
	@echo "Getting extra tags"
	@bash extraTags

$(LOOPALLSO):  $(Objs)
	@echo "Linking"
	@$(LD) $(SOFLAGS) $(LDFLAGS) $(ROOTLIBS)  $(Objs) $(OutPutOpt) $(LOOPALLSO)
	@echo "$(LOOPALLSO) done"

LoopAllDict.$(SrcSuf): $(MainHead) $(SubPkgsHead)
	@echo "Generating dictionary $@"
	@rootcint -v4 -f $@ -c -I$(ROOFIT_BASE)/include -I$(CMSSW_BASE)/src  -I$(CMSSW_RELEASE_BASE)/src $(Dicts)

dict.$(SrcSuf): $(LinkDef)
	@echo "Generating dictionary $@"
	@rootcint -f dict.cc -c -p $(LinkDef)

%.$(ObjSuf): $(ExtPacks)

.$(SrcSuf).$(ObjSuf): $(ExtPacks)
	@echo "Compiling $<"
	@$(CXX) $(CXXFLAGS) -M -c $< -o $(patsubst %.$(ObjSuf), %.$(DepSuf), $@)
	@sed -i "s|.*:|$*.o: Makefile $(ExtPacks)|"  $(patsubst %.$(ObjSuf), %.$(DepSuf), $@)
	@$(CXX) $(CXXFLAGS) -g -c $< -o $@

-include $(Deps)

