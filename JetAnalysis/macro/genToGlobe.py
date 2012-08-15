#!/bin/env cmsRun

##
## cmsRun configuration to create globe ntuples from GEN
## 
## run as:
## cmsRun genToGlobe.py maxEvents=-1 files=<list of input files> output=<output file>
##

##
## For instructions to download the Hgg analyzer see:
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/Higgs2GAnalyzer#Checking_Out_H2gAnalyzer
## 

import FWCore.ParameterSet.Config as cms

process = cms.Process("Globe")

##
## Options through VarParsing
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#VarParsing_Example
## 

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')
options.setDefault ('files', 'file:hgg_spin0_VBF.root')
options.setDefault ('output','')
options.setDefault ('maxEvents',10)

options.register ('runGlobe',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,          # string, int, or float
                  "runGlobe")

options.parseArguments()

output = options.output
firstInput = options.files[0].replace(".root","")
if output == "" or output == ".root":
    output = "%s.root" % firstInput
elif output.startswith("_"):
    output = "%s%s" % ( firstInput, output )
annotation = firstInput

print output

## Source
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring(options.files)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

## Control printouts
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

## Globe ntuplizer
process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_GEN_cfi")
process.h2ganalyzer.RootFileName = output.replace(".root","_globe.root")

## Output
process.p = cms.Path(
    process.h2ganalyzer
    )
