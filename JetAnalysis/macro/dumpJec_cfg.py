import FWCore.ParameterSet.Config as cms


process = cms.Process("jectxt")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')
options.register ('runOnMc',
                  False, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,          # string, int, or float
                  "runOnMc")

options.parseArguments()

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# define your favorite global tag
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

if(options.runOnMc) :
    gt="GR_P_V42_AN2::All"
else        :
    gt = "START53_V11::All"

process.GlobalTag.globaltag=cms.string(gt)

process.source = cms.Source("EmptySource")
process.readAK5PF    = cms.EDAnalyzer('JetCorrectorDBReader',  
                                      # below is the communication to the database 
                                      payloadName    = cms.untracked.string('AK5PF'),
                                      # this is used ONLY for the name of the printed txt files. You can use any name that you like, 
                                      # but it is recommended to use the GT name that you retrieved the files from.
                                      ## globalTag      = cms.untracked.string('GR_R_42_V23'),
                                      globalTag      = cms.untracked.string(gt),
                                      printScreen    = cms.untracked.bool(False),
                                      createTextFile = cms.untracked.bool(True)
                                      )
process.p = cms.Path(process.readAK5PF)
