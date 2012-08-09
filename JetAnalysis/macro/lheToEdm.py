#!/bin/env cmsRun

##
## cmsRun configuration to convert LHE files to EDM
## 
## Creates LHEEventProduct collection starting from LHE file.
## It also converts LHE events to genParticles and make genJets.
## Optionally it runs the Higgs to gamma gamma globe ntuplizer
##
## run as:
## cmsRun lheToEdm.py maxEvents=-1 files=<list of input files> output=<output file>
##

##
## For instructions to download the Hgg analyzer see:
## https://twiki.cern.ch/twiki/bin/viewauth/CMS/Higgs2GAnalyzer#Checking_Out_H2gAnalyzer
## 

##
## Note: the HEPMC class has an annoying printout in the constructor
## To remove it
## addpkg SimDataFormats/GeneratorProducts 
## sed 's%cout << "Contructing HepMCProduct" << endl;%%' -i SimDataFormats/GeneratorProducts/src/HepMCProduct.cc
## scramb -j 10 SimDataFormats/GeneratorProducts
##

import FWCore.ParameterSet.Config as cms

process = cms.Process("LHE")

##
## Options through VarParsing
## https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideAboutPythonConfigFile#VarParsing_Example
## 
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('standard')
options.setDefault ('files', 'file:hgg_spin0_VBF.lhe')
options.setDefault ('output','')
options.setDefault ('maxEvents',10)

options.register ('runGlobe',
                  True, # default value
                  VarParsing.VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.VarParsing.varType.bool,          # string, int, or float
                  "runGlobe")

options.parseArguments()

output = options.output
firstInput = options.files[0].replace(".lhe","")
if output == "" or output == ".root":
    output = "%s.root" % firstInput
elif output.startswith("_"):
    output = "%s%s" % ( firstInput, output )
annotation = firstInput


##
## LHE source
## 
process.source = cms.Source("LHESource",
	fileNames = cms.untracked.vstring(options.files)
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('alpha'),
	name = cms.untracked.string('LHEF input'),
	annotation = cms.untracked.string(annotation)
)

## Control printouts
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

## Convert LHE to HepMC and then GenParticles
from GeneratorInterface.LHEInterface.lhe2HepMCConverter_cfi import lhe2HepMCConverter
process.generator = lhe2HepMCConverter.clone()

process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")

process.load("RecoJets.JetProducers.ak5GenJets_cfi") 
process.load("RecoJets.Configuration.GenJetParticles_cff")

process.p0 = cms.Path(
	process.generator *
	process.genParticles *
	process.genParticlesForJets *
        process.ak5GenJets
)

## Dump first events
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printList = cms.EDAnalyzer("ParticleListDrawer",
                                 src = cms.InputTag("genParticles"),
                                 maxEventsToPrint = cms.untracked.int32(5)
                                 )

process.p = cms.Path(
    process.printList
    )

## Globe ntuplizer
if options.runGlobe:
    process.load("HiggsAnalysis.HiggsTo2photons.h2ganalyzer_GEN_cfi")
    process.h2ganalyzer.RootFileName = output.replace(".root","_globe.root")
    process.p.insert(1,process.h2ganalyzer)

## Output 
process.LHE = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('keep *_source_*_*',
                                                                      'keep *_genParticles_*_*',
                                                                      'keep *_ak5GenJets_*_*'),
                               dataset = cms.untracked.PSet(dataTier = cms.untracked.string('LHE')),
                               fileName = cms.untracked.string(output)
)

## Done
process.outpath = cms.EndPath(process.LHE)
process.schedule = cms.Schedule(process.p0, process.p, process.outpath)
