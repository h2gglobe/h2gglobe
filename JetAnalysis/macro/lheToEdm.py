#!/bin/env cmsRun
import FWCore.ParameterSet.Config as cms

## process = cms.Process("Gen")
process = cms.Process("LHE")

process.source = cms.Source("LHESource",
	fileNames = cms.untracked.vstring('file:hgg_spin0_VBF.lhe')
)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))

process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('alpha'),
	name = cms.untracked.string('LHEF input'),
	annotation = cms.untracked.string('hgg spin0 VBF')
)

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'

## Options and Output Report
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
# reduce verbosity
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

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

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printList = cms.EDAnalyzer("ParticleListDrawer",
                                 src = cms.InputTag("genParticles"),
                                 maxEventsToPrint = cms.untracked.int32(5)
                                 )

process.p = cms.Path(
    process.printList
    )

process.LHE = cms.OutputModule("PoolOutputModule",
                               outputCommands = cms.untracked.vstring('keep *_source_*_*',
                                                                      'keep *_genParticles_*_*',
                                                                      'keep *_ak5GenJets_*_*'),
                               dataset = cms.untracked.PSet(dataTier = cms.untracked.string('LHE')),
                               fileName = cms.untracked.string('LHE.root')
)

process.outpath = cms.EndPath(process.LHE)

process.schedule = cms.Schedule(process.p0, process.p, process.outpath)
