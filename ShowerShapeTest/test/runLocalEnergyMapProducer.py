# Import configurations
import FWCore.ParameterSet.Config as cms


# set up process
process = cms.Process("EGAMMA")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis') 
options.register('isMiniAOD',True,options.multiplicity.singleton,options.varType.bool," whether we are running on MiniAOD or not")
options.parseArguments()


process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(options.inputFiles),  
                          )


# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(5000),
    limit = cms.untracked.int32(10000000)
)
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.load("Configuration.StandardSequences.Services_cff")

# set the number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(options.outputFile)
)

process.mapProducer = cms.EDAnalyzer("LocalEnergyMapProducer",
                                     fillFromEles = cms.bool(True),
                                     elesTag = cms.InputTag("slimmedElectrons"),
                                     phosTag = cms.InputTag("slimmedPhotons"),
                                     ebRecHitsTag = cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                     eeRecHitsTag = cms.InputTag("reducedEgamma","reducedEERecHits"),
                                     )

if not options.isMiniAOD:
    process.mapProducer.elesTag = cms.InputTag("gedGsfElectrons")
    process.mapProducer.phosTag = cms.InputTag("gedPhotons")
    process.mapProducer.ebRecHitsTag = cms.InputTag("reducedEcalRecHitsEB")
    process.mapProducer.eeRecHitsTag = cms.InputTag("reducedEcalRecHitsEE")


process.p = cms.Path(process.mapProducer)
