import FWCore.ParameterSet.Config as cms

process = cms.Process("DUMMYTESTLIGHT")

process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("Configuration.StandardSequences.Services_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Geometry_cff") 
#process.load("Configuration.StandardSequences.Reconstruction_cff")
#process.load("RecoTracker.CAtracker.CAGlobalTracking_cff")
process.load("RecoTracker.CAtracker.CAGlobalTracking00_cff")


# source
readFiles = cms.untracked.vstring()

source = cms.Source ("PoolSource",fileNames = readFiles)

###SingleMuPt1000
#readFiles.extend([ 'file:///afs/cern.ch/work/g/gigicap/TrakingPOG/2014/CMSSW_7_0_0_pre3_CellularAutomata/src/CAtracker/CAtracker/test/sampleFiles/SingleMuPt100_867A06F3-F098-E311-8118-003048FEB906.root'
#			]);

#ttbar
readFiles.extend( [
     '/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/42783A1F-1550-E411-B888-0025905B8582.root',
'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/A233AB64-0C50-E411-A954-0025905A60A0.root',
'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/C60B9C19-1550-E411-8517-002618FDA248.root',
'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/DCA629A7-0C50-E411-9AC6-002618943849.root',
'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/F20A822B-0B50-E411-99AC-0025905A60B4.root'
      ] );
#readFiles.extend( [
#       '/store/relval/CMSSW_7_1_0/RelValTTbar_13/GEN-SIM-RECO/PU25ns_POSTLS171_V15-v1/00000/12CAFE96-F9FE-E311-B68B-0025905964CC.root'
#                  ] )



process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'POSTLS171_V15::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')


process.evtInfo = cms.OutputModule("AsciiOutputModule")

process.localreco = cms.Sequence(process.siPixelRecHits*process.siStripMatchedRecHits)
#CAtracks
process.tracks = cms.Sequence(process.MeasurementTrackerEvent*process.offlineBeamSpot*process.CAglobaltracks)
#std tracks
#process.tracks = cms.Sequence(process.MeasurementTrackerEvent*process.offlineBeamSpot*process.iterTracking)


process.p = cms.Path(process.localreco*process.siPixelClusterShapeCache*process.standalonemuontracking*process.caloTowersRec*process.tracks)

process.ep = cms.EndPath(process.evtInfo)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
