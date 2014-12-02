import FWCore.ParameterSet.Config as cms
import sys

process = cms.Process("CASEEDER")

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


# source
readFiles = cms.untracked.vstring()

source = cms.Source ("PoolSource",fileNames = readFiles)

#SingleMuPt10
readFiles.extend( [
    '/store/relval/CMSSW_7_2_0_pre8/RelValSingleMuPt10_UP15/GEN-SIM-RECO/PRE_LS172_V15-v1/00000/480204D4-5650-E411-9A8E-0026189438ED.root',
    '/store/relval/CMSSW_7_2_0_pre8/RelValSingleMuPt10_UP15/GEN-SIM-RECO/PRE_LS172_V15-v1/00000/929F4678-5650-E411-9C6F-0025905B858C.root'] );

#ttbar
#readFiles.extend( [
#'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/42783A1F-1550-E411-B888-0025905B8582.root',
#'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/A233AB64-0C50-E411-A954-0025905A60A0.root',
#'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/C60B9C19-1550-E411-8517-002618FDA248.root',
#'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/DCA629A7-0C50-E411-9AC6-002618943849.root',
#'/store/relval/CMSSW_7_2_0_pre8/RelValTTbar_13/GEN-SIM-RECO/PU25ns_PRE_LS172_V15-v1/00000/F20A822B-0B50-E411-99AC-0025905A60B4.root'
#      ] );


process.source = source
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')


# hit building

process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
process.load("RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff")

#triplets otf
process.load("RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi")

# seeding
process.load("RecoTracker.CAtracker.GlobalSeedsFromMultiplets_cff")


#from RecoTracker.CAtracker.GlobalSeedsFromMultiplets_cff import *
from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock

process.CASeedingStep = process.GlobalSeedsFromMultipletsForDebuggers.clone(
  RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
       ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
       RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
       ptMin = 0.6,
       originRadius = 0.02,
       nSigmaZ = 4.0
       )
       )
)

import RecoTracker.CAtracker.CAPentuplets_cfi

process.PentupletLayers = RecoTracker.CAtracker.CAPentuplets_cfi.CAPentupletsAllLay.clone()

process.CASeedingStep.OrderedHitsFactoryPSet.SeedingLayers = 'PentupletLayers'

#external input 
val_EtaCut=float(sys.argv[2])
process.CASeedingStep.OrderedHitsFactoryPSet.EtaCut=cms.double(val_EtaCut)
val_Debug=int32(sys.argv[3])
process.CASeedingStep.OrderedHitsFactoryPSet.debug = cms.untracked.int32(val_Debug)
process.CASeedingStep.OrderedHitsFactoryPSet.maketrees = cms.untracked.bool(True)

process.evtInfo = cms.OutputModule("AsciiOutputModule")

val_OutName=sys.argv[4]
process.TFileService = cms.Service("TFileService", fileName = cms.string(val_OutName) )

process.s = cms.Sequence(process.siPixelRecHits*process.siStripMatchedRecHits*process.siPixelClusterShapeCache*process.PentupletLayers*process.CASeedingStep)
process.p = cms.Path(process.s)

process.ep = cms.EndPath(process.evtInfo)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
