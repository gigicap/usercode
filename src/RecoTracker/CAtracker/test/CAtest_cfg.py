import FWCore.ParameterSet.Config as cms
#import sys

process = cms.Process("CASEEDER")

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


# source
readFiles = cms.untracked.vstring()

source = cms.Source ("PoolSource",fileNames = readFiles)

###TTbar

###SingleMuPt1000
#readFiles.extend( [
#     '/store/relval/CMSSW_7_1_0/RelValSingleMuPt1000_UP15/GEN-SIM-RECO/POSTLS171_V15-v1/00000/F03572A3-B0FB-E311-BA24-0025905A6084.root' ] );
###TTbar+PU

#Local files
#mupt100
#readFiles.extend([ 'file:///afs/cern.ch/work/g/gigicap/TrakingPOG/2014/CMSSW_7_0_0_pre3_CellularAutomata/src/CAtracker/CAtracker/test/sampleFiles/SingleMuPt100_867A06F3-F098-E311-8118-003048FEB906.root'
#			]);
#ttbar
#readFiles.extend([ 
#'file:///afs/cern.ch/work/g/gigicap/TrakingPOG/2014/CMSSW_7_0_0_pre3_CellularAutomata/src/CAtracker/CAtracker/test/sampleFiles/TTbar13_6EA83C8D-C498-E311-B2A8-02163E00E993.root'
#			]);
			

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
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000))

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.GlobalTag.globaltag = 'POSTLS172_V1::All'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')


# hit building

process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
process.load("RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff")

#from RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff import *
#from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *

#triplets otf
process.load("RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi")

# seeding
process.load("RecoTracker.CAtracker.GlobalSeedsFromMultiplets_cff")


#from RecoTracker.CAtracker.GlobalSeedsFromMultiplets_cff import *
from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock

#from RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi import *
#from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock

#process.initialStepSeedsGigi = RecoTracker.CAtracker.GlobalSeedsFromMultiplets_cff.globalSeedsFromMultiplets.clone()



#process.CASeedingStep = RecoTracker.CAtracker.GlobalSeedsFromMultiplets_cff.GlobalSeedsFromMultiplets.clone(
process.CASeedingStep = process.GlobalSeedsFromMultiplets.clone(
  RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
       ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
       RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
       ptMin = 0.6,
       originRadius = 0.02,
       nSigmaZ = 4.0
       )
       )
)

#import RecoTracker.CAtracker.BarrelPentuplets_cfi
import RecoTracker.CAtracker.CAPentuplets_cfi
#import RecoTracker.CAtracker.CAPentupletsAllLay_cfi

process.PentupletLayers = RecoTracker.CAtracker.CAPentuplets_cfi.CAPentupletsAllLay.clone()
#process.PentupletLayers = RecoTracker.CAtracker.CAPentupletsAllLay_cfi.CAPentupletsAll.clone()

process.CASeedingStep.OrderedHitsFactoryPSet.SeedingLayers = 'PentupletLayers'

#external input 
#val_EtaCut=float(sys.argv[2])
#process.CASeedingStep.OrderedHitsFactoryPSet.EtaCut=cms.double(val_EtaCut)
process.CASeedingStep.OrderedHitsFactoryPSet.EtaCut=cms.double(0.0256)
process.CASeedingStep.OrderedHitsFactoryPSet.debug = cms.untracked.int32(2)
process.CASeedingStep.OrderedHitsFactoryPSet.maketrees = cms.untracked.bool(True)

process.evtInfo = cms.OutputModule("AsciiOutputModule")

process.TFileService = cms.Service("TFileService", fileName = cms.string("histotrip.root") )

process.s = cms.Sequence(process.siPixelRecHits*process.siStripMatchedRecHits*process.siPixelClusterShapeCache*process.PentupletLayers*process.CASeedingStep)
process.p = cms.Path(process.s)

process.ep = cms.EndPath(process.evtInfo)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
