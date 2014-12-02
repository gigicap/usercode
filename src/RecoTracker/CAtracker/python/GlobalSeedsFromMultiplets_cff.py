import FWCore.ParameterSet.Config as cms

#standard imports
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEParmError_cfi import *
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi import *
from RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi import *
from TrackingTools.MaterialEffects.MaterialPropagator_cfi import *

#import triplet producer
#from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
from RecoTracker.CAtracker.CAHitsGenerator_cfi import CAHitsGenerator as CellularAutomaton

import RecoTracker.CAtracker.CAPentuplets_cfi
PentupletLayers = RecoTracker.CAtracker.CAPentuplets_cfi.CAPentupletsAllLay.clone()

#globalSeedsFromMultiplets
import RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi
GlobalSeedsFromMultiplets = RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi.seedGeneratorFromRegionHitsEDProducer.clone(
OrderedHitsFactoryPSet = cms.PSet(
          CellularAutomaton,
          ComponentName = cms.string('CAHitsGenerator'), 
          SeedingLayers = cms.InputTag('PentupletLayers')
     )
)


import RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi

GlobalSeedsFromMultiplets.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet = cms.PSet(
 ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
         FilterAtHelixStage = cms.bool(False),
         FilterPixelHits = cms.bool(True),
         FilterStripHits = cms.bool(True),
         ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
         ClusterShapeCacheSrc = cms.InputTag('siPixelClusterShapeCache')
     )


#fordebuggers
from RecoTracker.CAtracker.CAHitsGenerator_cfi import CAHitsGeneratorForDebuggers as CellularAutomatonForDebuggers

GlobalSeedsFromMultipletsForDebuggers = RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi.seedGeneratorFromRegionHitsEDProducer.clone(
    OrderedHitsFactoryPSet = cms.PSet(
          CellularAutomatonForDebuggers,
          ComponentName = cms.string('CAHitsGeneratorForDebuggers'),
          SeedingLayers = cms.InputTag('PentupletLayers')
    )
 )

GlobalSeedsFromMultipletsForDebuggers.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet = cms.PSet(
       ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
       FilterAtHelixStage = cms.bool(False),
       FilterPixelHits = cms.bool(True),
       FilterStripHits = cms.bool(True),
       ClusterShapeHitFilterName = cms.string('ClusterShapeHitFilter'),
       ClusterShapeCacheSrc = cms.InputTag('siPixelClusterShapeCache')
)

