import FWCore.ParameterSet.Config as cms

#standard imports (from GlobalSeedsFromTriplets_cff.py)
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEParmError_cfi import *
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi import *
from RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi import *
from TrackingTools.MaterialEffects.MaterialPropagator_cfi import *

#import triplet producer an cahitsgenerator
from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
from RecoTracker.CAtracker.CAHitsGenerator_cfi import CAHitsGenerator as CellularAutomaton

import RecoTracker.CAtracker.BarrelPentuplets_cfi
PentupletLayers = RecoTracker.CAtracker.BarrelPentuplets_cfi.BarrelPentuplets.clone()

#globalSeedsFromMultiplets
import RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi
GlobalSeedsFromMultiplets = RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi.seedGeneratorFromRegionHitsEDProducer.clone(
OrderedHitsFactoryPSet = cms.PSet(
          CellularAutomaton,
          ComponentName = cms.string('CAHitsGenerator'), 
          SeedingLayers = cms.InputTag('PentupletLayers')
                                  #GeneratorPSet = cms.PSet(CellularAutomaton)                            
     )
)


import RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi
CATripletHitFilter  = RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi.ClusterShapeHitFilterESProducer.clone(
     ComponentName = cms.string('CATripletHitFilter'),
     PixelShapeFile= cms.string('RecoPixelVertexing/PixelLowPtUtilities/data/pixelShape.par'),
     minGoodStripCharge = cms.double(2069)
     )

GlobalSeedsFromMultiplets.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet = cms.PSet(
	 ComponentName = cms.string('PixelClusterShapeSeedComparitor'),
         FilterAtHelixStage = cms.bool(False),
         FilterPixelHits = cms.bool(True),
         FilterStripHits = cms.bool(True),
         ClusterShapeHitFilterName = cms.string('CATripletHitFilter'),
         ClusterShapeCacheSrc = cms.InputTag('siPixelClusterShapeCache')
     )
