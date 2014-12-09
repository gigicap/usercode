import FWCore.ParameterSet.Config as cms

#standard imports
from RecoLocalTracker.SiStripRecHitConverter.StripCPEfromTrackAngle_cfi import *
from RecoLocalTracker.SiStripRecHitConverter.SiStripRecHitMatcher_cfi import *
from RecoLocalTracker.SiPixelRecHits.PixelCPEParmError_cfi import *
from RecoTracker.TransientTrackingRecHit.TransientTrackingRecHitBuilder_cfi import *
from RecoTracker.MeasurementDet.MeasurementTrackerESProducer_cfi import *
from TrackingTools.MaterialEffects.MaterialPropagator_cfi import *

from RecoTracker.CAtracker.CAHitsGenerator_cfi import CAGeneratorFromPairs as CellularAutomaton

import RecoTracker.CAtracker.CAPentuplets_cfi
PentupletLayers = RecoTracker.CAtracker.CAPentuplets_cfi.CAPentupletsACPairs.clone()

#globalSeedsFromMultiplets
import RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi
GlobalSeedsFromCAPairs = RecoTracker.TkSeedGenerator.SeedGeneratorFromRegionHitsEDProducer_cfi.seedGeneratorFromRegionHitsEDProducer.clone(
OrderedHitsFactoryPSet = cms.PSet(
          CellularAutomaton,
          ComponentName = cms.string('CAHitsGenerator'), 
          SeedingLayers = cms.InputTag('PentupletLayers')
     )
)

