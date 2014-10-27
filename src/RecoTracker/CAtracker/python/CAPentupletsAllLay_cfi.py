import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *

CAPentuplets = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+TIB1',
			#'BPix3+TIB1+TIB2',
			#'BPix2+BPix3+TIBC1',
			#'BPix3+TIBC1+TID1_pos',
			'BPix3+TIBC1+TID1_neg'
			),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
),
TIBC =  cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        MinAbsZ = cms.double(45.0)
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(2)
)
)




