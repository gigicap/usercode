import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *
BarrelPentuplets = seedingLayersEDProducer.clone()

BarrelPentuplets.layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+TIB1',
			'BPix3+TIB1+TIB2'
)
#try to add TID
#BarrelPentuplets.layerList = cms.vstring('BPix1+BPix2+BPix3','BPix2+BPix3+TIB1','BPix3+TIB1+TIB2','BPix3+TIB1+TID1_pos','BPix3+TIB1+TID1_neg','BPix1+BPix2+FPix1_neg','BPix1+BPix2+FPix1_pos','FPix1_neg+TID1_neg+TID2_neg','FPix1_pos+TID1_pos+TID2_pos')


#BarrelPentuplets.layerList = cms.vstring('BPix1+BPix2+FPix1_pos','BPix2+FPix1_pos+TID1_pos','FPix1_pos+TID1_pos+TID2_pos')


BarrelPentuplets.BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
)

#BarrelPentuplets.FPix = cms.PSet(
#        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
#        HitProducer = cms.string('siPixelRecHits')
#)

BarrelPentuplets.TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
)

#BarrelPentuplets.TID = cms.PSet(
#         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
#         useRingSlector = cms.bool(True),
#         TTRHBuilder = cms.string('WithTrackAngle'),
#         maxRing = cms.int32(2),
#         minRing = cms.int32(1),
#)

