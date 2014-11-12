import FWCore.ParameterSet.Config as cms

from RecoTracker.TkSeedingLayers.seedingLayersEDProducer_cfi import *

CAPentupletsA = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+TIB1',
			'BPix3+TIB1+TIB2'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TIB = cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit")
)
)

CAPentupletsB_pos = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+FPix1_pos',
			'BPix2+FPix1_pos+TID1_pos',
			'FPix1_pos+TID1_pos+TID2_pos'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
FPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(1),
)
)

CAPentupletsB_neg = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+FPix1_neg',
			'BPix2+FPix1_neg+TID1_neg',
			'FPix1_neg+TID1_neg+TID2_neg'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
FPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(1)  
)
)

CAPentupletsC_pos = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+TIB1',
			'BPix3+TIB1+TID1_pos'
			
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TIB = cms.PSet(
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

CAPentupletsC_neg = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+TIB1',
			'BPix3+TIB1+TID1_neg',
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TIB = cms.PSet(
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


CAPentupletsB = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+FPix1_pos',
			'BPix1+BPix2+FPix1_neg',
			'BPix2+FPix1_pos+TID1_pos',
			'BPix2+FPix1_neg+TID1_neg',
			'FPix1_pos+TID1_pos+TID2_pos',
			'FPix1_neg+TID1_neg+TID2_neg'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
FPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(1),
)
)

#brand new BD
CAPentupletsBD = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix1+BPix2+FPix1_pos',
			'BPix1+BPix2+FPix1_neg',
			'BPix2+BPix3+FPix1_pos',
			'BPix2+BPix3+FPix1_neg',
			'BPix2+FPix1_pos+TID1_pos',
			'BPix2+FPix1_neg+TID1_neg',
			'BPix3+FPix1_pos+TID1_pos',
			'BPix3+FPix1_neg+TID1_neg',
			'FPix1_pos+TID1_pos+TID2_pos',
			'FPix1_neg+TID1_neg+TID2_neg'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
FPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(1),
)
)


#brand new D
CAPentupletsD = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+FPix1_pos',
			'BPix2+BPix3+FPix1_neg',
			'BPix3+FPix1_pos+TID1_pos',
			'BPix3+FPix1_neg+TID1_neg'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
FPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(1),
)
)


CAPentupletsD_pos = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+FPix1_pos',
			'BPix3+FPix1_pos+TID1_pos'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
FPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(1),
)
)

CAPentupletsD_neg = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+FPix1_neg',
			'BPix3+FPix1_neg+TID1_neg'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
FPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(1),
)
)


CAPentupletsC = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+TIB1',
			'BPix3+TIB1+TID1_neg',
			'BPix3+TIB1+TID1_pos'
),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
TIB = cms.PSet(
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



CAPentupletsAL = seedingLayersEDProducer.clone(
layerList = cms.vstring('BPix1+BPix2+BPix3',
			'BPix2+BPix3+TIB1',
			'BPix3+TIB1+TIB2',
			'BPix3+CTIB1+TID1_pos',
	        'BPix3+CTIB1+TID1_neg'
			),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
CTIB =  cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        MinAbsZ = cms.double(45.0)
),
TIB = cms.PSet(
       TTRHBuilder = cms.string('WithTrackAngle'),
       matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
       #MinAbsZ = cms.double(45.0)
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(2)
)
)

#A+C+D layers
CAPentupletsAllLay = seedingLayersEDProducer.clone(
layerList = cms.vstring(
			'BPix1+BPix2+BPix3',
			'BPix2+BPix3+TIB1',
			'BPix2+BPix3+FPix1_pos',
			'BPix2+BPix3+FPix1_neg',
			'BPix3+TIB1+TIB2',
			'BPix3+CTIB1+CTID1_pos',
		    'BPix3+CTIB1+CTID1_neg',
			'BPix3+FPix1_pos+TID1_pos',
			'BPix3+FPix1_neg+TID1_neg',
			#try to add also the B layers ... 
			'BPix1+BPix2+FPix1_pos',
             'BPix1+BPix2+FPix1_neg',
             'BPix2+FPix1_pos+TID1_pos',
             'BPix2+FPix1_neg+TID1_neg',
             'FPix1_pos+TID1_pos+TID2_pos',
             'FPix1_neg+TID1_neg+TID2_neg'
			),
BPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
FPix = cms.PSet(
        TTRHBuilder = cms.string('TTRHBuilderWithoutAngle4PixelTriplets'),
        HitProducer = cms.string('siPixelRecHits')
),
CTIB =  cms.PSet(
        TTRHBuilder = cms.string('WithTrackAngle'),
        matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
        MinAbsZ = cms.double(45.0)
),
TIB = cms.PSet(
       TTRHBuilder = cms.string('WithTrackAngle'),
       matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
),
CTID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(2)
),
TID = cms.PSet(
         matchedRecHits = cms.InputTag("siStripMatchedRecHits","matchedRecHit"),
         useRingSlector = cms.bool(True),
         TTRHBuilder = cms.string('WithTrackAngle'),
         maxRing = cms.int32(2),
         minRing = cms.int32(1)
)
)


