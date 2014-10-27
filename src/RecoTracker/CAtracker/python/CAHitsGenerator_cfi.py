import FWCore.ParameterSet.Config as cms

from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
#from RecoTracker.TkTrackingRegions.GlobalTrackingRegion_cfi import *
#from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import *


CAHitsGenerator = cms.PSet(
 	debug = cms.untracked.int32(0),
 	maketrees = cms.untracked.bool(False),
   	Builder = cms.untracked.string('WithAngleAndTemplate'),
   	EtaCut = cms.double(0.0256),
   	#ComponentName = cms.string('CAHitsGenerator'),
   	GeneratorPSet = cms.PSet(PixelTripletHLTGenerator)
	#GeneratorPSet = cms.PSet(PixelTripletHLTGeneratorWithFilter)
)


