import FWCore.ParameterSet.Config as cms

from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *

CAHitsGenerator = cms.PSet(
 	debug = cms.untracked.int32(0),
   	Builder = cms.untracked.string('WithAngleAndTemplate'),
   	EtaCut = cms.double(0.0256),
    MakeSecondBackward = cms.untracked.bool(True),
    MakeFinalForward = cms.untracked.bool(True),
    GeneratorPSet = cms.PSet(PixelTripletHLTGenerator)
)

CAGeneratorFromPairs = cms.PSet(
    debug = cms.untracked.int32(0),
    Builder = cms.untracked.string('WithAngleAndTemplate'),
    EtaCut = cms.double(0.0256),
    MakeSecondBackward = cms.untracked.bool(True),
    MakeFinalForward = cms.untracked.bool(True)
)


CAHitsGeneratorForDebuggers = cms.PSet(
    debug = cms.untracked.int32(1),
    Builder = cms.untracked.string('WithAngleAndTemplate'),
    EtaCut = cms.double(0.0256),
    maketrees = cms.untracked.bool(True),
    MakeSecondBackward = cms.bool(True),
    MakeFinalForward = cms.bool(True),
    GeneratorPSet = cms.PSet(PixelTripletHLTGenerator)
)
