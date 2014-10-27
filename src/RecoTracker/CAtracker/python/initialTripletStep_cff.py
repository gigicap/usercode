import FWCore.ParameterSet.Config as cms

### STEP 00 ###

initialTripletStepClusters = cms.EDProducer("TrackClusterRemover",
    clusterLessSolution = cms.bool(True),
    trajectories = cms.InputTag("initialStepTracksCA"),
    overrideTrkQuals = cms.InputTag('initialStepCA'),
    TrackQuality = cms.string('highPurity'),
    minNumberOfLayersWithMeasBeforeFiltering = cms.int32(0),
    pixelClusters = cms.InputTag('siPixelClusters'),
    stripClusters = cms.InputTag('siStripClusters'),
    Common = cms.PSet(
        maxChi2 = cms.double(9.0),
    )
)


# SEEDING LAYERS
import RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi
initialTripletStepSeedLayers = RecoTracker.TkSeedingLayers.PixelLayerTriplets_cfi.PixelLayerTriplets.clone()

#can be used only in forward... 
initialTripletStepSeedLayers.layerList = cms.vstring(
    'BPix1+BPix2+FPix1_pos', 
    'BPix1+BPix2+FPix1_neg', 
    'BPix1+FPix1_pos+FPix2_pos', 
    'BPix1+FPix1_neg+FPix2_neg'
)

initialTripletStepSeedLayers.BPix.skipClusters = cms.InputTag('initialTripletStepClusters')
initialTripletStepSeedLayers.FPix.skipClusters = cms.InputTag('initialTripletStepClusters')


# seeding
from RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff import *
from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock
initialTripletStepSeeds = RecoTracker.TkSeedGenerator.GlobalSeedsFromTriplets_cff.globalSeedsFromTriplets.clone(
    RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
    ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
    RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
    ptMin = 0.6,
    originRadius = 0.02,
    nSigmaZ = 4.0
    )
    )
    )
initialTripletStepSeeds.OrderedHitsFactoryPSet.SeedingLayers = 'initialTripletStepSeedLayers'

from RecoPixelVertexing.PixelLowPtUtilities.ClusterShapeHitFilterESProducer_cfi import *
import RecoPixelVertexing.PixelLowPtUtilities.LowPtClusterShapeSeedComparitor_cfi
initialTripletStepSeeds.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet = RecoPixelVertexing.PixelLowPtUtilities.LowPtClusterShapeSeedComparitor_cfi.LowPtClusterShapeSeedComparitor

# building
import TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff
initialTripletStepTrajectoryFilter = TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff.CkfBaseTrajectoryFilter_block.clone(
    minimumNumberOfHits = 3,
    minPt = 0.2
    )

import TrackingTools.KalmanUpdators.Chi2ChargeMeasurementEstimatorESProducer_cfi
initialTripletStepChi2Est = TrackingTools.KalmanUpdators.Chi2ChargeMeasurementEstimatorESProducer_cfi.Chi2ChargeMeasurementEstimator.clone(
    ComponentName = cms.string('initialTripletStepChi2Est'),
    nSigma = cms.double(3.0),
    MaxChi2 = cms.double(30.0),
    minGoodStripCharge = cms.double(1724),
    pTChargeCutThreshold = cms.double(15.)
)

import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi
initialTripletStepTrajectoryBuilder = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi.GroupedCkfTrajectoryBuilder.clone(
    trajectoryFilter = cms.PSet(refToPSet_ = cms.string('initialTripletStepTrajectoryFilter')),
    alwaysUseInvalidHits = True,
    maxCand = 3,
    estimator = cms.string('initialTripletStepChi2Est'),
    maxDPhiForLooperReconstruction = cms.double(2.0),
    maxPtForLooperReconstruction = cms.double(0.7)
    )

import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
initialTripletStepTrackCandidates = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
    src = cms.InputTag('initialTripletStepSeeds'),
    ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
    numHitsForSeedCleaner = cms.int32(50),
    onlyPixelHitsForSeedCleaner = cms.bool(True),
    TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('initialTripletStepTrajectoryBuilder')),
    doSeedingRegionRebuilding = True,
    useHitsSplitting = True
    )

# fitting
import RecoTracker.TrackProducer.TrackProducer_cfi
initialTripletStepTracks = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
    src = 'initialTripletStepTrackCandidates',
    AlgorithmName = cms.string('iter0'),
    Fitter = cms.string('FlexibleKFFittingSmoother')
    )

# Final selection
import RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi
from RecoTracker.IterativeTracking.DetachedTripletStep_cff import detachedTripletStepSelector
initialTripletStepSelector = RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.multiTrackSelector.clone(
    src='initialTripletStepTracks',
    useAnyMVA = cms.bool(True),
    GBRForestLabel = cms.string('MVASelectorIter0_13TeV_v0'),
    trackSelectors= cms.VPSet(
    RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
        name = 'initialTripletStepLoose',
        ), #end of pset
    RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.tightMTS.clone(
        name = 'initialTripletStepTight',
        preFilterName = 'initialTripletStepLoose',
        ),
    RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.highpurityMTS.clone(
        name = 'initialTripletStepV1',
        preFilterName = 'initialTripletStepTight',
        ),
    detachedTripletStepSelector.trackSelectors[4].clone(
        name = 'initialTripletStepV2',
        preFilterName=cms.string(''),
        keepAllTracks = cms.bool(False)
        ),
    detachedTripletStepSelector.trackSelectors[5].clone(
        name = 'initialTripletStepV3',
        preFilterName=cms.string(''),
        keepAllTracks = cms.bool(False)
        )
    ) #end of vpset
)#end of clone
import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
initialTripletStep = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = cms.VInputTag(cms.InputTag('initialTripletStepTracks'),
                                   cms.InputTag('initialTripletStepTracks'),
                                   cms.InputTag('initialTripletStepTracks')),
    hasSelector=cms.vint32(1,1,1),
    shareFrac = cms.double(0.99),
    indivShareFrac=cms.vdouble(1.0,1.0,1.0),
    selectedTrackQuals = cms.VInputTag(cms.InputTag("initialTripletStepSelector","initialTripletStepV1"),
                                       cms.InputTag("initialTripletStepSelector","initialTripletStepV2"),
                                       cms.InputTag("initialTripletStepSelector","initialTripletStepV3")),
    setsToMerge = cms.VPSet(cms.PSet( tLists=cms.vint32(0,1,2), pQual=cms.bool(True) )),
    writeOnlyTrkQuals=cms.bool(True)
    )

# Final sequence
InitialTripletStep = cms.Sequence(initialTripletStepClusters*initialTripletStepSeedLayers*
                           initialTripletStepSeeds*
                           initialTripletStepTrackCandidates*
                           initialTripletStepTracks*
                           initialTripletStepSelector*
                           initialTripletStep)