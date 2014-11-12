import FWCore.ParameterSet.Config as cms
### Build seeds from CA pentuplets (Gigi's Code)### 
### STEP 0 ###
 
# hit building
from RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff import *
from RecoTracker.TransientTrackingRecHit.TTRHBuilders_cff import *

#import triplet producer an cahitsgenerator
#from RecoPixelVertexing.PixelTriplets.PixelTripletHLTGenerator_cfi import *
#from RecoTracker.CAtracker.CAHitsGenerator_cfi import CAHitsGenerator as CellularAutomaton

#import seeding layers
#import RecoTracker.CAtracker.BarrelPentuplets_cfi
#import RecoTracker.CAtracker.CAPentupletsAllLay_cfi
import RecoTracker.CAtracker.CAPentuplets_cfi
PentupletLayers = RecoTracker.CAtracker.CAPentuplets_cfi.CAPentupletsAllLay.clone()

from RecoTracker.TkTrackingRegions.GlobalTrackingRegionFromBeamSpot_cfi import RegionPsetFomBeamSpotBlock
from RecoTracker.CAtracker.GlobalSeedsFromMultiplets_cff import *

initialStepSeedsCA = RecoTracker.CAtracker.GlobalSeedsFromMultiplets_cff.GlobalSeedsFromMultiplets.clone(
RegionFactoryPSet = RegionPsetFomBeamSpotBlock.clone(
       ComponentName = cms.string('GlobalRegionProducerFromBeamSpot'),
       RegionPSet = RegionPsetFomBeamSpotBlock.RegionPSet.clone(
       		ptMin = 0.6,
       		originRadius = 0.02,
       		nSigmaZ = 4.0
       		)
       )
)
initialStepSeedsCA.OrderedHitsFactoryPSet.SeedingLayers = 'PentupletLayers'

#from RecoPixelVertexing.PixelLowPtUtilities.ClusterHitFilterESProducer_cfi import *
#import RecoPixelVertexing.PixelLowPtUtilities.LowPtClusterShapeSeedComparitor_cfi
#initialStepSeedsCA.OrderedHitsFactoryPSet.GeneratorPSet.SeedComparitorPSet = RecoPixelVertexing.PixelLowPtUtilities.LowPtClusterShapeSeedComparitor_cfi.LowPtClusterShapeSeedComparitor   #?


# building  
#from TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff import CkfBaseTrajectoryFilter_block
#import RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi
#MeasurementTrackerEvent = RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi.MeasurementTrackerEvent.clone()
import TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff
initialStepTrajectoryFilterCA = TrackingTools.TrajectoryFiltering.TrajectoryFilter_cff.CkfBaseTrajectoryFilter_block.clone(
    minimumNumberOfHits = 5, #???
    minPt = 0.2
) 
#initialStepTrajectoryFilterCA = CkfBaseTrajectoryFilter_block.clone(
#	minimumNumberOfHits = 5, #???
#    minPt = 0.2
#)
  
import TrackingTools.KalmanUpdators.Chi2ChargeMeasurementEstimatorESProducer_cfi
initialStepChi2Est = TrackingTools.KalmanUpdators.Chi2ChargeMeasurementEstimatorESProducer_cfi.Chi2ChargeMeasurementEstimator.clone(
      ComponentName = cms.string('initialStepChi2Est'),
      nSigma = cms.double(3.0),
      MaxChi2 = cms.double(30.0),
      minGoodStripCharge = cms.double(1724),
      pTChargeCutThreshold = cms.double(15.)
)


from RecoTracker.CkfPattern.CkfTrackCandidates_cff import *

import RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi
initialStepTrajectoryBuilderCA = RecoTracker.CkfPattern.GroupedCkfTrajectoryBuilder_cfi.GroupedCkfTrajectoryBuilder.clone(
      trajectoryFilter = cms.PSet(refToPSet_ = cms.string('initialStepTrajectoryFilterCA')),
      alwaysUseInvalidHits = True,
      maxCand = 3,
      estimator = cms.string('initialStepChi2Est'),
      maxDPhiForLooperReconstruction = cms.double(2.0),
      maxPtForLooperReconstruction = cms.double(0.7)
)

#navigation school
#from RecoTracker.TkNavigation.NavigationSchoolESProducer_cff import *
#AnyDirectionAnalyticalPropagator = cms.ESProducer( "AnalyticalPropagatorESProducer",
#  MaxDPhi = cms.double( 1.6 ),
#  ComponentName = cms.string( "AnyDirectionAnalyticalPropagator" ),
#  PropagationDirection = cms.string( "anyDirection" )
#)

import RecoTracker.CkfPattern.CkfTrackCandidates_cfi
initialStepTrackCandidatesCA = RecoTracker.CkfPattern.CkfTrackCandidates_cfi.ckfTrackCandidates.clone(
     src = cms.InputTag('initialStepSeedsCA'),
     ### these two parameters are relevant only for the CachingSeedCleanerBySharedInput
     #numHitsForSeedCleaner = cms.int32(50),
     #onlyPixelHitsForSeedCleaner = cms.bool(True),
 	 #onlyPixelHitsForSeedCleaner = cms.bool(False),
     TrajectoryBuilderPSet = cms.PSet(refToPSet_ = cms.string('initialStepTrajectoryBuilderCA')),
     doSeedingRegionRebuilding = True,
     useHitsSplitting = True
)
 
# fitting
import RecoTracker.TrackProducer.TrackProducer_cfi
#initialStepTracksCA = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
#     src = 'initialStepTrackCandidatesCA',
#     AlgorithmName = cms.string('iter0'),
#     Fitter = cms.string('FlexibleKFFittingSmoother')
#     )


from TrackingTools.RecoGeometry.GlobalDetLayerGeometryESProducer_cfi import *

from TrackingTools.TrackFitters.KFTrajectoryFitterESProducer_cfi import *
#fitting stuff
##### NB: What are the rightparameters? (these are from HLT_FULL_cff.py)
RungeKuttaTrackerPropagator = cms.ESProducer( "PropagatorWithMaterialESProducer",
	PropagationDirection = cms.string( "alongMomentum" ),
	ComponentName = cms.string( "RungeKuttaTrackerPropagator" ),
	 Mass = cms.double( 0.105 ),
     ptMin = cms.double( -1.0 ),
     MaxDPhi = cms.double( 1.6 ),
	useRungeKutta = cms.bool( True )
                                              )

# TRACK FITTING AND SMOOTHING OPTIONS
import TrackingTools.TrackFitters.RungeKuttaFitters_cff
CAStepFitterSmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.KFFittingSmootherWithOutliersRejectionAndRK.clone(
     ComponentName = 'CAStepFitterSmoother',
     Fitter = cms.string('CAStepRKFitter'),
     Smoother = cms.string('CAStepRKSmoother')
 )

CAStepFitterSmootherForLoopers = CAStepFitterSmoother.clone(
     ComponentName = 'CAStepFitterSmootherForLoopers',
     Fitter = cms.string('CAStepRKFitterForLoopers'),
     Smoother = cms.string('CAStepRKSmootherForLoopers')
 )
 
  # Also necessary to specify minimum number of hits after final track fit
CAStepRKTrajectoryFitter = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectoryFitter.clone(
    ComponentName = cms.string('CAStepRKFitter'),
    minHits = 8
 )
 
CAStepRKTrajectoryFitterForLoopers = CAStepRKTrajectoryFitter.clone(
     ComponentName = cms.string('CAStepRKFitterForLoopers'),
     Propagator = cms.string('PropagatorWithMaterialForLoopers')
    )

CAStepRKTrajectorySmoother = TrackingTools.TrackFitters.RungeKuttaFitters_cff.RKTrajectorySmoother.clone(
    ComponentName = cms.string('CAStepRKSmoother'),
    errorRescaling = 10.0,
	minHits = 8    
 )
 
CAStepRKTrajectorySmootherForLoopers = CAStepRKTrajectorySmoother.clone(
    ComponentName = cms.string('CAStepRKSmootherForLoopers'),
   # Propagator = cms.string('PropagatorWithMaterialForLoopers'),
    )

 
 
import TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi
CAFlexibleKFFittingSmoother = TrackingTools.TrackFitters.FlexibleKFFittingSmoother_cfi.FlexibleKFFittingSmoother.clone(
	ComponentName = cms.string('CAFlexibleKFFittingSmoother'),
 	standardFitter = cms.string('CAStepFitterSmoother'),
	looperFitter = cms.string('CAStepFitterSmootherForLoopers'),
    )

#finally fitting tracks
initialStepTracksCA = RecoTracker.TrackProducer.TrackProducer_cfi.TrackProducer.clone(
     src = 'initialStepTrackCandidatesCA',
     AlgorithmName = cms.string('iter0'),
     Fitter = cms.string('CAFlexibleKFFittingSmoother'),
     Propagator = cms.string('RungeKuttaTrackerPropagator')
     #Propagator = cms.string('PropagatorWithMaterial')
    )
   
 
 # Final selection
import RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi
from RecoTracker.IterativeTracking.DetachedTripletStep_cff import detachedTripletStepSelector
initialStepSelectorCA = RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.multiTrackSelector.clone(
      src='initialStepTracksCA',
      useAnyMVA = cms.bool(True),
      GBRForestLabel = cms.string('MVASelectorIter0_13TeV_v0'),
      trackSelectors= cms.VPSet(
      RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.looseMTS.clone(
          name = 'initialStepLoose',
         ), #end of pset
      RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.tightMTS.clone(
          name = 'initialStepTight',
          preFilterName = 'initialStepLoose',
         ),
      RecoTracker.FinalTrackSelectors.multiTrackSelector_cfi.highpurityMTS.clone(
          name = 'initialStepV1',
          preFilterName = 'initialStepTight',
         ),
      detachedTripletStepSelector.trackSelectors[4].clone(
          name = 'initialStepV2',
          preFilterName=cms.string(''),
          keepAllTracks = cms.bool(False)
         ),
      detachedTripletStepSelector.trackSelectors[5].clone(
         name = 'initialStepV3',
         preFilterName=cms.string(''),
         keepAllTracks = cms.bool(False)
         )
     ) #end of vpset
 )#end of clone
 
 #I don't know how to handle it!
initialStepCA = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = cms.VInputTag(cms.InputTag('initialStepTracksCA'),
                                    cms.InputTag('initialStepTracksCA'),
                                    cms.InputTag('initialStepTracksCA')),
     hasSelector=cms.vint32(1,1,1),
     shareFrac = cms.double(0.99),
     indivShareFrac=cms.vdouble(1.0,1.0,1.0),
     selectedTrackQuals = cms.VInputTag(cms.InputTag("initialStepSelectorCA","initialStepV1"),
                                        cms.InputTag("initialStepSelectorCA","initialStepV2"),
                                        cms.InputTag("initialStepSelectorCA","initialStepV3")),
     setsToMerge = cms.VPSet(cms.PSet( tLists=cms.vint32(0,1,2), pQual=cms.bool(True) )),
     writeOnlyTrkQuals=cms.bool(True)
 )
 
 # Final sequence
InitialStepCAPentuplets = cms.Sequence(PentupletLayers*
                                       initialStepSeedsCA*
                                       #MeasurementTrackerEvent*
                                       initialStepTrackCandidatesCA*
                                       initialStepTracksCA*
                                       initialStepSelectorCA*
                                       initialStepCA)
