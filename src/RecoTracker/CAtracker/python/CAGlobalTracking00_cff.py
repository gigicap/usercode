import FWCore.ParameterSet.Config as cms


from Configuration.StandardSequences.Reconstruction_cff import * 
#from RecoTracker.CAtracker.initialStepCAPentuplets_cff import * 
from RecoTracker.CAtracker.initialStepCAMultiStep_cff import *
#add an intermediate step 00
from RecoTracker.CAtracker.initialTripletStep_cff import *


#modify Detached to take CA as input

#detachedTripletStepClustersCA = RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepClusters.clone()
detachedTripletStepClusters.trajectories = cms.InputTag('initialTripletStepTracks')
detachedTripletStepClusters.overrideTrkQuals = cms.InputTag('initialTripletStep')
detachedTripletStepClusters.pixelClusters = cms.InputTag('siPixelClusters')
detachedTripletStepClusters.stripClusters = cms.InputTag('siStripClusters')


#early general tracks requires to be updated in order to include InitialStepCAPentuplets
#earlyGeneralTracksCA = RecoTracker.FinalTrackSelectors.earlyGeneralTracks_cfi.earlyGeneralTracks.clone()
earlyGeneralTracks.TrackProducers = (cms.InputTag('initialStepTracksCA'),
                       cms.InputTag('jetCoreRegionalStepTracks'),
                       cms.InputTag('initialTripletStepTracks'),
                       cms.InputTag('lowPtTripletStepTracks'),
                       cms.InputTag('pixelPairStepTracks'),
                       cms.InputTag('detachedTripletStepTracks'),
                       cms.InputTag('mixedTripletStepTracks'),
                       cms.InputTag('pixelLessStepTracks'),
                       cms.InputTag('tobTecStepTracks')
                       )
                       
earlyGeneralTracks.hasSelector=cms.vint32(1,1,1,1,1,1,1,1,1)
earlyGeneralTracks.indivShareFrac=cms.vdouble(1.0,0.19,0.19,0.16,0.19,0.13,0.11,0.11,0.09)   #what value shall I use for the intermediate step?
earlyGeneralTracks.setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3,4,5,6,7,8), pQual=cms.bool(True) ))

earlyGeneralTracks.selectedTrackQuals = cms.VInputTag(cms.InputTag("initialStepCA"),
                     cms.InputTag("jetCoreRegionalStepSelector","jetCoreRegionalStep"),
                     cms.InputTag("initialTripletStep"),
                     cms.InputTag("lowPtTripletStepSelector","lowPtTripletStep"),
                     cms.InputTag("pixelPairStepSelector","pixelPairStep"),
                     cms.InputTag("detachedTripletStep"),
                     cms.InputTag("mixedTripletStep"),
                     cms.InputTag("pixelLessStep"),
                     cms.InputTag("tobTecStepSelector","tobTecStep")
                     )
                     
                    
#for the jetcoreregionalstep
#?
#iter0TrackRefsForJets = cms.InputTag('initialStepTracksCA')

#use CAtracking as Initialstep
iterTrackingCA = cms.Sequence(
                             InitialStepCAPentuplets*
                             InitialTripletStep*			   
                             DetachedTripletStep*
                             LowPtTripletStep*
                             PixelPairStep*
                             MixedTripletStep*
                             PixelLessStep*
                             TobTecStep*
                 	     #JetCoreRegionalStep*   
                             earlyGeneralTracks*
                             muonSeededStep*
                             preDuplicateMergingGeneralTracks*
                             generalTracksSequence*
                             ConvStep*
                             conversionStepTracks
                             )
                             
# From electronseeds (still something to be fixed...
initialStepSeedClusterMask.trajectories = cms.InputTag("initialStepSeedsCA")

initialTripletStepSeedClusterMask = seedClusterRemover.clone(
    trajectories = cms.InputTag("initialTripletStepSeeds"),
    oldClusterRemovalInfo = cms.InputTag("initialStepSeedClusterMask")
)


newCombinedSeeds.seedCollections = cms.VInputTag(
      cms.InputTag('initialStepSeedsCA'),
      cms.InputTag('initialTripletStepSeeds'),
      cms.InputTag('pixelPairStepSeeds'),
      cms.InputTag('mixedTripletStepSeeds'),
      cms.InputTag('pixelLessStepSeeds'),
      cms.InputTag('tripletElectronSeeds'),
      cms.InputTag('pixelPairElectronSeeds'),
      cms.InputTag('stripPairElectronSeeds')
      )
          
electronSeedsSeq = cms.Sequence( initialStepSeedClusterMask*
								  initialTripletStepSeedClusterMask*
                                  pixelPairStepSeedClusterMask*
                                  mixedTripletStepSeedClusterMask*
                                  pixelLessStepSeedClusterMask*
                                  tripletElectronSeedLayers*
                                  tripletElectronSeeds*
                                  tripletElectronClusterMask*
                                  pixelPairElectronSeedLayers*
                                  pixelPairElectronSeeds*
                                  stripPairElectronSeedLayers*
                                  stripPairElectronSeeds*
                                  newCombinedSeeds)                   
                           
#to be used instead of trackingGlobalReco
CAglobaltracks = cms.Sequence(iterTrackingCA*electronSeedsSeq*doAlldEdXEstimators*trackExtrapolator)