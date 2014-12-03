import FWCore.ParameterSet.Config as cms


from Configuration.StandardSequences.Reconstruction_cff import * 
from RecoTracker.CAtracker.initialStepCAPentuplets_cff import * 

#modify Detached to take CA as input

#detachedTripletStepClustersCA = RecoTracker.IterativeTracking.DetachedTripletStep_cff.detachedTripletStepClusters.clone()
detachedTripletStepClusters.trajectories = cms.InputTag('initialStepTracksCA')
detachedTripletStepClusters.overrideTrkQuals = cms.InputTag('initialStepCA')
detachedTripletStepClusters.pixelClusters = cms.InputTag('siPixelClusters')
detachedTripletStepClusters.stripClusters = cms.InputTag('siStripClusters')

#correct vertex collection
detachedTripletStepSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")

lowPtTripletStepSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")

pixelPairStepSeeds.RegionFactoryPSet.RegionPSet.VertexCollection = cms.InputTag("firstStepPrimaryVerticesCA")
pixelPairStepSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")

mixedTripletStepSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")

pixelLessStepSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")
tobTecStepSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")

ak4CaloJetsForTrk.srcPVs = cms.InputTag('firstStepPrimaryVerticesCA')
firstStepGoodPrimaryVertices.src=cms.InputTag('firstStepPrimaryVerticesCA')
jetCoreRegionalStepSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")

muonSeededTracksInOutSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")
muonSeededTracksOutInSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")

duplicateTrackSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")

convStepSelector.vertices    = cms.InputTag("firstStepPrimaryVerticesCA")
photonConvTrajSeedFromSingleLeg.primaryVerticesTag = cms.InputTag('firstStepPrimaryVerticesCA')

pixelPairElectronSeeds.RegionFactoryPSet.RegionPSet.VertexCollection = cms.InputTag("firstStepPrimaryVerticesCA")

#early general tracks requires to be updated in order to include InitialStepCAPentuplets
#earlyGeneralTracksCA = RecoTracker.FinalTrackSelectors.earlyGeneralTracks_cfi.earlyGeneralTracks.clone()
earlyGeneralTracks.TrackProducers = (cms.InputTag('initialStepTracksCA'),
                       cms.InputTag('jetCoreRegionalStepTracks'),
                       cms.InputTag('lowPtTripletStepTracks'),
                       cms.InputTag('pixelPairStepTracks'),
                       cms.InputTag('detachedTripletStepTracks'),
                       cms.InputTag('mixedTripletStepTracks'),
                       cms.InputTag('pixelLessStepTracks'),
                       cms.InputTag('tobTecStepTracks')
                       )

earlyGeneralTracks.selectedTrackQuals = cms.VInputTag(cms.InputTag("initialStepCA"),
                     cms.InputTag("jetCoreRegionalStepSelector","jetCoreRegionalStep"),
                     cms.InputTag("lowPtTripletStepSelector","lowPtTripletStep"),
                     cms.InputTag("pixelPairStepSelector","pixelPairStep"),
                     cms.InputTag("detachedTripletStep"),
                     cms.InputTag("mixedTripletStep"),
                     cms.InputTag("pixelLessStep"),
                     cms.InputTag("tobTecStepSelector","tobTecStep")
                     )
                     
                    
#for the jetcoreregionalstep
initialStepTrackRefsForJets.src = cms.InputTag('initialStepTracksCA')

#use CAtracking as Initialstep
iterTrackingCA = cms.Sequence(
                             InitialStepCAPentuplets*			   
                             DetachedTripletStep*
                             LowPtTripletStep*
                             PixelPairStep*
                             MixedTripletStep*
                             PixelLessStep*
                             TobTecStep*
                 	         JetCoreRegionalStep*
                             earlyGeneralTracks*
                             muonSeededStep*
                             preDuplicateMergingGeneralTracks*
                             generalTracksSequence*
                             ConvStep*
                             conversionStepTracks
                             )
                             
initialStepSeedClusterMask.trajectories = cms.InputTag("initialStepSeedsCA")

newCombinedSeeds.seedCollections = cms.VInputTag(
      cms.InputTag('initialStepSeedsCA'),
      cms.InputTag('pixelPairStepSeeds'),
      cms.InputTag('mixedTripletStepSeeds'),
      cms.InputTag('pixelLessStepSeeds'),
      cms.InputTag('tripletElectronSeeds'),
      cms.InputTag('pixelPairElectronSeeds'),
      cms.InputTag('stripPairElectronSeeds')
      )
                             
                           
#to be used instead of trackingGlobalReco
CAglobaltracks = cms.Sequence(iterTrackingCA*electronSeedsSeq*doAlldEdXEstimators*trackExtrapolator)