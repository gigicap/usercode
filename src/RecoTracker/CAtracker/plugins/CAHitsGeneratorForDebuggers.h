#ifndef CAHitsGeneratorForDebuggers_h
#define CAHitsGeneratorForDebuggers_h

// cmssw includes
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2DCollection.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"


#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGenerator.h"
#include "RecoTracker/TkSeedingLayers/interface/OrderedSeedingHits.h"
#include "RecoTracker/TkSeedGenerator/interface/MultiHitGenerator.h"
#include "RecoTracker/TkSeedingLayers/interface/OrderedMultiHits.h"

#include "RecoTracker/TkSeedGenerator/interface/SeedCreatorFactory.h"
#include "RecoTracker/TkSeedGenerator/interface/FastHelix.h"

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducerFactory.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGeneratorFromPairAndLayersFactory.h"
#include "RecoTracker/TkSeedGenerator/interface/MultiHitGeneratorFromPairAndLayers.h"
#include "RecoTracker/TkSeedGenerator/interface/MultiHitGeneratorFromPairAndLayersFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGeneratorFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"        
#include "RecoTracker/TkSeedGenerator/plugins/SeedFromConsecutiveHitsCreator.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/BaseTrackerRecHit.h"


#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSets.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"
#include "RecoTracker/TkHitPairs/interface/HitPairGeneratorFromLayerPair.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGeneratorFromPairAndLayers.h"
#include "RecoPixelVertexing/PixelTriplets/interface/HitTripletGeneratorFromPairAndLayersFactory.h"
#include "RecoPixelVertexing/PixelTriplets/interface/LayerTriplets.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "RecoPixelVertexing/PixelTriplets/interface/OrderedHitTriplets.h"
#include "RecoPixelVertexing/PixelTriplets/plugins/PixelTripletHLTGenerator.h"
#include "RecoTracker/TkSeedGenerator/plugins/MultiHitGeneratorFromChi2.h"
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"

#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducerFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegion.h"
#include "RecoTracker/TkTrackingRegions/interface/GlobalTrackingRegion.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Utilities/interface/EDGetToken.h"

//root
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TVector3.h"

#include <unordered_map>


#include "LayerClassification.h"
#include "CACell.h"

//
// CAHitsGeneratorForDebuggers class declaration
//


class CAHitsGeneratorForDebuggers : public MultiHitGenerator {
   public:
      using ConstRecHitPointer = BaseTrackerRecHit const *;
    
      CAHitsGeneratorForDebuggers(const edm::ParameterSet& ,  edm::ConsumesCollector& iC);
      ~CAHitsGeneratorForDebuggers();
    
      void init(const edm::Event& ev, const edm::EventSetup& es);
      void hitSets(const TrackingRegion& , OrderedMultiHits& , const edm::Event& , const edm::EventSetup& );
	  
	//Generate CA cells from the produced triplets
    void CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder>& builder, SeedingHitSet seed,  int seednumber);
	//Fit neighborly cells
	void fitTripletSeeds(CAcell *trip, edm::ESHandle<MagneticField> & mf);
	//Forward step of the CA algo
    int ForwardCA();
	//Backward step of the CA algo
	void BackwardCA(int m_status);
	//mark a cell as "used" by the CA
	CAcell* define_used(CAcell *cell, int id,std::vector<ConstRecHitPointer>* multicontainer);
	//finalforward cleaner
	void FinalForwardCleaning();
	/**Some useful functions**/
 	//Intersect the neighbor lists of two cells (used in the Neighborhood map definition)
        std::vector<CAcell *> ListIntersect(std::list<CAcell *>, std::list<CAcell *>);
	//Used to give a unique Id to every hit
        unsigned int OmniRef(const BaseTrackerRecHit*);
        //Delta eta
        double DeltaEta(GlobalPoint, GlobalPoint);
        //avoid duplicates
        bool IsValidSeq(int seq);
		int IsAtLeft(int identif);
		std::vector<int> IsAtRight(int identif);


//maximum CA_status
int max_status; //Now set to 10, in general to build n-tuplets set to > n-2


//temp storage variables for debug tree 
double tmp_FitChiSquared;
int tmp_FitNdof;
int tmp_event;
int tmp_FitFoundHits;
double tmp_LocalMomentum;
double tmp_hqbp[3];
double tmp_hdxdz[3];
double tmp_hdydz[3];
int tmp_hitid0;
int tmp_hitid1;
int tmp_hitid2;
unsigned int tmp_rawid0;
unsigned int tmp_rawid1;
unsigned int tmp_rawid2;
double tmp_dEta;
int tmp_identifier;
int tmp_neighborly;
int tmp_CAstatus;
int tmp_isused;

int tmp_procTriplets;  
int tmp_prodCells;
int tmp_prodLay1Cells;
int tmp_prodCells123;
int tmp_prodCells123eta;
int tmp_prodCells234;
int tmp_prodCells345;

int tmp_prodCells10121;
int tmp_prodCells11312;
int tmp_prodCells11531;
int tmp_prodCells10343;
int tmp_prodCells_10079;


int tmp_neigCells; 
double tmp_etaCValue;  
int tmp_fitCells;  
int tmp_nSteps;  
int tmp_nSeeds; 
int tmp_nSeeds5h; 
int tmp_nSeeds4h; 

//for time measurements 
long int    tmp_t_total;
long int    tmp_t_triplets;
long int    tmp_t_casetup;
long int    tmp_t_caforward;
long int    tmp_t_cabackward;
long int    tmp_t_cabackward4;
long int    tmp_t_secondforward;
long int    tmp_t_multprod;
int invalidseqcounter; 

std::vector<double> penteta;

   private:
    
        //Debug switch (set in the _cfi  ==0 no debud output, ==1 text output, ==2 text+debug tree, ==3 +time measurements)
        int m_debug;
        bool m_tree;
        bool m_secondbkg;
        bool do_finalfor;

 	
        //CAcell collection
        std::vector<CAcell> tripletCollection;
        std::vector<CAcell*> fittedTripletCollection;
    
       
        //TrajectorySeedCollection seedcollection;       //to be set
    
    
       // edm::ParameterSet m_regionProducerPSet;

    
        std::string m_builderName;
        std::string m_fitterName;

        edm::ESHandle<TrajectoryFitter> m_fitter;
        edm::ESHandle<TrackerGeometry> m_tracker;
    
    //triplet generator stuff
    edm::ParameterSet oPset;
    std::string oName;
    
    double dEta_cut;

	TTree *triptree;
	TTree *evtree;
	TH1D *seedhist;
	int eventcounter;

    std::vector<std::vector<ConstRecHitPointer>>multiplets;
	std::vector<int>final_multi_collection;		//store only the index

 
      	void fittreegenerator();
        void fillfittree(CAcell *t);
 
    
        //data structure for the neighborhood maps
        std::unordered_map<unsigned int, std::list<CAcell *>> hitUsage;


    std::unique_ptr<HitTripletGeneratorFromPairAndLayers> theGenerator;

	edm::ESHandle<TrackerGeometry> tracker;
	edm::ESHandle<TransientTrackingRecHitBuilder> builder;
	edm::ESHandle<MagneticField> theMF;

	std::string layerBuilderName;
	edm::ESHandle<SeedingLayerSetsBuilder> layerBuilder;
	
	bool initialised;
	bool makehistos;
	typedef LayerHitMapCache  LayerCacheType;
 	LayerCacheType            theLayerCache;
    
     edm::EDGetTokenT<SeedingLayerSetsHits> theSeedingLayerToken;
      // ----------member data ---------------------------
};


#endif
