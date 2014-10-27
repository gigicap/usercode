// -*- C++ -*-
// Package:    CAHitsGenerator
// Class:      CAHitsGenerator
// 
/**\class CAHitsGenerator CAHitsGenerator.cc RecoTracker/CAtracker/plugins/CAHitsGenerator.cc*/
//
//  Author:  Gigi Cappello
//  Created:  Wed, 22 Lug 2014 
//
/*This is a debug version of the code, to be used with CAtest_cfg
Two different flags (to be set in the cfg) can be used:
-  debug controls the verbose level of the output 
	* 0 = no output
	* 1 = minimal information (summary of the production)
	* 2 = summary and and some milestones
	* >=3 = very verbose outup (cell per cell info). Use carefully!
- maketrees triggers the production of a root file with useful infos about the cells and triplet production
	* false = no file
	* true = a file with default name histotrip.root is created     
For a standard usage maintain both the flags set to the default value (=0)  
*/

// system include files
#include "CAHitsGenerator.h"  //class CAcell defined here
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <chrono>
#include <math.h>


using namespace ctfseeding;
    

CAHitsGenerator::CAHitsGenerator(const edm::ParameterSet& cfg , edm::ConsumesCollector& iC) :
theSeedingLayerToken(iC.consumes<SeedingLayerSetsHits>(cfg.getParameter<edm::InputTag>("SeedingLayers")))
    {

     //get parameters from cfg
     m_debug = cfg.getUntrackedParameter<int>("debug");
     m_tree = cfg.getUntrackedParameter<bool>("maketrees");
     m_builderName = cfg.getUntrackedParameter <std::string> ("Builder");
     dEta_cut = cfg.getParameter<double>("EtaCut");
        
    //triplet producer definition
    oPset = cfg.getParameter<edm::ParameterSet>("GeneratorPSet");
    oName = oPset.getParameter<std::string>("ComponentName");

    //theGenerator.reset(MultiHitGeneratorFromPairAndLayersFactory::get()->create(oName, oPset));
    theGenerator.reset(HitTripletGeneratorFromPairAndLayersFactory::get()->create(oName, oPset, iC));
    theGenerator->init(HitPairGeneratorFromLayerPair( 0, 1, &theLayerCache), &theLayerCache);
        
	//Store tree generation
	if(m_tree){
	edm::Service < TFileService > fs;
	triptree = fs->make<TTree> ("triptree","");
	evtree = fs->make<TTree> ("evtree","");
	seedhist = fs->make<TH1D> ("seedhist","#eta of the produced seeds",100,-4.0,4.0);
    fittreegenerator();
	}
	
    eventcounter = 0;
    max_status = 3;
	
    initialised = false;
}


CAHitsGenerator::~CAHitsGenerator()
{
   // do anything here that needs to be done at desctruction time
   // Anything more to deallocate? 
   //resultCA.clear();
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void CAHitsGenerator::init(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

    if (m_debug >1)
    	std::cout<<"Begin init---"<<std::endl;
    
    
    iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
    iSetup.get<TransientRecHitRecord>().get(m_builderName, builder);
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);
    //many differences with the 7_0_0 version of the code
    
initialised = true;
    
return;
}


void CAHitsGenerator::hitSets(const TrackingRegion& region, OrderedMultiHits & resultCA, const edm::Event& iEvent, const edm::EventSetup& iSetup)
//void CAHitsGenerator::hitSets(OrderedMultiHits& resultCA, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
        
   
      //for time measurements (developer version)
     
   auto t_total = std::chrono::system_clock::now();

	if(!initialised) init(iEvent,iSetup);
	
	//BASE ONE 
	OrderedHitTriplets *scoll = new OrderedHitTriplets();
	
	//resultCA.reserve(10000);
	scoll->reserve(10000);
    //scoll->reset();
    //scoll = std::unique_ptr<OrderedHitTriplets>(new OrderedHitTriplets());
    //for time measurements (developer version)
    auto t_triplets = std::chrono::system_clock::now();
    
    //std::auto_ptr<OrderedHitTriplets> scoll(new OrderedHitTriplets);

     
    if (m_debug > 1)  std::cout<<"Generating triplets"<<std::endl;
    
    edm::Handle<SeedingLayerSetsHits> hlayers;
    iEvent.getByToken(theSeedingLayerToken, hlayers);
        
    const SeedingLayerSetsHits& layers = *hlayers;

    if(layers.numberOfLayersInSet() != 3)
        throw cms::Exception("Configuration") << "CombinedMultiHitGenerator expects SeedingLayerSetsHits::numberOfLayersInSet() to be 3, got " << layers.numberOfLayersInSet();
        
            
    std::vector<LayerTriplets::LayerSetAndLayers> trilayers = LayerTriplets::layers(layers);
    int trilaycount = 0;

    for(const auto& setAndLayers: trilayers) {
    	    // std::cout<<"trilaycont = "<<trilaycount<<std::endl;
    	     
        theGenerator->setSeedingLayers(setAndLayers.first, setAndLayers.second);
        theGenerator->hitTriplets(region, *scoll, iEvent, iSetup);
        
        if (m_debug > 1)  std::cout<<"scoll size : "<<scoll->size()<<std::endl;
        trilaycount++;
        }

     
    //produce triplets

   theLayerCache.clear();


    //the definition of scoll differs from the previous version of the code
    
    tripletCollection.reserve(scoll->size());
    fittedTripletCollection.reserve(scoll->size());
    resultCA.reserve(scoll->size());

    
    auto t_triplets_end = std::chrono::system_clock::now();
    auto t_tripletsE = std::chrono::duration_cast<std::chrono::microseconds>(t_triplets_end - t_triplets);
    auto t_casetup = std::chrono::system_clock::now();

    

    for (size_t is = 0; is<scoll->size(); is++) {
        CAcellGenerator(builder, (*scoll)[is], is);
            if (m_debug > 1)  std::cout << "Cells have been generated " << std::endl;
        
        hitUsage[tripletCollection[is].rawId0].push_back(&tripletCollection[is]);
        hitUsage[tripletCollection[is].rawId1].push_back(&tripletCollection[is]);
        hitUsage[tripletCollection[is].rawId2].push_back(&tripletCollection[is]);
        
    }
    
  
// std::cout<<"NUMBER OF TRIPLETS: "<<tripletCollection.size()<<std::endl;
    int  prodCells123 = 0;
    int  prodCells123eta = 0;
    int  prodCells234 = 0;
    int  prodCells345 = 0;
    int prodCells10343 = 0;
	int  prodCells10121 = 0;
	int  prodCells11312 = 0;
	int prodCells11531 = 0;

	invalidseqcounter = 0;

    
//Loop over the triplets to build the neighborhood map
//i.e. fill the "HasNeighbor" var and IsNeighborly
int ncont = 0;
for (size_t it = 0; it<tripletCollection.size(); it++) {
    int IsNeighbor = 0;
    
    
        std::list<CAcell *> list_h0 = hitUsage[tripletCollection[it].rawId0];
        std::list<CAcell *> list_h1 = hitUsage[tripletCollection[it].rawId1];
        std::list<CAcell *> list_h2 = hitUsage[tripletCollection[it].rawId2];
    
    
        std::vector<CAcell *> j_list01 = ListIntersect(list_h0, list_h1);
        std::vector<CAcell *> j_list12 = ListIntersect(list_h1, list_h2);
    
    
    if(m_debug>2){
    std::cout<<"list_h0 size == "<<list_h0.size()<<std::endl;
    std::cout<<"list_h1 size == "<<list_h1.size()<<std::endl;
    std::cout<<"list_h2 size == "<<list_h2.size()<<std::endl;
	std::cout<<"----------------------------------"<<std::endl;
   std::cout<<"j_list01 size == "<<j_list01.size()<<std::endl;
   std::cout<<"j_list12 size == "<<j_list12.size()<<std::endl;
    
    }


     for (size_t il = 0; il<j_list01.size(); il++){
         if((j_list01[il] != &tripletCollection[it]) && (j_list01[il]->tripletIdentifier != tripletCollection[it].tripletIdentifier) && fabs(j_list01[il]->dEta - tripletCollection[it].dEta)< dEta_cut)
             tripletCollection[it].left_neighbors.push_back(j_list01[il]);
     }
    
    
    for (size_t il = 0; il<j_list12.size(); il++){
        if(j_list12[il] != &tripletCollection[it] && (j_list12[il]->tripletIdentifier != tripletCollection[it].tripletIdentifier) && fabs(j_list12[il]->dEta - tripletCollection[it].dEta)< dEta_cut)
            tripletCollection[it].right_neighbors.push_back(j_list12[il]);
    }
    
    if (m_debug > 2 ){ 
    std::cout<<"left_list size == "<<tripletCollection[it].left_neighbors.size()<<std::endl;
    std::cout<<"right_list size == "<<tripletCollection[it].right_neighbors.size()<<std::endl;
    }
 //   std::cout<<"done "<<std::endl;

    
    IsNeighbor = tripletCollection[it].left_neighbors.size()+tripletCollection[it].right_neighbors.size();
    
    tripletCollection[it].IsNeighborly = IsNeighbor;
    
    if(m_debug> 2)
        std::cout<<"triplet has = "<<tripletCollection[it].IsNeighborly<<" neighbors. Is in lays: "<<tripletCollection[it].tripletIdentifier<<std::endl;
    
    
    list_h0.clear();
    list_h1.clear();
    list_h2.clear();
    j_list01.clear();
    j_list12.clear();
    
	//End of the neighborhood map definition. Info stored in the cells
    
    //the following variables are used for debug only, up to now
    if(tripletCollection[it].tripletIdentifier == 321){
        prodCells123++;
	if(fabs(tripletCollection[it].dEta) < 1.3)
		         prodCells123eta++;
		}
	if(tripletCollection[it].tripletIdentifier == 432) prodCells234++;
	if(tripletCollection[it].tripletIdentifier == 543) prodCells345++;
	//Add the same info for triplets in the EndCaps 
	if(tripletCollection[it].tripletIdentifier == 10343) prodCells10343++;
	if(tripletCollection[it].tripletIdentifier == 10121) prodCells10121++;
	if(tripletCollection[it].tripletIdentifier == 11312) prodCells11312++;
	if(tripletCollection[it].tripletIdentifier == 11531) prodCells11531++;


    //fit only the triplets with neighbors
    if (tripletCollection[it].IsNeighborly != 0) {
        ncont++;
    // fast helix fitter
        fitTripletSeeds(&tripletCollection[it], theMF);
 	//store the neighborly and fitted triplets
    	if (tripletCollection[it].FitSuccessful) {
            tripletCollection[it].CAstatus = 1;     //initialize CAstatus
            fittedTripletCollection.push_back(&tripletCollection[it]);
        }
    }
        

}


//for time measurements (developer version)
    auto t_casetup_end 	= std::chrono::system_clock::now();
    auto t_casetupE 	= std::chrono::duration_cast<std::chrono::microseconds>(t_casetup_end - t_casetup);
    auto t_caforward 	= std::chrono::system_clock::now();

//Loops on fittedtriplets -> Forward CA
    int n_fiter = -1;
    n_fiter = ForwardCA();
    
    //for time measurements (developer version)
//if (m_debug > 2){
    auto t_caforward_end = std::chrono::system_clock::now();
    auto t_caforwardE = std::chrono::duration_cast<std::chrono::microseconds>(t_caforward_end - t_caforward);
    auto t_cabackward = std::chrono::system_clock::now();
//}
//Connect triplets into multiplets Backward CA
    BackwardCA();
    
//for time measurements (developer version)
//if (m_debug > 2){
    auto t_cabackward_end = std::chrono::system_clock::now();
    auto t_cabackwardE = std::chrono::duration_cast<std::chrono::microseconds>(t_cabackward_end - t_cabackward);
    auto t_multprod = std::chrono::system_clock::now();
//}

int nSeeds = 0;

for(size_t im = 0; im < multiplets.size(); im++){
	//std::vector<TransientTrackingRecHit::ConstRecHitPointer> multiHitPointer;
    std::vector<ConstRecHitPointer>  multiHitPointer;
	for(size_t ii = 0; ii < multiplets[im].size(); ii++)
		multiHitPointer.push_back(multiplets[im][ii]);
	if(multiHitPointer.size()!=5)
		std::cout<<"WARNING: multiset size = "<<multiHitPointer.size()<<std::endl;
	else{
    		SeedingHitSet multiset(multiHitPointer[4], multiHitPointer[3], multiHitPointer[2], multiHitPointer[1], multiHitPointer[0]);
			resultCA.push_back(multiset);
			nSeeds++;
			GlobalPoint p_last = multiHitPointer[0]->globalPosition();
			double pl_eta = p_last.eta();
			penteta.push_back(pl_eta);
	}
}

//fill the penteta tree
if(m_tree){ 
	for (size_t itree = 0; itree<penteta.size(); itree++) {
       seedhist->Fill(penteta[itree]);
    }

 }	
//for time measurements (developer version)
//if (m_debug > 2){
   auto t_multprod_end = std::chrono::system_clock::now();
   auto t_multprodE = std::chrono::duration_cast<std::chrono::microseconds>(t_multprod_end - t_multprod);
   auto t_total_end = std::chrono::system_clock::now();
   auto t_totalE = std::chrono::duration_cast<std::chrono::microseconds>(t_total_end - t_total);
//}

//fill the debug tree 
if(m_tree){    //std::cout<<"Generated!"<<std::endl;
    for (size_t itree = 0; itree<fittedTripletCollection.size(); itree++) {
        fillfittree(fittedTripletCollection[itree]);
    }
//if(m_tree){    //TO CHANGE BACK
//    for (size_t itree = 0; itree<tripletCollection.size(); itree++) {
//        fillfittree(&tripletCollection[itree]);
//    }

    
//fill event tree
tmp_procTriplets = scoll->size();  
tmp_prodCells = tripletCollection.size();
tmp_prodCells123 = prodCells123;
tmp_prodCells123eta = prodCells123eta;
tmp_prodCells234 = prodCells234;
tmp_prodCells345 = prodCells345;
tmp_prodCells10343 = prodCells10343;
tmp_prodCells10121 = prodCells10121;
tmp_prodCells11312 = prodCells11312;
tmp_prodCells11531 = prodCells11531;

tmp_neigCells = ncont; 
tmp_etaCValue = dEta_cut;  
tmp_fitCells = fittedTripletCollection.size();  
tmp_nSteps = n_fiter;
tmp_nSeeds = nSeeds; 
  
//fill timing tree
if(m_debug > 2){
	tmp_t_total = t_totalE.count();
	tmp_t_triplets = t_tripletsE.count();
	tmp_t_casetup = t_casetupE.count();
	tmp_t_caforward = t_caforwardE.count();
	tmp_t_cabackward = t_cabackwardE.count();
	tmp_t_multprod = t_multprodE.count();
	}

evtree->Fill();
}

std::cout<<"Final number of multiSeeds = "<<resultCA.size()<<std::endl;

if(m_debug > 0){
std::cout<<"================================"<<std::endl;
std::cout<<"Processed "<<scoll->size()<<" seeds"<<std::endl;
std::cout<<"produced  "<<tripletCollection.size()<<" CAcells, of whom:"<<std::endl;
std::cout<<ncont<<" have neighbors"<<std::endl;
std::cout<<fittedTripletCollection.size()<<" have been successfully fitted!"<<std::endl;
std::cout<<"Forward iterations took n_steps = "<<n_fiter<<std::endl;
std::cout<<"Final number of multiSeeds = "<<resultCA.size()<<std::endl;
std::cout<<"================================"<<std::endl;
}

//a brand new debug
//std::cout<<tripletCollection.size()<<"\t"<<fittedTripletCollection.size()<<"\t"<<n_fiter<<"\t"<<resultCA.size()<<std::endl;

/*********************************************/        
//Clear all  /*********************************************/        
//Clear all    

	
    //try with swap to actually un-allocate memory 
    delete scoll;
    tripletCollection.clear();
	std::vector<CAcell>().swap(tripletCollection);
    fittedTripletCollection.clear();
	std::vector<CAcell*>().swap(fittedTripletCollection);
    hitUsage.clear();

    multiplets.clear();
	//std::vector<TransientTrackingRecHit::RecHitContainer>().swap(multiplets);
        resultCA.clear();
	penteta.clear();

eventcounter++;
//std::cout<<"last row: event counter #"<<eventcounter<<std::endl;
}

void CAHitsGenerator::CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder> & builder, SeedingHitSet sset,  int seednumber){
    //Build the CA cells, i.e. triplets (objects containing RecHits and fit parameters
        if (m_debug>2) std::cout << " building a ca cell "<< std::endl;

   	CAcell trip;
    //translate hits into a SeedingHitSet
    
    size_t i = 0;
    
    std::vector<const BaseTrackerRecHit*> starting_triplet;
    
    while(i < sset.size()){
    	if(sset[i]!=nullptr)  starting_triplet.push_back(sset.get(i));
    	i++;
    }
      
   if (m_debug>2) std::cout << " starting triplet size:  "<<starting_triplet.size()<< std::endl;

    
	trip.hits = starting_triplet;
    
    //TrackinkgRecHits for fast Id access
    
    if (m_debug>2) std::cout << " setting hi0/1/2 "<< std::endl;
    
    //const TrackingRecHit* hi0 = sset[0]->hit();
    //const TrackingRecHit* hi1 = sset[1]->hit();
    //const TrackingRecHit* hi2 = sset[2]->hit();
    
    //const BaseTrackerRecHit* hi0 = sset.get(0);
    //const BaseTrackerRecHit* hi1 = sset.get(1);
    //const BaseTrackerRecHit* hi2 = sset.get(2);
    
    auto hi0 = starting_triplet[0];
    auto hi1 = starting_triplet[1];
    auto hi2 = starting_triplet[2];


/*    unsigned int sdidin =hi0->geographicalId().subdetId();
        std::cout<<"subdetids = "<<sdidin<<std::endl;

    unsigned int sdidmid =hi1->geographicalId().subdetId();
        std::cout<<"subdetids = "<<sdidmid<<std::endl;

    unsigned int sdidout =hi2->geographicalId().subdetId();
    std::cout<<"subdetids = "<<sdidout<<std::endl;*/
    
      if (m_debug>2) std::cout << " setting rawId's "<< std::endl;

      trip.rawId0 = OmniRef(hi0);
      trip.rawId1 = OmniRef(hi1);
      trip.rawId2 = OmniRef(hi2);

      if (m_debug>2) std::cout << " setting global points / 1"<< std::endl;


    GlobalPoint p0 = hi0->globalPosition();
    GlobalPoint p1 = hi1->globalPosition();
    GlobalPoint p2 = hi2->globalPosition();
    
          if (m_debug>2) std::cout << " setting global points / 2"<< std::endl;


    trip.p0 = p0;	
    trip.p1 = p1;	   
    trip.p2 = p2;

    trip.dEta = DeltaEta(p0,p2);
    
          if (m_debug>2) std::cout << " setting detId's "<< std::endl;

    
    DetIDClassification hit0cl(hi0->geographicalId());
    DetIDClassification hit1cl(hi1->geographicalId());
    DetIDClassification hit2cl(hi2->geographicalId());
    trip.hitId0 = hit0cl.getLayer();
    trip.hitId1 = hit1cl.getLayer();
    trip.hitId2 = hit2cl.getLayer();


    trip.eventNumber = eventcounter;
    trip.tripletIdentifier = trip.hitId0+10*trip.hitId1+100*trip.hitId2;

    
    if (m_debug>2) {
        std::cout << " layer 0 =  " << trip.hitId0 << std::endl;
        std::cout << " layer 1 =  " << trip.hitId1 << std::endl;
        std::cout << " layer 2 =  " << trip.hitId2 << std::endl;
        std::cout << " identifier = "<< trip.tripletIdentifier <<std::endl;
    }
    
    
    trip.TripletNumber = seednumber;
	trip.FitSuccessful = false;
	trip.LocalMomentum = 0;

	//Is Used (for the final Backwards CA)
	trip.IsUsed = 0;


	//check unrecognized triplets (with hits from layer sequencies different from those admitted. 
	//This should include also problems related to two-sided layers)
	//if(IsValidSeq(trip.tripletIdentifier))
    				tripletCollection.push_back(trip);
	//else{
	//	std::cout<<"Invalid triplet sequence: "<<trip.tripletIdentifier<<std::endl;
	//	invalidseqcounter++;
	//}

	return; 
}

void CAHitsGenerator::fitTripletSeeds(CAcell *trip, edm::ESHandle<MagneticField> & mf) {

 
    
    double nomField = mf->nominalValue();

    FastHelix fit(trip->p2, trip->p1, trip->p0, nomField , mf.product());
    
    GlobalTrajectoryParameters params = fit.stateAtVertex();
    
    if (fit.isValid()){
    trip->FitSuccessful = true;
	trip->LocalMomentum = params.momentum().mag();
    }

return;
}



int CAHitsGenerator::ForwardCA(){
    int Stop = -1;
    int step_iterator = 0;
    
    int it_max = 20; //Maximum number of iteraitons
        
    while (Stop < (int)fittedTripletCollection.size() && step_iterator<it_max) {
        
        step_iterator++;
        Stop = 0;
        
        
        for (size_t t =0; t<fittedTripletCollection.size(); t++) {
            //std::cout<<" ________ "<<std::endl;
            //std::cout<<"triplet no. "<<t<<std::endl;
            if (fittedTripletCollection[t]->left_neighbors.size()!=0) {
                for (size_t il = 0; il<fittedTripletCollection[t]->left_neighbors.size(); il++) {
                    bool en_cont = false;   //useful variable
			if(m_debug > 2){
                    	std::cout<<"neighbor list no = "<<il<<std::endl;
                    	std::cout<<"fittedTriplet status = "<<fittedTripletCollection[t]->CAstatus<<"left_neighbors status = "<<fittedTripletCollection[t]->left_neighbors[il]->CAstatus <<std::endl;
			}
                    if (fittedTripletCollection[t]->CAstatus <= fittedTripletCollection[t]->left_neighbors[il]->CAstatus && en_cont == false){
                        (fittedTripletCollection[t]->CAstatus)++;
                        en_cont = true;
                    }
                if (en_cont == false)
                    Stop++;
                
                }
              
            }
            else
                Stop++;
        }
        
        //std::cout<<"Stop value = "<<Stop<<std::endl;
    }
    
    //debug
if(m_debug>2){
    for (size_t t =0; t<fittedTripletCollection.size(); t++) {
        std::cout<<" ._._._._._._._._ "<<std::endl;
        std::cout<<" cell number "<<t<<std::endl;
        std::cout<<" processed with status: "<<fittedTripletCollection[t]->CAstatus<<std::endl;
        
    }
    }
    
    if(step_iterator>=it_max)
        std::cout<<"WARNING: too many CA iterations!"<<std::endl;
    
    return step_iterator;
}


void CAHitsGenerator::BackwardCA(){

int multi_id = 1;

for (size_t t =0; t<fittedTripletCollection.size(); t++) {
//start only with non-used triplets with maximum CAstatus
if(fittedTripletCollection[t]->CAstatus >= max_status && fittedTripletCollection[t]->IsUsed==0){
	CAcell *current_cell = fittedTripletCollection[t];
	//TransientTrackingRecHit::RecHitContainer multicontainer;
	std::vector<ConstRecHitPointer> multicontainer;
	while(current_cell){
		current_cell = define_used(current_cell, multi_id, &multicontainer);
		}
		multi_id++;	
	multiplets.push_back(multicontainer);
	multicontainer.clear();
	}
}

return;
}

CAcell* CAHitsGenerator::define_used(CAcell *cell, int id, std::vector<ConstRecHitPointer>*multicontainer){
cell->IsUsed = id;
double pt_tmp = 99999.0;

if(cell->CAstatus==1){
	multicontainer->push_back(cell->hits[2]);
	multicontainer->push_back(cell->hits[1]);
	multicontainer->push_back(cell->hits[0]);
return NULL;	
}

else{
CAcell *next = new CAcell();
multicontainer->push_back(cell->hits[2]);
	for(size_t il = 0; il< cell->left_neighbors.size(); il++){
		if(cell->CAstatus - cell->left_neighbors[il]->CAstatus == 1){
		double pt_diff = fabs(cell->getLocalMomentum() - cell->left_neighbors[il]->getLocalMomentum());
		if(pt_diff <= pt_tmp){
			pt_tmp = pt_diff;
			next = cell->left_neighbors[il];
			} 
		}
	}
return next;
}

}



//Useful tools to store triplets in root files
void CAHitsGenerator::fittreegenerator(){

triptree->Branch("Event",&tmp_event,"Event/I");
//triptree->Branch("FitChiSquared",&tmp_FitChiSquared,"FitChiSquared/D");
//triptree->Branch("FitNDof",&tmp_FitNdof,"FitNdof/I");
//triptree->Branch("FitFoundHits",&tmp_FitFoundHits,"FitFoundHits/I");
triptree->Branch("LocalMomentum",&tmp_LocalMomentum,"LocalMomentum/D");

//triptree->Branch("hqbp",tmp_hqbp,"hqbp[3]/D");
//triptree->Branch("hdxdz",tmp_hdxdz,"hdxdz[3]/D");
//triptree->Branch("hdydz",tmp_hdydz,"hdydz[3]/D");

triptree->Branch("hitId0",&tmp_hitid0,"hitId0/I");
triptree->Branch("hitId1",&tmp_hitid1,"hitId1/I");
triptree->Branch("hitId2",&tmp_hitid2,"hitId2/I");
    
triptree->Branch("rawId0",&tmp_rawid0,"rawId0/I");
triptree->Branch("rawId1",&tmp_rawid1,"rawId1/I");
triptree->Branch("rawId2",&tmp_rawid2,"rawId2/I");
    
triptree->Branch("Eta",&tmp_dEta,"Eta/D");
    
    
triptree->Branch("Neighborly",&tmp_neighborly,"Neighborly/I");
triptree->Branch("CAstatus",&tmp_CAstatus,"CAstatus/I");
    
    
triptree->Branch("Identifier",&tmp_identifier,"Identifier/I");

triptree->Branch("IsUsed",&tmp_isused,"IsUsed/I");


//event tree
evtree->Branch("procTriplets",&tmp_procTriplets,"procTriplets/I");
evtree->Branch("prodCells",&tmp_prodCells,"prodCells/I");
evtree->Branch("prodCells123",&tmp_prodCells123,"prodCells123/I");
evtree->Branch("prodCells123eta",&tmp_prodCells123eta,"prodCells123eta/I");
evtree->Branch("prodCells234",&tmp_prodCells234,"prodCells234/I");
evtree->Branch("prodCells345",&tmp_prodCells345,"prodCells345/I");
evtree->Branch("prodCells10343",&tmp_prodCells10343,"prodCells10343/I");
evtree->Branch("prodCells10121",&tmp_prodCells10121,"prodCells10121/I");
evtree->Branch("prodCells11312",&tmp_prodCells11312,"prodCells11312/I");
evtree->Branch("prodCells11531",&tmp_prodCells11531,"prodCells11531/I");
evtree->Branch("neigCells",&tmp_neigCells,"neigCells/I");
evtree->Branch("etaCValue",&tmp_etaCValue,"etaCValue/D");
evtree->Branch("fitCells",&tmp_fitCells,"fitCells/I");
evtree->Branch("nSteps",&tmp_nSteps,"nSteps/I");
evtree->Branch("nSeeds",&tmp_nSeeds,"nSeeds/I");
evtree->Branch("nInvalid",&invalidseqcounter,"nInvalid/I");

    
//timetree
    evtree->Branch("t_total",&tmp_t_total , "t_total/I");
    evtree->Branch("t_triplets",&tmp_t_triplets , "t_triplets/I");
    evtree->Branch("t_casetup",&tmp_t_casetup, "t_casetup/I");
    evtree->Branch("t_caforward",&tmp_t_caforward  , "t_caforward/I");
    evtree->Branch("t_cabackward",&tmp_t_cabackward , "t_cabackward/I");
    evtree->Branch("t_multiprod",&tmp_t_multprod  , "t_multprod/I");

    
return;
}


void CAHitsGenerator::fillfittree(CAcell *t){

tmp_LocalMomentum = t->getLocalMomentum();

tmp_event = t->geteventNumber();

tmp_hitid0 = t->gethitId0();
tmp_hitid1 = t->gethitId1();
tmp_hitid2 = t->gethitId2();
    
tmp_rawid0 = t->rawId0;
tmp_rawid1 = t->rawId1;
tmp_rawid2 = t->rawId2;
    
tmp_dEta = t->dEta;
    
tmp_neighborly = t->IsNeighborly;
tmp_identifier = t->tripletIdentifier;

tmp_CAstatus = t->CAstatus;

tmp_isused = t->IsUsed;
    
triptree->Fill();
    

return;
}


//Other functions

std::vector<CAcell *> CAHitsGenerator::ListIntersect(std::list<CAcell *> list_h1, std::list<CAcell *> list_h2){
    
    std::vector<CAcell *> temp;
    
    std::list<CAcell *>::const_iterator L1itr = list_h1.begin();
    std::list<CAcell *>::const_iterator L2itr = list_h2.begin();
    
    
    while(L1itr != list_h1.end() && L2itr != list_h2.end()){
        
        //check to see if they're equal, add to our temp list
        if(*L1itr == *L2itr){
            temp.push_back(*L1itr);
            L1itr++;
            L2itr++;
        }
        else if(*L1itr < *L2itr){
            L1itr++;
        }
        else{
            L2itr++;
        }
    }
    return temp;
}



unsigned int CAHitsGenerator::OmniRef(const BaseTrackerRecHit* rhit){

    unsigned int ocref = 0;
    
    int subdetid = rhit->geographicalId().subdetId();
    if (subdetid==PixelSubdetector::PixelBarrel||subdetid==PixelSubdetector::PixelEndcap) {
        const SiPixelRecHit* pRHit = dynamic_cast<const SiPixelRecHit*>(rhit);
        if (!pRHit->cluster().isNonnull())
            std::cout<< ">>> RecHit does not have an associated cluster!" << std::endl;
        ocref = pRHit->omniClusterRef().rawIndex();
    }
    
    else if (subdetid==SiStripDetId::TIB||subdetid==SiStripDetId::TOB||subdetid==SiStripDetId::TID||subdetid==SiStripDetId::TEC) {
        const std::type_info &tid = typeid(*rhit);
        
        if (tid == typeid(SiStripMatchedRecHit2D)) {
            const SiStripMatchedRecHit2D* sMatchedRHit = dynamic_cast<const SiStripMatchedRecHit2D*>(rhit);
            if (!sMatchedRHit->monoHit().cluster().isNonnull() || !sMatchedRHit->stereoHit().cluster().isNonnull())
								std::cout<< ">>> RecHit does not have an associated cluster!" << std::endl;
            ocref = sMatchedRHit->monoClusterRef().rawIndex() + sMatchedRHit->stereoClusterRef().rawIndex();
        }
        else if (tid == typeid(SiStripRecHit2D)) {
            const SiStripRecHit2D* sRHit = dynamic_cast<const SiStripRecHit2D*>(rhit);
            if (!sRHit->cluster().isNonnull())
                std::cout<< ">>> RecHit does not have an associated cluster!" << std::endl;
            ocref = sRHit->omniClusterRef().rawIndex();
        }
        else if (tid == typeid(SiStripRecHit1D)) {
            const SiStripRecHit1D* sRHit = dynamic_cast<const SiStripRecHit1D*>(rhit);
            if (!sRHit->cluster().isNonnull())
                std::cout<< ">>> RecHit does not have an associated cluster!" << std::endl;
            ocref = sRHit->omniClusterRef().rawIndex();
        }
        
        else {
            std::cout << ">>> getMatchedClusters: TrackingRecHit not associated to any SiStripCluster! subdetid = " << subdetid <<std::endl;
        }
        
    }
    else {
        std::cout <<  ">>> getMatchedClusters: TrackingRecHit not associated to any cluster! subdetid = " << subdetid <<std::endl;
    }
    
return ocref;
}


double CAHitsGenerator::DeltaEta(GlobalPoint p1, GlobalPoint p2){
    
    GlobalVector gV(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
    double eta = gV.eta();
    
    return eta;
}

//a dummy check on sequencies... to be removed (or get smarter) in the future
bool CAHitsGenerator::IsValidSeq(int seq){
	int valid_sequencies[5] = {543 , 432 ,321 , 10343 , -10257};
	int i = 0;
	while (i<5 && seq!=valid_sequencies[i]) i++;
	if(i==5) 
		return false;
	else 
		return true;	
}

//that's all folks!!!
