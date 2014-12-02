// -*- C++ -*-
// Package:    CAHitsGeneratorForDebuggers
// Class:      CAHitsGeneratorForDebuggers
// 
/**\class CAHitsGeneratorForDebuggers CAHitsGeneratorForDebuggers.cc RecoTracker/CAtracker/plugins/CAHitsGeneratorForDebuggers.cc*/
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
#include "CAHitsGeneratorForDebuggers.h"  //class CAcell defined here
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <chrono>
#include <math.h>


using namespace ctfseeding;
    

CAHitsGeneratorForDebuggers::CAHitsGeneratorForDebuggers(const edm::ParameterSet& cfg , edm::ConsumesCollector& iC) :
theSeedingLayerToken(iC.consumes<SeedingLayerSetsHits>(cfg.getParameter<edm::InputTag>("SeedingLayers")))
    {

     //get parameters from cfg
     m_debug = cfg.getUntrackedParameter<int>("debug");
     m_tree = cfg.getUntrackedParameter<bool>("maketrees");
     m_builderName = cfg.getUntrackedParameter <std::string> ("Builder");
     dEta_cut = cfg.getParameter<double>("EtaCut");
     m_secondbkg = cfg.getParameter<bool>("MakeSecondBackward");
     do_finalfor = cfg.getParameter<bool>("MakeFinalForward");

        
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


CAHitsGeneratorForDebuggers::~CAHitsGeneratorForDebuggers()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void CAHitsGeneratorForDebuggers::init(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

    if (m_debug >1)
    	std::cout<<"Begin init---"<<std::endl;
    
    
    iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
    iSetup.get<TransientRecHitRecord>().get(m_builderName, builder);
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);
    
initialised = true;
    
return;
}


void CAHitsGeneratorForDebuggers::hitSets(const TrackingRegion& region, OrderedMultiHits & resultCA, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   using namespace edm;
        
   
      //for time measurements (developer version)
     
   auto t_total = std::chrono::system_clock::now();

	if(!initialised) init(iEvent,iSetup);
	
	//BASE ONE 
	OrderedHitTriplets *scoll = new OrderedHitTriplets();

    //to be changed
	scoll->reserve(10000);
    auto t_triplets = std::chrono::system_clock::now();
    
     
    if (m_debug > 1)  std::cout<<"Generating triplets"<<std::endl;
    
    edm::Handle<SeedingLayerSetsHits> hlayers;
    iEvent.getByToken(theSeedingLayerToken, hlayers);
        
    const SeedingLayerSetsHits& layers = *hlayers;

    if(layers.numberOfLayersInSet() != 3)
        throw cms::Exception("Configuration") << "CombinedMultiHitGenerator expects SeedingLayerSetsHits::numberOfLayersInSet() to be 3, got " << layers.numberOfLayersInSet();
        
            
    std::vector<LayerTriplets::LayerSetAndLayers> trilayers = LayerTriplets::layers(layers);
    int trilaycount = 0;

    for(const auto& setAndLayers: trilayers) {
        theGenerator->setSeedingLayers(setAndLayers.first, setAndLayers.second);
        theGenerator->hitTriplets(region, *scoll, iEvent, iSetup);
        
        if (m_debug > 1)  std::cout<<"scoll size : "<<scoll->size()<<std::endl;
        trilaycount++;
        }

     
    //produce triplets

   theLayerCache.clear();

    
    tripletCollection.reserve(scoll->size());
    fittedTripletCollection.reserve(scoll->size());
    
    auto t_triplets_end = std::chrono::system_clock::now();
    auto t_tripletsE = std::chrono::duration_cast<std::chrono::microseconds>(t_triplets_end - t_triplets);
    auto t_casetup = std::chrono::system_clock::now();

    
//Neighborhood map definition
    for (size_t is = 0; is<scoll->size(); is++) {
        CAcellGenerator(builder, (*scoll)[is], is);
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
	int prodCells_10079 = 0;

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
     	//added to avoid duplicates in layer superimpositions
     	if((j_list01[il] != &tripletCollection[it]) && (j_list01[il]->tripletIdentifier == IsAtLeft(tripletCollection[it].tripletIdentifier))){ 		
         if((j_list01[il]->tripletIdentifier != tripletCollection[it].tripletIdentifier) && fabs(j_list01[il]->dEta - tripletCollection[it].dEta)< dEta_cut)
             tripletCollection[it].left_neighbors.push_back(j_list01[il]);
     	}
     }
    
    
    for (size_t il = 0; il<j_list12.size(); il++){  
        if(j_list12[il] != &tripletCollection[it] && (j_list12[il]->tripletIdentifier != tripletCollection[it].tripletIdentifier) && fabs(j_list12[il]->dEta - tripletCollection[it].dEta)< dEta_cut)
            tripletCollection[it].right_neighbors.push_back(j_list12[il]);
    }
    
    if (m_debug > 4 ){ 
    std::cout<<"left_list size == "<<tripletCollection[it].left_neighbors.size()<<std::endl;
    std::cout<<"right_list size == "<<tripletCollection[it].right_neighbors.size()<<std::endl;
    }

    
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
	if(fabs(tripletCollection[it].dEta) < 2.0)
		         prodCells123eta++;
		}
	if(tripletCollection[it].tripletIdentifier == 432) prodCells234++;
	if(tripletCollection[it].tripletIdentifier == 543) prodCells345++;
	//Add the same info for triplets in the EndCaps 
	if(tripletCollection[it].tripletIdentifier == 10343) prodCells10343++;
	if(tripletCollection[it].tripletIdentifier == 10121) prodCells10121++;
	if(tripletCollection[it].tripletIdentifier == 11312) prodCells11312++;
	if(tripletCollection[it].tripletIdentifier == 11531) prodCells11531++;
	if(tripletCollection[it].tripletIdentifier == -10079) prodCells_10079++;

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
    
    if (m_debug > 1) std::cout<<" -- > ForwardCA done!"<<std::endl;
    
    
//for time measurements (developer version)
    auto t_caforward_end = std::chrono::system_clock::now();
    auto t_caforwardE = std::chrono::duration_cast<std::chrono::microseconds>(t_caforward_end - t_caforward);
    auto t_cabackward = std::chrono::system_clock::now();


//Connect triplets into multiplets Backward CA (longest seeds)
    BackwardCA(max_status);
	if (m_debug > 1) std::cout<<" -- > BackwardCA for 5 hits done!"<<std::endl;
	
	auto t_cabackward_end = std::chrono::system_clock::now();
    auto t_cabackwardE = std::chrono::duration_cast<std::chrono::microseconds>(t_cabackward_end - t_cabackward);

	
//Connect the remaining triplets into multiplets with one hit less
    auto t_cabackward4 = std::chrono::system_clock::now();

if (m_secondbkg){
    BackwardCA(max_status -1);
	if (m_debug > 1) std::cout<<" -- > BackwardCA for 4 hits done!"<<std::endl;
}
     auto t_cabackward4_end = std::chrono::system_clock::now();
    auto t_cabackward4E = std::chrono::duration_cast<std::chrono::microseconds>(t_cabackward4_end - t_cabackward4);
    
//Final Forward Cleaning (again based on pt)
    auto t_secondforward = std::chrono::system_clock::now();
    
    
if(do_finalfor == true)
    FinalForwardCleaning();
else{			//if not final_multi_collection becomes a dummy copy of multiplets
	for(size_t im = 0; im < multiplets.size(); im++)
		final_multi_collection.push_back(im);
}  
    
    if (m_debug > 1) std::cout<<" -- > Second forward done!"<<std::endl;

    
    auto t_secondforward_end = std::chrono::system_clock::now();
    auto t_secondforwardE = std::chrono::duration_cast<std::chrono::microseconds>(t_secondforward_end - t_secondforward);
    

    auto t_multprod = std::chrono::system_clock::now();


int nSeeds = 0;
int nSeeds5h = 0;
int nSeeds4h = 0;
 	
 	

for(size_t im = 0; im < final_multi_collection.size(); im++){
	int m_index = final_multi_collection[im];
    std::vector<ConstRecHitPointer>  multiHitPointer;
	for(size_t ii = 0; ii < multiplets[m_index].size(); ii++)
		multiHitPointer.push_back(multiplets[m_index][ii]);
	if(multiHitPointer.size()==5)
	{
    		SeedingHitSet multiset(multiHitPointer[4], multiHitPointer[3], multiHitPointer[2], multiHitPointer[1], multiHitPointer[0]);
			resultCA.push_back(multiset);
			nSeeds++;
			nSeeds5h++;
			GlobalPoint p_last = multiHitPointer[0]->globalPosition();
			double pl_eta = p_last.eta();
			penteta.push_back(pl_eta);
	}
	else 	if(multiHitPointer.size()==4)
	{	
    		SeedingHitSet multiset(multiHitPointer[3], multiHitPointer[2], multiHitPointer[1], multiHitPointer[0]);
			resultCA.push_back(multiset);
			nSeeds++;
			nSeeds4h++;			
			GlobalPoint p_last = multiHitPointer[0]->globalPosition();
			double pl_eta = p_last.eta();
			penteta.push_back(pl_eta);		
	}

}


if (m_debug > 1) std::cout<<" -- > Multiset done!"<<std::endl;


//fill the penteta tree
if(m_tree){ 
	for (size_t itree = 0; itree<penteta.size(); itree++) {
       seedhist->Fill(penteta[itree]);
    }

 }	
//for time measurements (developer version)

   auto t_multprod_end = std::chrono::system_clock::now();
   auto t_multprodE = std::chrono::duration_cast<std::chrono::microseconds>(t_multprod_end - t_multprod);
   auto t_total_end = std::chrono::system_clock::now();
   auto t_totalE = std::chrono::duration_cast<std::chrono::microseconds>(t_total_end - t_total);


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
tmp_prodCells_10079 = prodCells_10079;
tmp_prodLay1Cells = prodCells123eta + prodCells10121 + prodCells_10079;
tmp_neigCells = ncont; 
tmp_etaCValue = dEta_cut;  
tmp_fitCells = fittedTripletCollection.size();  
tmp_nSteps = n_fiter;
tmp_nSeeds = nSeeds; 
tmp_nSeeds5h = nSeeds5h; 
tmp_nSeeds4h = nSeeds4h; 
  
//fill timing tree

	tmp_t_total = t_totalE.count();
	tmp_t_triplets = t_tripletsE.count();
	tmp_t_casetup = t_casetupE.count();
	tmp_t_caforward = t_caforwardE.count();
	tmp_t_cabackward = t_cabackwardE.count();
	tmp_t_cabackward4 = t_cabackward4E.count();
	tmp_t_secondforward = t_secondforwardE.count();
	tmp_t_multprod = t_multprodE.count();

evtree->Fill();
}

if(m_debug > 0){
std::cout<<"================================"<<std::endl;
std::cout<<"Processed "<<scoll->size()<<" seeds"<<std::endl;
std::cout<<"produced  "<<tripletCollection.size()<<" CAcells, of whom:"<<std::endl;
std::cout<<ncont<<" have neighbors"<<std::endl;
std::cout<<fittedTripletCollection.size()<<" have been successfully fitted!"<<std::endl;
std::cout<<"Forward iterations took n_steps = "<<n_fiter<<std::endl;
std::cout<<"Multiset before final cleaning = "<<multiplets.size()<<std::endl;
std::cout<<"Multiset after final cleaning = "<<final_multi_collection.size()<<" (Ratio = "<< final_multi_collection.size()/(double)multiplets.size()<<")"<<std::endl;
std::cout<<"Final number of multiSeeds = "<<resultCA.size()<<std::endl;
std::cout<<"================================"<<std::endl;
}


//clear all
    
delete scoll;
tripletCollection.clear();
std::vector<CAcell>().swap(tripletCollection);
fittedTripletCollection.clear();
std::vector<CAcell*>().swap(fittedTripletCollection);
hitUsage.clear();

multiplets.clear();
final_multi_collection.clear();

penteta.clear();

eventcounter++;
}

void CAHitsGeneratorForDebuggers::CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder> & builder, SeedingHitSet sset,  int seednumber){

   	CAcell trip;
    //translate hits into a SeedingHitSet
    
    size_t i = 0;
    
    std::vector<const BaseTrackerRecHit*> starting_triplet;
    
    while(i < sset.size()){
    	if(sset[i]!=nullptr)  starting_triplet.push_back(sset.get(i));
    	i++;
    }
      

    
	trip.hits = starting_triplet;
    
    auto hi0 = starting_triplet[0];
    auto hi1 = starting_triplet[1];
    auto hi2 = starting_triplet[2];


      trip.rawId0 = OmniRef(hi0);
      trip.rawId1 = OmniRef(hi1);
      trip.rawId2 = OmniRef(hi2);


    GlobalPoint p0 = hi0->globalPosition();
    GlobalPoint p1 = hi1->globalPosition();
    GlobalPoint p2 = hi2->globalPosition();

    trip.p0 = p0;	
    trip.p1 = p1;	   
    trip.p2 = p2;

    trip.dEta = DeltaEta(p0,p2);
    
    LayerClassification hit0cl(hi0->geographicalId());
    LayerClassification hit1cl(hi1->geographicalId());
    LayerClassification hit2cl(hi2->geographicalId());
    trip.hitId0 = hit0cl.getLayer();
    trip.hitId1 = hit1cl.getLayer();
    trip.hitId2 = hit2cl.getLayer();


    trip.eventNumber = eventcounter;
    trip.tripletIdentifier = trip.hitId0+10*trip.hitId1+100*trip.hitId2;

    
    if (m_debug>4) {
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

    tripletCollection.push_back(trip);


return;
}

void CAHitsGeneratorForDebuggers::fitTripletSeeds(CAcell *trip, edm::ESHandle<MagneticField> & mf) {

 
    
    double nomField = mf->nominalValue();

    FastHelix fit(trip->p2, trip->p1, trip->p0, nomField , mf.product());
    
    GlobalTrajectoryParameters params = fit.stateAtVertex();
    
    if (fit.isValid()){
    trip->FitSuccessful = true;
	trip->LocalMomentum = params.momentum().mag();
    }

return;
}



int CAHitsGeneratorForDebuggers::ForwardCA(){
    int Stop = -1;
    int step_iterator = 0;
    
    int it_max = 20; //Maximum number of iteraitons
        
    while (Stop < (int)fittedTripletCollection.size() && step_iterator<it_max) { 
        step_iterator++;
        Stop = 0;
        
        
        for (size_t t =0; t<fittedTripletCollection.size(); t++) {
            
            if (fittedTripletCollection[t]->left_neighbors.size()!=0) {
               bool en_cont = false;                       //useful variable - moved from 2 lines behind 
               for (size_t il = 0; il<fittedTripletCollection[t]->left_neighbors.size(); il++) {
                    if (fittedTripletCollection[t]->CAstatus <= fittedTripletCollection[t]->left_neighbors[il]->CAstatus && en_cont == false){
                        	(fittedTripletCollection[t]->CAstatus)++;
                      	 	en_cont = true;
                    			}
                }
                if (en_cont == false)  Stop++;
            }
            else  Stop++;		//i.e. if fittedTripletCollection[t] does not have left neighbors
            
        }   //end of loop on triplets      
        
    }		//end of while
    
    //debug
if(m_debug>2){
    for (size_t t =0; t<fittedTripletCollection.size(); t++) {
        std::cout<<" ._._._._._._._._ "<<std::endl;
        std::cout<<" cell number "<<t<<std::endl;
        std::cout<<" has = "<<fittedTripletCollection[t]->IsNeighborly<<" neighbors. Is in lays: "<<fittedTripletCollection[t]->tripletIdentifier<<std::endl;
        std::cout<<"left_list size == "<<fittedTripletCollection[t]->left_neighbors.size()<<std::endl;
       	std::cout<<"left neighbours (Id, eta) == {";
        for(size_t k=0; k< fittedTripletCollection[t]->left_neighbors.size(); k++){
        	std::cout<<"(";
        	std::cout<<(fittedTripletCollection[t]->left_neighbors[k])->tripletIdentifier;
        	std::cout<<(fittedTripletCollection[t]->left_neighbors[k])->dEta;        	
        	std::cout<<"),";
        }
        std::cout<<"}"<<std::endl;
    	std::cout<<"right_list size == "<<fittedTripletCollection[t]->right_neighbors.size()<<std::endl;
    	std::cout<<"right neighbours (Id, eta) == {";
        for(size_t k=0; k< fittedTripletCollection[t]->right_neighbors.size(); k++){
        	std::cout<<"(";
        	std::cout<<(fittedTripletCollection[t]->right_neighbors[k])->tripletIdentifier;
        	std::cout<<(fittedTripletCollection[t]->right_neighbors[k])->dEta;        	
        	std::cout<<"),";
        }
        std::cout<<"}"<<std::endl;
        std::cout<<" processed with status: "<<fittedTripletCollection[t]->CAstatus<<std::endl;
    }
    }
    
    if(step_iterator>=it_max)
        std::cout<<"WARNING: too many CA iterations!"<<std::endl;
    
    return step_iterator;
}


void CAHitsGeneratorForDebuggers::BackwardCA(int m_status){

int multi_id = 1;

for (size_t t =0; t<fittedTripletCollection.size(); t++) {
//start only with non-used triplets with maximum CAstatus
if (m_debug > 2) std::cout<<"Processing triplet no. "<<t<<" for the backwardCA"<<std::endl;
if(fittedTripletCollection[t]->CAstatus == m_status && fittedTripletCollection[t]->IsUsed==0){
	CAcell *current_cell = fittedTripletCollection[t];
	std::vector<ConstRecHitPointer> multicontainer;
	while(current_cell){
		current_cell = define_used(current_cell, multi_id, &multicontainer);
		}
		multi_id++;	
	if (m_debug > 2) std::cout<<"Multicontainer size =  "<<multicontainer.size()<<std::endl;
	multiplets.push_back(multicontainer);
	fittedTripletCollection[t]->multilabel = multiplets.size() -1;
	multicontainer.clear();
	}
}

return;
}

CAcell* CAHitsGeneratorForDebuggers::define_used(CAcell *cell, int id, std::vector<ConstRecHitPointer>*multicontainer){
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


void CAHitsGeneratorForDebuggers::FinalForwardCleaning(){
	
for (size_t t =0; t<fittedTripletCollection.size(); t++) {
    //std::cout<<"in fittedTripletCollection with status = "<<fittedTripletCollection[t]->CAstatus<<" and used by "<<fittedTripletCollection[t]->IsUsed<<std::endl;
	if (fittedTripletCollection[t]->CAstatus == 1 && fittedTripletCollection[t]->IsUsed!=0){ //if it starts a multiplet...
		CAcell *current_cell = fittedTripletCollection[t];
		while(current_cell->multilabel==-1){
            double pt_tmp = 99999.0;
			CAcell *next = new CAcell();
			for(size_t il =0; il<current_cell->right_neighbors.size(); il++){
					if(current_cell->tripletIdentifier != IsAtLeft(current_cell->right_neighbors[il]->tripletIdentifier)) 				continue;
					double pt_diff = fabs(current_cell->getLocalMomentum() - current_cell->right_neighbors[il]->getLocalMomentum());
					if(pt_diff <= pt_tmp){
								pt_tmp = pt_diff;
								next = current_cell->right_neighbors[il];	 
					}
			}
			current_cell = next;
			}	//out of the while
        final_multi_collection.push_back(current_cell->multilabel);
			}
		}//out of the loop
return;
}


//Useful tools to store triplets in root files
void CAHitsGeneratorForDebuggers::fittreegenerator(){

//triplets tree
triptree->Branch("Event",&tmp_event,"Event/I");
triptree->Branch("LocalMomentum",&tmp_LocalMomentum,"LocalMomentum/D");

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
evtree->Branch("prodLay1Cells",&tmp_prodLay1Cells,"prodLay1Cells/I");
evtree->Branch("prodCells123",&tmp_prodCells123,"prodCells123/I");
evtree->Branch("prodCells123eta",&tmp_prodCells123eta,"prodCells123eta/I");
evtree->Branch("prodCells234",&tmp_prodCells234,"prodCells234/I");
evtree->Branch("prodCells345",&tmp_prodCells345,"prodCells345/I");
evtree->Branch("prodCells10343",&tmp_prodCells10343,"prodCells10343/I");
evtree->Branch("prodCells10121",&tmp_prodCells10121,"prodCells10121/I");
evtree->Branch("prodCells11312",&tmp_prodCells11312,"prodCells11312/I");
evtree->Branch("prodCells11531",&tmp_prodCells11531,"prodCells11531/I");
evtree->Branch("prodCells_10079",&tmp_prodCells_10079,"prodCells_10079/I");
evtree->Branch("neigCells",&tmp_neigCells,"neigCells/I");
evtree->Branch("etaCValue",&tmp_etaCValue,"etaCValue/D");
evtree->Branch("fitCells",&tmp_fitCells,"fitCells/I");
evtree->Branch("nSteps",&tmp_nSteps,"nSteps/I");
evtree->Branch("nSeeds",&tmp_nSeeds,"nSeeds/I");
evtree->Branch("nSeeds5h",&tmp_nSeeds5h,"nSeeds5h/I");
evtree->Branch("nSeeds4h",&tmp_nSeeds4h,"nSeeds4h/I");
evtree->Branch("nInvalid",&invalidseqcounter,"nInvalid/I");

    
//time variables in event tree
evtree->Branch("t_total",&tmp_t_total , "t_total/I");
evtree->Branch("t_triplets",&tmp_t_triplets , "t_triplets/I");
evtree->Branch("t_casetup",&tmp_t_casetup, "t_casetup/I");
evtree->Branch("t_caforward",&tmp_t_caforward  , "t_caforward/I");
evtree->Branch("t_cabackward",&tmp_t_cabackward , "t_cabackward/I");
evtree->Branch("t_cabackward4",&tmp_t_cabackward4 , "t_cabackward4/I");
evtree->Branch("t_secondforward",&tmp_t_secondforward , "t_secondforward/I");
evtree->Branch("t_multiprod",&tmp_t_multprod  , "t_multprod/I");

    
return;
}


void CAHitsGeneratorForDebuggers::fillfittree(CAcell *t){

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

std::vector<CAcell *> CAHitsGeneratorForDebuggers::ListIntersect(std::list<CAcell *> list_h1, std::list<CAcell *> list_h2){
    
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



unsigned int CAHitsGeneratorForDebuggers::OmniRef(const BaseTrackerRecHit* rhit){

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


double CAHitsGeneratorForDebuggers::DeltaEta(GlobalPoint p1, GlobalPoint p2){
    
    GlobalVector gV(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
    double eta = gV.eta();
    
    return eta;
}

//check is two triplets can be left-right joined 
//NB: It should be changed anlytime the Layer configuration is changed (This works for A-B-C-D or subsets)
int CAHitsGeneratorForDebuggers::IsAtLeft(int identif){
	int ret_val = 0;
	switch (identif){
		case 321:
			ret_val = 0;
			break;
		case 10121:
			ret_val = 0;
			break;
		case -10079:
			ret_val = 0;
			break;
		case 432:
			ret_val = 321;
			break;
		case 543:
			ret_val = 432;
			break;
		case 10132:
			ret_val = 321;
			break;
		case 10343:
			ret_val = 432;
			break;
		case 11313:
			ret_val = 10132;
			break;
		case -10068:
			ret_val = 321;
			break;
		case -10257:
			ret_val = 432;
			break;
		case -11307:
			ret_val = -10068;
			break;
		case 11312:
			ret_val = 10121;
			break;
		case -11308:
			ret_val = -10079;
			break;
		case 11531:
			ret_val = 11312;
			break;
		case -11531:
			ret_val = -11308;
			break;
		default:
			std::cout<<"WARNING: triplet Id: "<<identif<<" not recognized!"<<std::endl;	
			ret_val = -1;
			break;			
	}
	return ret_val;
}


//map of the possible right neighbors (not used up to now)
std::vector<int> CAHitsGeneratorForDebuggers::IsAtRight(int identif){
std::vector<int> ret_val;
//as for lefties, but right neighbors can be more than one
	switch (identif){
		case 321:
			ret_val.push_back(432);
			ret_val.push_back(10132);
			ret_val.push_back(-10068);
			break;
		case 10121:
			ret_val.push_back(11312);
			break;
		case -10079:
			ret_val.push_back(-11308);		
			break;
		case 432:
			ret_val.push_back(543);
			ret_val.push_back(10343);
			ret_val.push_back(-10257);
			break;
		case 543:
			ret_val.push_back(0);
			break;
		case 10132:
			ret_val.push_back(11313);
			break;
		case 10343:
			ret_val.push_back(0);
			break;
		case 11313:
			ret_val.push_back(0);
			break;
		case -10068:
			ret_val.push_back(-11307);
			break;
		case -10257:
			ret_val.push_back(0);
			break;
		case -11307:
			ret_val.push_back(0);
			break;
		case 11312:
			ret_val.push_back(11531);		
			break;
		case -11308:
			ret_val.push_back(-11531);		
			break;
		case 11531:
			ret_val.push_back(0);		
			break;
		case -11531:
			ret_val.push_back(0);
			break;
		default:
			std::cout<<"WARNING: triplet Id: "<<identif<<" not recognized!"<<std::endl;	
			ret_val.push_back(-1);
			break;			
	}

return ret_val;	
}
