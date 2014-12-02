// -*- C++ -*-
// Package:    CAHitsGenerator
// Class:      CAHitsGenerator
//
/**\class CAHitsGenerator CAHitsGenerator.cc RecoTracker/CAtracker/plugins/CAHitsGenerator.cc*/
//
//  Author:  Gigi Cappello
//  Created:  Wed, 22 Lug 2014
//

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
    m_builderName = cfg.getUntrackedParameter <std::string> ("Builder");
    dEta_cut = cfg.getParameter<double>("EtaCut");
    m_secondbkg = cfg.getParameter<bool>("MakeSecondBackward");
    do_finalfor = cfg.getParameter<bool>("MakeFinalForward");
    
    //triplet producer definition
    oPset = cfg.getParameter<edm::ParameterSet>("GeneratorPSet");
    oName = oPset.getParameter<std::string>("ComponentName");
    
    theGenerator.reset(HitTripletGeneratorFromPairAndLayersFactory::get()->create(oName, oPset, iC));
    theGenerator->init(HitPairGeneratorFromLayerPair( 0, 1, &theLayerCache), &theLayerCache);
    
    eventcounter = 0;
    max_status = 3;
    
    initialised = false;
}


CAHitsGenerator::~CAHitsGenerator()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void CAHitsGenerator::init(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
    iSetup.get<TransientRecHitRecord>().get(m_builderName, builder);
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);
    
    initialised = true;
    
    return;
}


void CAHitsGenerator::hitSets(const TrackingRegion& region, OrderedMultiHits & resultCA, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    using namespace edm;
    
    
    //for time measurements (developer version)
    
    
    if(!initialised) init(iEvent,iSetup);
    
    //BASE ONE
    OrderedHitTriplets *scoll = new OrderedHitTriplets();
    
    //to be changed
    scoll->reserve(10000);
    
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
        
        trilaycount++;
    }
    
    
    //produce triplets
    
    theLayerCache.clear();
    
    
    tripletCollection.reserve(scoll->size());
    fittedTripletCollection.reserve(scoll->size());
    
    //Neighborhood map definition
    for (size_t is = 0; is<scoll->size(); is++) {
        CAcellGenerator(builder, (*scoll)[is], is);
        hitUsage[tripletCollection[is].rawId0].push_back(&tripletCollection[is]);
        hitUsage[tripletCollection[is].rawId1].push_back(&tripletCollection[is]);
        hitUsage[tripletCollection[is].rawId2].push_back(&tripletCollection[is]);
        
    }
    
    
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
        
        
        IsNeighbor = tripletCollection[it].left_neighbors.size()+tripletCollection[it].right_neighbors.size();
        
        tripletCollection[it].IsNeighborly = IsNeighbor;
        
        list_h0.clear();
        list_h1.clear();
        list_h2.clear();
        j_list01.clear();
        j_list12.clear();
        
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
    

    
    //Loops on fittedtriplets -> Forward CA
    int n_fiter = -1;
    n_fiter = ForwardCA();
    
    if (m_debug > 1) std::cout<<" -- > ForwardCA done!"<<std::endl;

    //Connect triplets into multiplets Backward CA (longest seeds)
    BackwardCA(max_status);
    if (m_debug > 1) std::cout<<" -- > BackwardCA for 5 hits done!"<<std::endl;
    
    
    //Connect the remaining triplets into multiplets with one hit less
    if (m_secondbkg){
    BackwardCA(max_status -1);
    if (m_debug > 1) std::cout<<" -- > BackwardCA for 4 hits done!"<<std::endl;
    }
    
    //Final Forward Cleaning (again based on pt)
    
    
    if(do_finalfor == true)
    FinalForwardCleaning();
    else{			//if not final_multi_collection becomes a dummy copy of multiplets
        for(size_t im = 0; im < multiplets.size(); im++)
        final_multi_collection.push_back(im);
    }
    
    if (m_debug > 1 && do_finalfor == true) std::cout<<" -- > Second forward done!"<<std::endl;
    
    
    for(size_t im = 0; im < final_multi_collection.size(); im++){
        int m_index = final_multi_collection[im];
        std::vector<ConstRecHitPointer>  multiHitPointer;
        for(size_t ii = 0; ii < multiplets[m_index].size(); ii++)
        multiHitPointer.push_back(multiplets[m_index][ii]);
        if(multiHitPointer.size()==5)
        {
            SeedingHitSet multiset(multiHitPointer[4], multiHitPointer[3], multiHitPointer[2], multiHitPointer[1], multiHitPointer[0]);
            resultCA.push_back(multiset);
   
         }
        else 	if(multiHitPointer.size()==4)
        {
            SeedingHitSet multiset(multiHitPointer[3], multiHitPointer[2], multiHitPointer[1], multiHitPointer[0]);
            resultCA.push_back(multiset);
           
            
        }
        
    }
    
    
    if (m_debug > 1) std::cout<<" -- > Multiset done!"<<std::endl;
    
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
    
    eventcounter++;
}

void CAHitsGenerator::CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder> & builder, SeedingHitSet sset,  int seednumber){
    
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


void CAHitsGenerator::BackwardCA(int m_status){
    
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


void CAHitsGenerator::FinalForwardCleaning(){
    
    for (size_t t =0; t<fittedTripletCollection.size(); t++) {
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

//check is two triplets can be left-right joined 
//NB: It should be changed anlytime the Layer configuration is changed (This works for A-B-C-D or subsets)
int CAHitsGenerator::IsAtLeft(int identif){
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
std::vector<int> CAHitsGenerator::IsAtRight(int identif){
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
