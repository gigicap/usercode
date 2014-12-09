// -*- C++ -*-
// Package:    CAGeneratorFromPairs
// Class:      CAGeneratorFromPairs (Pair version)
//
/**\class CAGeneratorFromPairs CAGeneratorFromPairs.cc RecoTracker/CAtracker/plugins/CAGeneratorFromPairs.cc*/
//
//  Author:  Gigi Cappello
//  Created:  Wed, 22 Lug 2014
//

// system include files
#include "CAGeneratorFromPairs.h"  //class CAcell defined here
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include <chrono>
#include <math.h>


using namespace ctfseeding;


CAGeneratorFromPairs::CAGeneratorFromPairs(const edm::ParameterSet& cfg , edm::ConsumesCollector& iC) :
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
    
    theGenerator.reset(new HitPairGeneratorFromLayerPair(0, 1, &theLayerCache, 0, theMaxElement));
    
    eventcounter = 0;
    max_status = 3;
    //max_status = number of cells to join
    //max_status = n - m + 1 (where m = number of hits in a cell, n = number of hits in the final seed)
    
    initialised = false;
}


CAGeneratorFromPairs::~CAGeneratorFromPairs()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void CAGeneratorFromPairs::init(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    
    iSetup.get<TrackerDigiGeometryRecord>().get(tracker);
    iSetup.get<TransientRecHitRecord>().get(m_builderName, builder);
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);
    
    initialised = true;
    
    return;
}


void CAGeneratorFromPairs::hitSets(const TrackingRegion& region, OrderedMultiHits & resultCA, const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    using namespace edm;
    
    
    
    if(!initialised) init(iEvent,iSetup);
    
    //BASE ONE
    OrderedHitPairs *scoll = new OrderedHitPairs();
    
    //to be changed
    scoll->reserve(50000);
    
    //use the beamspot instead?
    vertex = region->origin();
    
    //For pairs
    edm::Handle<SeedingLayerSetsHits> hlayers;
    ev.getByToken(theSeedingLayerToken, hlayers);
    
    const SeedingLayerSetsHits& layers = *hlayers;
    
    if(layers.numberOfLayersInSet() != 2)
    throw cms::Exception("Configuration") << "CombinedHitPairGenerator expects SeedingLayerSetsHits::numberOfLayersInSet() to be 2, got " << layers.numberOfLayersInSet();
    
    
    int pairlaycount = 0;
    for(SeedingLayerSetsHits::SeedingLayerSet layerSet: layers) {
             theGenerator->setSeedingLayers(layerSet);
             theGenerator->hitPairs( region, *scoll, iEvent, iSetup);
        
            pairlaycount++;
           }
    
    
    theLayerCache.clear();
    
    
    pairCollection.reserve(scoll->size());
    fittedPairCollection.reserve(scoll->size());
    
    //Neighborhood map definition
    for (size_t is = 0; is<scoll->size(); is++) {
        CAcellGenerator(builder, (*scoll)[is], is);
        hitUsage[pairCollection[is].rawId0].push_back(&pairCollection[is]);
        hitUsage[pairCollection[is].rawId1].push_back(&pairCollection[is]);
    }
    
    
    //Loop over the triplets to build the neighborhood map
    //i.e. fill the "HasNeighbor" var and IsNeighborly
    int ncont = 0;
    for (size_t it = 0; it<pairCollection.size(); it++) {
        int IsNeighbor = 0;
        
        
        std::list<CAcell *> j_list0 = hitUsage[pairCollection[it].rawId0];
        std::list<CAcell *> j_list1 = hitUsage[pairCollection[it].rawId1];
        

        for (size_t il = 0; il<j_list0.size(); il++){
            //added to avoid duplicates in layer superimpositions
            if((j_list0[il] != &pairCollection[it]) && (j_list0[il]->doubletIdentifier == IsAtLeft(pairCollection[it].doubletIdentifier))){
                if((j_list0[il]->doubletIdentifier != pairCollection[it].doubletIdentifier) && fabs(j_list0[il]->dEta - pairCollection[it].dEta)< dEta_cut)
                pairCollection[it].left_neighbors.push_back(j_list0[il]);
            }
        }
        
        
        for (size_t il = 0; il<j_list1.size(); il++){
            if(j_list1[il] != &pairCollection[it] && (j_list1[il]->doubletIdentifier != pairCollection[it].doubletIdentifier) && fabs(j_list1[il]->dEta - pairCollection[it].dEta)< dEta_cut)
            pairCollection[it].right_neighbors.push_back(j_list1[il]);
        }
        
        
        IsNeighbor = pairCollection[it].left_neighbors.size()+pairCollection[it].right_neighbors.size();
        
        pairCollection[it].IsNeighborly = IsNeighbor;
        
        j_list0.clear();
        j_list1.clear();
        
        //fit only the triplets with neighbors
        if (pairCollection[it].IsNeighborly != 0) {
            ncont++;
            // fast helix fitter
            fitPair(&pairCollection[it], vertex, theMF);    //add the vertex calculation
            //store the neighborly and fitted triplets
            if (pairCollection[it].FitSuccessful) {
                pairCollection[it].CAstatus = 1;     //initialize CAstatus
                fittedPairCollection.push_back(&pairCollection[it]);
            }
        }
        
        
    }
    

    
    //Loops on fittedtriplets -> Forward CA
    int n_fiter = -1;
    n_fiter = ForwardCA();
    
    if (m_debug > 1) std::cout<<" -- > ForwardCA done!"<<std::endl;

    //Connect triplets into multiplets Backward CA (longest seeds)
    BackwardCA(max_status);
    if (m_debug > 1) std::cout<<" -- > BackwardCA for longest seeds done!"<<std::endl;
    
    
    //Connect the remaining triplets into multiplets with one hit less
    if (m_secondbkg){
    BackwardCA(max_status -1);
    if (m_debug > 1) std::cout<<" -- > BackwardCA for shorter seeds done!"<<std::endl;
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
        std::cout<<"Processed "<<scoll->size()<<" objs"<<std::endl;
        std::cout<<"produced  "<<pairCollection.size()<<" CAcells, of whom:"<<std::endl;
        std::cout<<ncont<<" have neighbors"<<std::endl;
        std::cout<<fittedPairCollection.size()<<" have been successfully fitted!"<<std::endl;
        std::cout<<"Forward iterations took n_steps = "<<n_fiter<<std::endl;
        std::cout<<"Multiset before final cleaning = "<<multiplets.size()<<std::endl;
        std::cout<<"Multiset after final cleaning = "<<final_multi_collection.size()<<" (Ratio = "<< final_multi_collection.size()/(double)multiplets.size()<<")"<<std::endl;
        std::cout<<"Final number of multiSeeds = "<<resultCA.size()<<std::endl;
        std::cout<<"================================"<<std::endl;
    }
    
    
    //clear all
    
    delete scoll;
    pairCollection.clear();
    std::vector<CAcell>().swap(pairCollection);
    fittedPairCollection.clear();
    std::vector<CAcell*>().swap(fittedpairCollection);
    hitUsage.clear();
    
    multiplets.clear();
    final_multi_collection.clear();
    
    eventcounter++;
}

void CAGeneratorFromPairs::CAcellGenerator(edm::ESHandle<TransientTrackingRecHitBuilder> & builder, SeedingHitSet sset,  int seednumber){
    
   	CAcell pr;
    //translate hits into a SeedingHitSet
    
    size_t i = 0;
    
    std::vector<const BaseTrackerRecHit*> starting_pair;
    
    while(i < sset.size()){
        if(sset[i]!=nullptr)  starting_pair.push_back(sset.get(i));
        i++;
    }
    
    
    
    pr.hits = starting_pair;
    
    auto hi0 = starting_pair[0];
    auto hi1 = starting_pair[1];
    
    
    pr.rawId0 = OmniRef(hi0);
    pr.rawId1 = OmniRef(hi1);
    
    
    GlobalPoint p0 = hi0->globalPosition();
    GlobalPoint p1 = hi1->globalPosition();
    
    pr.p0 = p0;
    pr.p1 = p1;
    
    pr.dEta = DeltaEta(p0,p1);
    
    LayerClassification hit0cl(hi0->geographicalId());
    LayerClassification hit1cl(hi1->geographicalId());
    pr.hitId0 = hit0cl.getLayer();
    pr.hitId1 = hit1cl.getLayer();
    
    
    pr.eventNumber = eventcounter;
    pr.doubletIdentifier = pr.hitId0+10*pr.hitId1;
    
    
    if (m_debug>4) {
        std::cout << " layer 0 =  " << pr.hitId0 << std::endl;
        std::cout << " layer 1 =  " << pr.hitId1 << std::endl;
        std::cout << " identifier = "<< pr.doubletIdentifier <<std::endl;
    }
    
    
    pr.DoubletNumber = seednumber;
    pr.FitSuccessful = false;
    pr.LocalMomentum = 0;
    
    //Is Used (for the final Backwards CA)
    pr.IsUsed = 0;
    
    pairCollection.push_back(pr);
    
    
    return;
}

//use the vertex for the pair
void CAGeneratorFromPairs::fitPair(CAcell *cell, GlobalPoint vertex, edm::ESHandle<MagneticField> & mf) {
    
    double nomField = mf->nominalValue();
    
    FastHelix fit(cell->p1, cell->p0, vertex, nomField , mf.product());
    
    GlobalTrajectoryParameters params = fit.stateAtVertex();
    
    if (fit.isValid()){
        pr->FitSuccessful = true;
        pr->LocalMomentum = params.momentum().mag();
    }
    
    return;
}



int CAGeneratorFromPairs::ForwardCA(){
    int Stop = -1;
    int step_iterator = 0;
    
    int it_max = 20; //Maximum number of iteraitons
    
    while (Stop < (int)fittedPairCollection.size() && step_iterator<it_max) {
        step_iterator++;
        Stop = 0;
        
        
        for (size_t t =0; t<fittedPairCollection.size(); t++) {
            
            if (fittedPairCollection[t]->left_neighbors.size()!=0) {
                bool en_cont = false;                       //useful variable - moved from 2 lines behind
                for (size_t il = 0; il<fittedPairCollection[t]->left_neighbors.size(); il++) {
                    if (fittedPairCollection[t]->CAstatus <= fittedPairCollection[t]->left_neighbors[il]->CAstatus && en_cont == false){
                        (fittedPairCollection[t]->CAstatus)++;
                        en_cont = true;
                    }
                }
                if (en_cont == false)  Stop++;
            }
            else  Stop++;		//i.e. if fittedPairCollection[t] does not have left neighbors
            
        }   //end of loop on pairs
        
    }		//end of while
    
    //debug
    if(m_debug>2){
        for (size_t t =0; t<fittedPairCollection.size(); t++) {
            std::cout<<" ._._._._._._._._ "<<std::endl;
            std::cout<<" cell number "<<t<<std::endl;
            std::cout<<" has = "<<fittedPairCollection[t]->IsNeighborly<<" neighbors. Is in lays: "<<fittedPairCollection[t]->doubletIdentifier<<std::endl;
            std::cout<<"left_list size == "<<fittedPairCollection[t]->left_neighbors.size()<<std::endl;
            std::cout<<"left neighbours (Id, eta) == {";
            for(size_t k=0; k< fittedPairCollection[t]->left_neighbors.size(); k++){
                std::cout<<"(";
                std::cout<<(fittedPairCollection[t]->left_neighbors[k])->doubletIdentifier;
                std::cout<<(fittedPairCollection[t]->left_neighbors[k])->dEta;
                std::cout<<"),";
            }
            std::cout<<"}"<<std::endl;
            std::cout<<"right_list size == "<<fittedPairCollection[t]->right_neighbors.size()<<std::endl;
            std::cout<<"right neighbours (Id, eta) == {";
            for(size_t k=0; k< fittedPairCollection[t]->right_neighbors.size(); k++){
                std::cout<<"(";
                std::cout<<(fittedPairCollection[t]->right_neighbors[k])->doubletIdentifier;
                std::cout<<(fittedPairCollection[t]->right_neighbors[k])->dEta;
                std::cout<<"),";
            }
            std::cout<<"}"<<std::endl;
            std::cout<<" processed with status: "<<fittedPairCollection[t]->CAstatus<<std::endl;
        }
    }
    
    if(step_iterator>=it_max)
    std::cout<<"WARNING: too many CA iterations!"<<std::endl;
    
    return step_iterator;
}


void CAGeneratorFromPairs::BackwardCA(int m_status){
    
    int multi_id = 1;
    
    for (size_t t =0; t<fittedPairCollection.size(); t++) {
        //start only with non-used triplets with maximum CAstatus
        if (m_debug > 2) std::cout<<"Processing triplet no. "<<t<<" for the backwardCA"<<std::endl;
        if(fittedPairCollection[t]->CAstatus == m_status && fittedPairCollection[t]->IsUsed==0){
            CAcell *current_cell = fittedPairCollection[t];
            std::vector<ConstRecHitPointer> multicontainer;
            while(current_cell){
                current_cell = define_used(current_cell, multi_id, &multicontainer);
            }
            multi_id++;
            if (m_debug > 2) std::cout<<"Multicontainer size =  "<<multicontainer.size()<<std::endl;
            multiplets.push_back(multicontainer);
            fittedPairCollection[t]->multilabel = multiplets.size() -1;
            multicontainer.clear();
        }
    }
    
    return;
}

CAcell* CAGeneratorFromPairs::define_used(CAcell *cell, int id, std::vector<ConstRecHitPointer>*multicontainer){
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


void CAGeneratorFromPairs::FinalForwardCleaning(){
    
    for (size_t t =0; t<fittedPairCollection.size(); t++) {
        if (fittedPairCollection[t]->CAstatus == 1 && fittedPairCollection[t]->IsUsed!=0){ //if it starts a multiplet...
            CAcell *current_cell = fittedPairCollection[t];
            while(current_cell->multilabel==-1){
                double pt_tmp = 99999.0;
                CAcell *next = new CAcell();
                for(size_t il =0; il<current_cell->right_neighbors.size(); il++){
                    if(current_cell->doubletIdentifier != IsAtLeft(current_cell->right_neighbors[il]->doubletIdentifier)) 				continue;
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


unsigned int CAGeneratorFromPairs::OmniRef(const BaseTrackerRecHit* rhit){
    
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


double CAGeneratorFromPairs::DeltaEta(GlobalPoint p1, GlobalPoint p2){
    
    GlobalVector gV(p2.x()-p1.x(),p2.y()-p1.y(),p2.z()-p1.z());
    double eta = gV.eta();
    
    return eta;
}

//check is two triplets can be left-right joined 
//NB: It should be changed anlytime the Layer configuration is changed (This works for A-B-C-D or subsets)
int CAGeneratorFromPairs::IsAtLeft(int identif){
    int ret_val = 0;
    switch (identif){
        case 21:
            ret_val = 0;
            break;
        case 32:
            ret_val = 21;
            break;
        case 43:
            ret_val = 32;
            break;
        case 54:
            ret_val = 43;
            break;
        case 1034:
            ret_val = 43;
            break;
        case -1026:
            ret_val = 43;
            break;
        default:
            std::cout<<"WARNING: Id: "<<identif<<" not recognized!"<<std::endl;
            ret_val = -1;
            break;			
    }
    return ret_val;
}


//map of the possible right neighbors (not used up to now)
std::vector<int> CAGeneratorFromPairs::IsAtRight(int identif){
    std::vector<int> ret_val;
    //as for lefties, but right neighbors can be more than one
    switch (identif){
        case 21:
            ret_val.push_back(32);
            break;
        case 32:
            ret_val.push_back(43);
            break;
        case 43:
            ret_val.push_back(54);
            ret_val.push_back(1034);
            ret_val.push_back(-1026);
            break;
        case 54:
            ret_val.push_back(0);
            break;
        case 1034:
            ret_val.push_back(0);
            break;
        case -1026:
            ret_val.push_back(0);
            break;
        default:
            std::cout<<"WARNING: Id: "<<identif<<" not recognized!"<<std::endl;
            ret_val.push_back(-1);
            break;			
    }
    
    return ret_val;	
}
