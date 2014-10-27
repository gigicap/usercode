#pragma once

// system include files
#include <memory>
#include <vector>
#include <map>
#include <deque>
#include <queue>
#include <iostream>
#include <set>
#include <limits>
#include <array>
#include <algorithm>

#include <math.h>

// cmssw includes
#include "DataFormats/GeometryVector/interface/GlobalVector.h"  //<--
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"   //<--
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"        //<--
#include "RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h"        //<--
#include "DataFormats/TrackerRecHit2D/interface/OmniClusterRef.h"

//root
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"
#include "TVector3.h"

#include <unordered_map>


//Class CAcell definition
class CAcell{
    
public:
	//CAcell hits (from seeds)
   std::vector<const BaseTrackerRecHit*> hits;

	//fit parameters
	// LocalTrajectoryParameters FitHit0Local;
	// LocalTrajectoryError FitHit0LocalError;
    
	// LocalTrajectoryParameters FitHit1Local;
	// LocalTrajectoryError FitHit1LocalError;
    
	// LocalTrajectoryParameters FitHit2Local;
	// LocalTrajectoryError FitHit2LocalError;
    
	bool FitSuccessful;
	//double FitChiSquared;
	//int FitNdof;
	//int FitFoundHits;
	double LocalMomentum;

	int eventNumber;
    
    //triplet identifier e.g. (1,2,3) =  123
    int tripletIdentifier;

	int hitId0;
	int hitId1;
	int hitId2;

    GlobalPoint p0 ,p1, p2;
    
    double dEta;
    
    unsigned int rawId0;
    unsigned int rawId1;
    unsigned int rawId2;
    
    //neighborMAP
    std::vector<CAcell *> left_neighbors;       //list of the cells sharing hits 0-1 with the cell
    std::vector<CAcell *> right_neighbors;      //list of the cells sharing hits 1-2 with the cell

    //How many neighbors does the cell have?
    int IsNeighborly;
    
    
    //just for convenience and debug purposes
    int TripletNumber;

    //for the backwardCA
    int IsUsed;
    
    //for fitting
    PTrajectoryStateOnDet seedStartState;
    
    CAcell(){
        TripletNumber = 0;
        CAstatus = 0;
        tripletIdentifier = 0;
        }

    
//double getFitChiSquared(){ return FitChiSquared;}
//int getFitNdof(){return FitNdof;}
//int getFitFoundHits(){return FitFoundHits;}
double getLocalMomentum(){ return LocalMomentum;}

/*double geth0qbp(){ return FitHit0Local.qbp();}
double geth0dxdz(){ return FitHit0Local.dxdz();}
double geth0dydz(){ return FitHit0Local.dydz();}

double geth1qbp(){ return FitHit0Local.qbp();}
double geth1dxdz(){ return FitHit0Local.dxdz();}
double geth1dydz(){ return FitHit0Local.dydz();}

double geth2qbp(){ return FitHit2Local.qbp();}
double geth2dxdz(){ return FitHit2Local.dxdz();}
double geth2dydz(){ return FitHit2Local.dydz();}*/

int geteventNumber(){return eventNumber;}
    
int gethitId0(){return hitId0;}
int gethitId1(){return hitId1;}
int gethitId2(){return hitId2;}

    
//CAstuff
    int CAstatus;
        
};
