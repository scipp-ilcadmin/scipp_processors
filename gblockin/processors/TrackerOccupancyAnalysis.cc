#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * TrackerOccupancyAnalysis.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
 *
 *
 */

#include "TrackerOccupancyAnalysis.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <algorithm>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;

TrackerOccupancyAnalysis TrackerOccupancyAnalysis;

static int _nEvt = 0;
static TFile* _rootfile;
static vector<float> barrelMomentumVals;
static vector<double> barrelPosVals;
static vector<int> barrelCell0Vals;
static vector<int> barrelCell1Vals;
static vector<float> barrelEdepVals;
static vector<float> endcapMomentumVals;
static vector<double> endcapPosVals;
static vector<int> endcapCell0Vals;
static vector<int> endcapCell1Vals;
static vector<float> endcapEdepVals;


TrackerOccupancyAnalysis::TrackerOccupancyAnalysis() : Processor("TrackerOccupancyAnalysis") 
{
  _description = "Protype Processor";
  //    registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void TrackerOccupancyAnalysis::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  //  _rootfile = new TFile("TrackerOccupancyAnalysis.root", "RECREATE");
  _nEvt = 0;

}

void TrackerOccupancyAnalysis::processRunHeader( LCRunHeader* run)
{

}

void TrackerOccupancyAnalysis::processEvent( LCEvent * evt)
{
  _nEvt++;
  int check_layer = 0;
  int hit_count_limit = 100;
  bool use_limit = false;
  bool reject_negative = false;
  int hit_count = 0;
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHits");

  for (int i =0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit=dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      float barrelMomentum = *hit->getMomentum();
      barrelMomentumVals.push_back(barrelMomentum);
      float barrelPos = *hit->getPosition();
      barrelPosVals.push_back(barrelPos);
      int barrelCell0 = hit->getCellID0();
      barrelCell0Vals.push_back(barrelCell0);
      int barrelCell1 = hit->getCellID1();
      barrelCell1Vals.push_back(barrelCell1);
      float barrelEdep = hit->getEDep();
      barrelEdepVals.push_back(barrelEdep);
    }
  for (int i=0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit2=dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      float endcapMomentum = *hit2->getMomentum();
      endcapMomentumVals.push_back(endcapMomentum);
      float endcapPos = *hit2->getPosition();
      endcapPosVals.push_back(endcapPos);
      int endcapCell0 = hit2->getCellID0();
      endcapCell0Vals.push_back(endcapCell0);
      int endcapCell1 = hit2->getCellID1();
      float endcapEdep = hit2->getEDep();
      endcapEdepVals.push_back(endcapEdep);
    }
  /* cout << " barrel_max_MOM::  " << *max_element(barrelMomentumVals.begin(), barrelMomentumVals.end()) << endl;
  cout << " barrel_min_MOM::  " << *min_element(barrelMomentumVals.begin(), barrelMomentumVals.end()) << endl;
  cout << " barrel_max_POS::  " << *max_element(barrelPosVals.begin(), barrelPosVals.end()) << endl;
  cout << " barrel_min_POS::  " << *min_element(barrelPosVals.begin(), barrelPosVals.end()) << endl;
  cout << " barrel_max_EDEP:: " << *max_element(barrelEdepVals.begin(), barrelEdepVals.end()) << endl;
  cout << " barrel_min_EDEP:: " << *min_element(barrelEdepVals.begin(), barrelEdepVals.end()) << endl << endl;
  cout << " endcap_max_MOM::  " << *max_element(endcapMomentumVals.begin(), endcapMomentumVals.end()) << endl;
  cout << " endcap_min_MOM::  " << *min_element(endcapMomentumVals.begin(), endcapMomentumVals.end()) << endl;
  cout << " endcap_max_POS::  " << *max_element(endcapPosVals.begin(), endcapPosVals.end()) << endl;
  cout << " endcap_min_POS::  " << *min_element(endcapPosVals.begin(), endcapPosVals.end()) << endl;
  cout << " endcap_max_EDEP:: " << *max_element(endcapEdepVals.begin(), endcapEdepVals.end()) << endl;
  cout << " endcap_min_EDEP:: " << *min_element(endcapEdepVals.begin(), endcapEdepVals.end()) << endl << endl;
  */
}


void TrackerOccupancyAnalysis::check( LCEvent * evt)
{

}

void TrackerOccupancyAnalysis::end()
{
  
}
