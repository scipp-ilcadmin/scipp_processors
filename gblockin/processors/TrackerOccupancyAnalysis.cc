#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * TrackerOccupancyAnalysis.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
 */

#include "TrackerOccupancyAnalysis.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <algorithm>

#include <UTIL/ILDConf.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;

TrackerOccupancyAnalysis TrackerOccupancyAnalysis;


static TFile* _rootfile;
static TH2D* _xyPos;
static TH3D* _xyzPos;
static int _nEvt = 0;

static vector<int> layers;
static vector<int> subDets;
static vector<int> modules;
static vector<int> sensors;
static vector<int> sides;
static vector<float> eDepVals;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static vector<int> cell0Vals;
static vector<int> cell1Vals;
static vector<int> nmccontsVals;
static vector<int> tileIDVals;
static vector<float> xmomVals;
static vector<float> ymomVals;
static vector<float> zmomVals;
template<typename T>
static T getMax(vector<T> &vec) 
{
  return *max_element(vec.begin(), vec.end());
}
template<typename T>
static T getMin(vector<T> &vec)
{
  return *min_element(vec.begin(), vec.end());
}

TrackerOccupancyAnalysis::TrackerOccupancyAnalysis() : Processor("TrackerOccupancyAnalysis") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void TrackerOccupancyAnalysis::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("TOA.root", "RECREATE");
  _xyPos = new TH2D("xypos", "xypos", 57, -100.0, 100.0, 57, -100.0, 100.0);
  _xyzPos = new TH3D("xyzpos", "xyzpos", 57, -100.0, 100.0, 57, -100.0, 100.0, 57, -185.0, 185.0);
  _nEvt = 0;

}

void TrackerOccupancyAnalysis::processRunHeader( LCRunHeader* run)
{

}

void TrackerOccupancyAnalysis::processEvent( LCEvent * evt)
{
  _nEvt++;
  //LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHits");
  for (int i = 0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( endcapHits );
      int layer = idDec( hit )[ILDCellID0::layer];
      layers.push_back(layer);
      int subdet = idDec( hit )[ILDCellID0::subdet];
      subDets.push_back(subdet);
      int module = idDec( hit )[ILDCellID0::module];
      modules.push_back(module);
      int sensor = idDec( hit )[ILDCellID0::sensor];
      sensors.push_back(sensor);
      int side = idDec( hit )[ILDCellID0::side];
      sides.push_back(side);
      double posx = hit->getPosition()[0]; // indecies for x,y,z components;
      posxVals.push_back(posx);
      double posy = hit->getPosition()[1];
      posyVals.push_back(posy);
      double posz = hit->getPosition()[2];
      poszVals.push_back(posz);
      int cell0 = hit->getCellID0();
      cell0Vals.push_back(cell0);
      int cell1 = hit->getCellID1();
      cell1Vals.push_back(cell1);
      float eDep = hit->getEDep();
      eDepVals.push_back(eDep);
      float xmom = hit->getMomentum()[0];
      xmomVals.push_back(xmom);
      float ymom = hit->getMomentum()[1];
      ymomVals.push_back(ymom);
      float zmom = hit->getMomentum()[2];
      zmomVals.push_back(zmom);
      if (side == 1) 
	{
	  //cout << module << endl;
	  _xyPos->Fill(posx, posy);
	  _xyzPos->Fill(posx, posy, posz);
	}
    }
  }


void TrackerOccupancyAnalysis::check( LCEvent * evt)
{

}

void TrackerOccupancyAnalysis::end()
{
  //cout << " max energy: " << getMax(eDepVals) << " MinEnergy: " << getMin(eDepVals) << endl;
  //cout << " max xmom: " << getMax(xmomVals) << " min xmom: " << getMin(xmomVals) << endl;
  //cout << " max ymom: " << getMax(ymomVals) << " min ymom: " <<getMin(ymomVals) << endl;
  //cout << " max zmom: " << getMax(xmomVals) << " min zmom: " <<getMin(zmomVals) << endl;
  _rootfile->Write();
}
