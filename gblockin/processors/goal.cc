#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * goal.cc
 * @author Gregory Blockinger
 * September 26th, 2018
 *
 */

#include "goal.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>
#include <string>

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
using Entry = std::pair<double, int>;

goal goal;

static TFile* _rootfile;

static TH1D* _bzPos;
static TH2D* _bxyPos;
static TH2D* _byzPos;
static TH2D* _byPos;
static TH3D* _bxyzPos;
static TH2D* _zandr;
static TH3D* _l1bxyzPos;
static TH1I* _blayers;
static TH1D* _l1radVals;
static TH1D* _l1thetas;
static TH1D* _zpix;
static int _nEvt = 0;

static vector<double> bposxVals;
static vector<double> bposyVals;
static vector<double> bposzVals;
static vector<pair<double, int>> data;
static vector<int> blayers;
static vector<double> xyradius;
static vector<double> l1vals;
static vector<double> thetavals;
static vector<double> modvals;
static vector<double> sidevals;
static vector<double> sensorvals;

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

goal::goal() : Processor("goal") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));
}

void goal::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("goal.root", "RECREATE");
  _bxyPos = new TH2D("bxypos", "Barrel (X,Y) Distribution", 500, -20, 20, 500, -20, 20);
  _bxyzPos = new TH3D("bxyzpos", "Barrel (X,Y,Z) Distribution", 124, -80.0, 80.0, 124, -80.0, 80.0, 124, -190.0, 190.0);
  _l1radVals = new TH1D("l1radVals", "EndCap Radial (X,Y) Distribution (All layers)", 100, -100, 100);
  _l1thetas = new TH1D("l1thetas", "l1thetas", 50, -20, 380);
  _zandr = new TH2D("zandr", "zandr", 131, -70, 70, 131, -20, 20);
  _zpix = new TH1D("zpix", "zpix", 131, -10, 370);
  _nEvt = 0;
 }

void goal::processRunHeader( LCRunHeader* run)
{

}

void goal::processEvent( LCEvent * evt)
{
  _nEvt++;
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  int size = barrelHits->getNumberOfElements();
  _nEvt++;
  for (int i = 0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( barrelHits );
      int layer = idDec(hit) [ILDCellID0::layer];
      if (layer == 1)
	{
	  double bposx = hit->getPosition()[0]; // indecies for x,y,z components;
	  double bposy = hit->getPosition()[1];
	  double bposz = hit->getPosition()[2];
	  double xyrad = sqrt( (bposx*bposx) + (bposy*bposy) );
	  int side = idDec(hit) [ILDCellID0::side];
	  int module = idDec(hit) [ILDCellID0::module];
	  int sensor = idDec(hit) [ILDCellID0::sensor];
	  Entry entry = std::make_pair(xyrad, layer);
	  double theta = (atan2(bposy, bposx) + M_PI) * 180/M_PI; //angle in radians ranging from 0->2Pi
	  double arc = xyrad * theta;
	  thetavals.push_back(theta);
	  bposxVals.push_back(bposx);
	  bposyVals.push_back(bposy);
	  bposzVals.push_back(bposz);
	  blayers.push_back(layer);
	  sidevals.push_back(side);
	  sensorvals.push_back(sensor);
	  _l1radVals->Fill(xyrad);
	  _bxyPos->Fill(bposx, bposy);
	  _zpix->Fill(theta);
	}

    }
}

void goal::check( LCEvent * evt)
{

}

void goal::end()
{
  cout << "MAX bPosZ: " << getMax(bposzVals) << " MIN bPosZ: " << getMin(bposzVals) <<  endl;
  cout << "MAX bPosY: " << getMax(bposyVals) << " MIN bPosY: " << getMin(bposyVals) <<  endl;
  cout << "MAX bPosX: " << getMax(bposxVals) << " MIN bPosX: " << getMin(bposxVals) <<  endl;
  //cout << *max_element(bposzVals.begin(), bposzVals.end()) << endl;
  //cout << *min_element(bposzVals.begin(), bposzVals.end()) << endl;

  _bxyPos->GetXaxis()->SetTitle("X (mm)");
  _bxyPos->GetYaxis()->SetTitle("Y (mm)");
  _bxyzPos->GetZaxis()->SetTitle("Z (mm)");

  //_l1radVals->SetFillColor(kRed);
  //_l1thetas->SetFillColor(kRed);
  _rootfile->Write();

}
