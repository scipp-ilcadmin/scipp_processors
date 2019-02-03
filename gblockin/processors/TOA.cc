#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * TOA.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
 */

#include "TOA.h"
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

TOA TOA;


static TFile* _rootfile;

static TH1D* _bxPos;
static TH1D* _byPos;
static TH2D* _bxyPos;
static TH2D* _bxzPos;
static TH2D* _byzPos;
static TH3D* _bxyzPos;
static TH1I* _blayers;

static TH1D* _exPos;
static TH1D* _eyPos;
static TH2D* _exyPos;
static TH2D* _exzPos;
static TH2D* _eyzPos;
static TH3D* _exyzPos;

static TH3D* _dpos;
static int _nEvt = 0;

static vector<double> bposxVals;
static vector<double> bposyVals;
static vector<double> bposzVals;
static vector<int> blayers;

static vector<double> eposxVals;
static vector<double> eposyVals;
static vector<double> eposzVals;

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

TOA::TOA() : Processor("TOA") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));
}

void TOA::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("TOA.root", "RECREATE");
  //_blayers = new TH1I("blayers", "blayers", 7, -1, 7);
  //_bxPos = new TH1D("bxposhits", "bxposhits", 57, -67.0, 67.0);
  //_byPos = new TH1D("byposhits", "byposhits", 57, -67.0, 67.0);
  _bxyPos = new TH2D("bxypos", "bxypos", 200, -67.0, 67.0, 200, -67.0, 67.0);
  //_bxzPos = new TH2D("bxzpos", "bxzpos", 5, -67.0, 67.0, 5, -65.0, 65.0);
  //_byzPos = new TH2D("byzpos", "byzpos", 100, -67.0, 67.0, 100, -65.0, 65.0);
  _bxyzPos = new TH3D("bxyzpos", "bxyzpos", 200, -62.0, 62.0, 200, -62.0, 62.0, 200, -65.0, 65.0);
  //_exPos = new TH1D("exposhits", "exposhits", 57, -62.0, 55.0);
  //_eyPos = new TH1D("eyposhits", "eyposhits", 57, -41.0, 48.0);
  _exyPos = new TH2D("exypos", "exypos", 200, -80.0, 80.0, 200, -80.0, 80.0);
  //_exzPos = new TH2D("exzpos", "exzpos", 100, -80.0, 80.0, 100, -185.0, 185.0);
  //_eyzPos = new TH2D("eyzpos", "eyzpos", 100, -80.0, 80.0, 100, -185.0, 185.0);
  _exyzPos = new TH3D("exyzpos", "exyzpos", 200, -80.0, 80.0, 200, -80.0, 80.0, 200, -185.0, 185.0);
  _dpos = new TH3D("_dpos", "_dpos", 200, -80, 80, 200, -80, 80, 200, -200, 200);
  _nEvt = 0;
}

void TOA::processRunHeader( LCRunHeader* run)
{

}

void TOA::processEvent( LCEvent * evt)
{
  _nEvt++;
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHits"); 

  _nEvt++;
  for (int i = 0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( barrelHits );

      double bposx = hit->getPosition()[0]; // indecies for x,y,z components;
      bposxVals.push_back(bposx);
      double bposy = hit->getPosition()[1];
      bposyVals.push_back(bposy);
      double bposz = hit->getPosition()[2];
      bposzVals.push_back(bposz);
      int layer = idDec( hit) [ILDCellID0::layer];
      blayers.push_back(layer);
      //_blayers->Fill(layer);
      //_bxPos->Fill(bposx);
      //_byPos->Fill(bposy);
      _bxyPos->Fill(bposx, bposy);
      _bxyzPos->Fill(bposx, bposy, bposz);
      _dpos->Fill(bposx, bposy, bposz);
      //_bxzPos->Fill(bposx, bposz);
      //_byzPos->Fill(bposy, bposz);
    }

  for (int i = 0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( endcapHits );
      
      double eposx = hit->getPosition()[0]; // indecies for x,y,z components;
      eposxVals.push_back(eposx);
      double eposy = hit->getPosition()[1];
      eposyVals.push_back(eposy);
      double eposz = hit->getPosition()[2];
      eposzVals.push_back(eposz);
      //_exPos->Fill(eposx);
      //_eyPos->Fill(eposy);
      _exyPos->Fill(eposx, eposy);
      _exyzPos->Fill(eposx, eposy, eposz);
      _dpos->Fill(eposx, eposy, eposz);
      //_exzPos->Fill(eposx, eposz);
      //_eyzPos->Fill(eposy, eposz);

    }
}

void TOA::check( LCEvent * evt)
{

}

void TOA::end()
{
  cout << "MAX bPosX: " << getMax(bposxVals) << " MIN bPosX: " << getMin(bposxVals) <<  endl;
  cout << "MAX bPosY: " << getMax(bposyVals) << " MIN bPosY: " << getMin(bposyVals) <<  endl;
  //cout << "MAX bPosZ: " << getMax(bposzVals) << " MIN bPosZ: " << getMin(bposzVals) <<  endl;
  //cout << "MAX ePosX: " << getMax(eposxVals) << " MIN ePosX: " << getMin(eposxVals) <<  endl;
  //cout << "MAX ePosY: " << getMax(eposyVals) << " MIN ePosY: " << getMin(eposyVals) <<  endl;
  //cout << "MAX ePosZ: " << getMax(eposzVals) << " MIN ePosZ: " << getMin(eposzVals) <<  endl;
  cout << _nEvt << endl;
  //_bxyzPos->GetXaxis()->SetTitle("X (mm)");
  //_bxyzPos->GetYaxis()->SetTitle("Y (mm)");
  //_bxyzPos->GetZaxis()->SetTitle("Z (mm)");
  //_exyzPos->GetXaxis()->SetTitle("X (mm)");
  //_exyzPos->GetYaxis()->SetTitle("Y (mm)");
  //_exyzPos->GetZaxis()->SetTitle("Z (mm)");

  _rootfile->Write();
  _rootfile->Close();
}
