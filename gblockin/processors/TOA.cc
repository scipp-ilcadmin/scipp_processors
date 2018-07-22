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
static TH3D* _bxyzPos;

static TH1D* _exPos;
static TH1D* _eyPos;
static TH2D* _exyPos;
static TH3D* _exyzPos;

static int _nEvt = 0;

static vector<double> bposxVals;
static vector<double> bposyVals;
static vector<double> bposzVals;

static vector<double> eposxVals;
static vector<double> eposyVals;
static vector<double> eposzVals;

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

  _bxPos = new TH1D("bxposhits", "bxposhits", 57, -62.0, 55.0);
  _byPos = new TH1D("byposhits", "byposhits", 57, -41.0, 48.0);
  _bxyPos = new TH2D("bxypos", "bxypos", 57, -100.0, 100.0, 57, -100.0, 100.0);
  _bxyzPos = new TH3D("bxyzpos", "bxyzpos", 57, -100.0, 100.0, 57, -100.0, 100.0, 57, -100.0, 100.0);
  _exPos = new TH1D("exposhits", "exposhits", 57, -62.0, 55.0);
  _eyPos = new TH1D("eyposhits", "eyposhits", 57, -41.0, 48.0);
  _exyPos = new TH2D("exypos", "exypos", 57, -100.0, 100.0, 57, -100.0, 100.0);
  _exyzPos = new TH3D("exyzpos", "exyzpos", 57, -100.0, 100.0, 57, -100.0, 100.0, 57, -100.0, 100.0);
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
      _bxPos->Fill(bposx);
      _byPos->Fill(bposy);
      _bxyPos->Fill(bposx, bposy);
      _bxyzPos->Fill(bposx, bposy, bposz);
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
      _exPos->Fill(eposx);
      _eyPos->Fill(eposy);
      _exyPos->Fill(eposx, eposy);
      _exyzPos->Fill(eposx, eposy, eposz);

    }
}

void TOA::check( LCEvent * evt)
{

}

void TOA::end()
{
  _rootfile->Write();
}