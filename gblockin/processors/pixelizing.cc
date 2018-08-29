#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * pixelizing.cc
 * @author Gregory Blockinger
 * August 26th, 2018
 *
 */

#include "pixelizing.h"
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

pixelizing pixelizing;


static TFile* _rootfile;

static TH1D* _bpos;
static TH1D* _bxPos;
static TH1D* _byPos;
static TH2D* _bxyPos;
static TH2D* _bxzPos;
static TH2D* _byzPos;
static TH3D* _bxyzPos;
static TH1I* _blayers;

static TH1D* _epos;
static TH1D* _exPos;
static TH1D* _eyPos;
static TH2D* _exyPos;
static TH2D* _exzPos;
static TH2D* _eyzPos;
static TH3D* _exyzPos;

static int _nEvt = 0;

static vector<double> bposxVals;
static vector<double> bposyVals;
static vector<double> bposzVals;
static vector<double> bposVals;
static vector<pair<double, int>> data;
static vector<int> blayers;

static vector<double> eposxVals;
static vector<double> eposyVals;
static vector<double> eposzVals;
static vector<double> eposVals;


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

pixelizing::pixelizing() : Processor("pixelizing") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));
}

void pixelizing::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("pixelizing.root", "RECREATE");
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
  _bpos = new TH1D("bpos", "bpos", 100, -100, 100);
  _epos = new TH1D("epos", "epos", 100, -100, 100); 
  _nEvt = 0;
}

void pixelizing::processRunHeader( LCRunHeader* run)
{

}

void pixelizing::processEvent( LCEvent * evt)
{
  _nEvt++;
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHits"); 
  int size = barrelHits->getNumberOfElements();
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
      double pos = sqrt( (bposx*bposx) + (bposy*bposy) + (bposz*bposz) );
      _bpos->Fill(pos); 
      int layer = idDec( hit) [ILDCellID0::layer];
      blayers.push_back(layer);
      //_blayers->Fill(layer);
      //_bxPos->Fill(bposx);
      //_byPos->Fill(bposy);
      _bxyPos->Fill(bposx, bposy);
      _bxyzPos->Fill(bposx, bposy, bposz);
      //_bxzPos->Fill(bposx, bposz);
      //_byzPos->Fill(bposy, bposz);
      Entry entry = std::make_pair(pos, layer);
      data.push_back(entry);
      std::sort(
		data.begin(), 
		data.end(), 
		[&](const Entry &a, const Entry &b)
		{ 
		  return a.first < b.first;
		});
        
      std::cout << "Sorted (by layer)" << std::endl;
      std::for_each(
		    data.begin(), 
		    data.end(),
		    [](const Entry &entry)
		    {
		      std::cout << entry.first << " " << entry.second << std::endl;
		    });

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
      //_exzPos->Fill(eposx, eposz);
      //_eyzPos->Fill(eposy, eposz);

    }
}

void pixelizing::check( LCEvent * evt)
{

}

void pixelizing::end()
{
  /*cout << "MAX bPosX: " << getMax(bposxVals) << " MIN bPosX: " << getMin(bposxVals) <<  endl;
  cout << "MAX bPosY: " << getMax(bposyVals) << " MIN bPosY: " << getMin(bposyVals) <<  endl;
  cout << "MAX bPosZ: " << getMax(bposzVals) << " MIN bPosZ: " << getMin(bposzVals) <<  endl;
  cout << "MAX ePosX: " << getMax(eposxVals) << " MIN ePosX: " << getMin(eposxVals) <<  endl;
  cout << "MAX ePosY: " << getMax(eposyVals) << " MIN ePosY: " << getMin(eposyVals) <<  endl;
  cout << "MAX ePosZ: " << getMax(eposzVals) << " MIN ePosZ: " << getMin(eposzVals) <<  endl;
  cout << _nEvt << endl;
  _bxyzPos->GetXaxis()->SetTitle("X (mm)");
  _bxyzPos->GetYaxis()->SetTitle("Y (mm)");
  _bxyzPos->GetZaxis()->SetTitle("Z (mm)");
  _exyzPos->GetXaxis()->SetTitle("X (mm)");
  _exyzPos->GetYaxis()->SetTitle("Y (mm)");
  _exyzPos->GetZaxis()->SetTitle("Z (mm)");
  */
  _rootfile->Write();
  _rootfile->Close();
}
