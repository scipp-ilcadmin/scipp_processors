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
  _bxyPos = new TH2D("bxypos", "bxypos", 500, 12, 17, 500, -2, 11);
  _bxyzPos = new TH3D("bxyzpos", "Barrel (X,Y,Z) Distribution", 124, -62.0, 62.0, 124, -62.0, 62.0, 124, -65.0, 65.0);
  _l1radVals = new TH1D("l1radVals", "l1radVals", 100, 13, 17);
  _l1thetas = new TH1D("l1thetas", "l1thetas", 50, -20, 380);
  _zandr = new TH2D("zandr", "zandr", 131, -70, 70, 131, -20, 20);
  _zpix = new TH1D("zpix", "zpix", 131, -65, 65);
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
  int pixels[20];
  for (int i = 0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( barrelHits );
      double bposx = hit->getPosition()[0]; // indecies for x,y,z components;
      double bposy = hit->getPosition()[1];
      double bposz = hit->getPosition()[2];
      int mcpart = hit->getMCParticle()->getPDG();
      //cout << mcpart << endl;
      double xyrad = sqrt( (bposx*bposx) + (bposy*bposy) );
      int layer = idDec(hit) [ILDCellID0::layer];
      int side = idDec(hit) [ILDCellID0::side];
      int module = idDec(hit) [ILDCellID0::module];
      int sensor = idDec(hit) [ILDCellID0::sensor];
      Entry entry = std::make_pair(xyrad, layer);
      double theta = (atan2(bposy, bposx) + M_PI); //angle in radians ranging from 0->2Pi
      double arc = xyrad * theta;
      thetavals.push_back(theta);
      //bposxVals.push_back(bposx);
      //bposyVals.push_back(bposy);
      //bposzVals.push_back(bposz);
      blayers.push_back(layer);
      sidevals.push_back(side);
      sensorvals.push_back(sensor);
      _bxyzPos->Fill(bposx,bposy,bposz);
      switch (entry.second)
	{
	case 1:
	  {
	    vector<double> vecx;
	    vector<double> vecy;
	    _l1radVals->Fill(entry.first);
            l1vals.push_back(entry.first);
            _l1thetas->Fill(theta);
	    //_l1bxyzPos->Fill(bposx, bposy, bposz);
      	  }
	  break;
	case 2:
	  break;
	case 3:
	  break;
	case 4:
	  break;
	case 5:
	  break;
	}
    
/*    data.push_back(entry);
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
*/
    }

}

void goal::check( LCEvent * evt)
{

}

void goal::end()
{
  //cout << "MAX module: " << getMax(modvals) << " MIN module: " << getMin(modvals) <<  endl;
  //cout << "MAX bPosY: " << getMax(bposyVals) << " MIN bPosY: " << getMin(bposyVals) <<  endl;
  //cout << "MAX bPosX: " << getMax(bposxVals) << " MIN bPosX: " << getMin(bposxVals) <<  endl;
  //_bxyzPos->GetXaxis()->SetTitle("X (mm)");
  //_bxyzPos->GetYaxis()->SetTitle("Y (mm)");
  _bxyzPos->GetZaxis()->SetTitle("Z (mm)");


  _l1radVals->SetFillColor(kRed);
  _l1thetas->SetFillColor(kRed);
  _rootfile->Write();

}
