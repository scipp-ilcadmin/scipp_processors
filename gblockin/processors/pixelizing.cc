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
static TH1D* _bzPos;
static TH2D* _bxyPos;
static TH2D* _byzPos;
static TH2D* _byPos;
static TH2D* _zandr;
static TH3D* _l1bxyzPos;
static TH3D* _l2bxyzPos;
static TH3D* _l3bxyzPos;
static TH3D* _l4bxyzPos;
static TH3D* _l5bxyzPos;
static TH1I* _blayers;
static TH1D* _l1radVals;
static TH1D* _l2radVals;
static TH1D* _l3radVals;
static TH1D* _l4radVals;
static TH1D* _l5radVals;
static TH1D* _l1thetas;
static TH1D* _l2thetas;
static TH1D* _l3thetas;
static TH1D* _l4thetas;
static TH1D* _l5thetas;
static TH1D* _zpix;
static int _nEvt = 0;

static vector<double> bposxVals;
static vector<double> bposyVals;
static vector<double> bposzVals;
static vector<double> bposVals;
static vector<pair<double, int>> data;
static vector<int> blayers;
static vector<double> xyradius;
static vector<double> l1vals;
static vector<double> l2vals;
static vector<double> l3vals;
static vector<double> l4vals;
static vector<double> l5vals;
static vector<double> thetavals;

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
  _bxyPos = new TH2D("bxypos", "Barrel (x,y) distribution ", 500, -65, 65, 500, -65, 65);
  _l1bxyzPos = new TH3D("l1bxyzpos", "l1bxyzpos", 124, -62.0, 62.0, 124, -62.0, 62.0, 124, -65.0, 65.0);
  _l2bxyzPos = new TH3D("l2bxyzpos", "l2bxyzpos", 124, -62.0, 62.0, 124, -62.0, 62.0, 124, -65.0, 65.0);
  _l3bxyzPos = new TH3D("l3bxyzpos", "l3bxyzpos", 124, -62.0, 62.0, 124, -62.0, 62.0, 124, -65.0, 65.0);
  _l4bxyzPos = new TH3D("l4bxyzpos", "l4bxyzpos", 124, -62.0, 62.0, 124, -62.0, 62.0, 124, -65.0, 65.0);
  _l5bxyzPos = new TH3D("l5bxyzpos", "l5bxyzpos", 124, -62.0, 62.0, 124, -62.0, 62.0, 124, -65.0, 65.0);
  _bpos = new TH1D("bpos", "bpos", 100, -100, 100);
  _l1radVals = new TH1D("l1radVals", "l1radVals", 100, 13, 17);
  _l2radVals = new TH1D("l2radVals", "l2radVals", 200, 21.5, 25);
  _l3radVals = new TH1D("l3radVals", "l3radVals", 200, 34.5, 37);
  _l4radVals = new TH1D("l4radVals", "l4radVals", 200, 47, 50);
  _l5radVals = new TH1D("l5radVals", "l5radVals", 200, 59.5, 62);
  _l1thetas = new TH1D("l1thetas", "l1thetas", 50, -20, 380);
  _l2thetas = new TH1D("l2thetas", "l2thetas", 50, -20, 380);
  _l3thetas = new TH1D("l3thetas", "l3thetas", 50, -20, 380);
  _l4thetas = new TH1D("l4thetas", "l4thetas", 50, -20, 380);
  _l5thetas = new TH1D("l5thetas", "l5thetas", 50, -20, 380);
  _zandr = new TH2D("zandr", "zandr", 131, -70, 70, 131, 0, 20);
  _zpix = new TH1D("zpix", "zpix", 131, -65, 65);
  _nEvt = 0;
}

void pixelizing::processRunHeader( LCRunHeader* run)
{

}

void pixelizing::processEvent( LCEvent * evt)
{
  _nEvt++;
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  int size = barrelHits->getNumberOfElements();
  _nEvt++;
  for (int i = 0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( barrelHits );
      double bposx = hit->getPosition()[0]; // indecies for x,y,z components;
      double bposy = hit->getPosition()[1];
      double bposz = hit->getPosition()[2];
      double pos = sqrt( (bposx*bposx) + (bposy*bposy) + (bposz*bposz) );
      double xyrad = sqrt( (bposx*bposx) + (bposy*bposy) );
      int layer = idDec(hit) [ILDCellID0::layer];
      Entry entry = std::make_pair(xyrad, layer);
      double theta = (atan2(bposy, bposx) + M_PI); //angle in radians ranging from 0->2Pi
      double arc = xyrad * theta;
      thetavals.push_back(theta);
      bposxVals.push_back(bposx);
      bposyVals.push_back(bposy);
      bposzVals.push_back(bposz);
      blayers.push_back(layer);
      _bxyPos->Fill(bposx, bposy);
      switch (entry.second)
	{
	case 1:
	  _l1radVals->Fill(entry.first);
          l1vals.push_back(entry.first);
          _l1thetas->Fill(theta);
	  _l1bxyzPos->Fill(bposx, bposy, bposz);
	  _zandr->Fill(bposz, xyrad);
	  _l1bxyzPos->Fill(bposx,bposy,bposz);
	  //	  _bxyPos->Fill(bposx, bposy);
	  break;
	case 2:
	  _l2radVals->Fill(entry.first);
	  l2vals.push_back(entry.first);
	  _l2thetas->Fill(theta);
	  _l2bxyzPos->Fill(bposx, bposy, bposz);
	  break;
	case 3:
	  _l3radVals->Fill(entry.first);
	  l3vals.push_back(entry.first);
	  _l3thetas->Fill(theta);
	  _l3bxyzPos->Fill(bposx, bposy, bposz);
	  break;
	case 4:
	  _l4radVals->Fill(entry.first);
	  l4vals.push_back(entry.first);
	  _l4thetas->Fill(theta);
	  _l4bxyzPos->Fill(bposx, bposy, bposz);
	  break;
	case 5:
	  _l5radVals->Fill(entry.first);
	  l5vals.push_back(entry.first);
	  _l5thetas->Fill(theta);
	  _l5bxyzPos->Fill(bposx, bposy, bposz);
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

void pixelizing::check( LCEvent * evt)
{

}

void pixelizing::end()
{
  //cout << "MAX bPosX: " << getMax(bposxVals) << " MIN bPosX: " << getMin(bposxVals) <<  endl;
  //cout << "MAX bPosY: " << getMax(bposyVals) << " MIN bPosY: " << getMin(bposyVals) <<  endl;
  //cout << "MAX bPosZ: " << getMax(bposzVals) << " MIN bPosZ: " << getMin(bposzVals) <<  endl;
  //cout << _nEvt << endl;
  _bxyPos->GetXaxis()->SetTitle("X (mm)");
  _bxyPos->GetYaxis()->SetTitle("Y (mm)");
  //_bxyzPos->GetZaxis()->SetTitle("Z (mm)");


  _l1radVals->SetFillColor(kRed);
  _l2radVals->SetFillColor(kBlack);
  _l3radVals->SetFillColor(kAzure);
  _l4radVals->SetFillColor(kTeal);
  _l5radVals->SetFillColor(kYellow);
  
  _l1thetas->SetFillColor(kRed);
  _l2thetas->SetFillColor(kBlack);
  _l3thetas->SetFillColor(kAzure);
  _l4thetas->SetFillColor(kTeal);
  _l5thetas->SetFillColor(kYellow);

  _rootfile->Write();

  //cout << "MAX theta: " << getMax(thetavals) << " MIN theta: " << getMin(thetavals) <<  endl;
  //cout << "MAX l2val: " << getMax(l2vals) << " MIN l2val: " << getMin(l2vals) <<  endl;
  //cout << "MAX l3val: " << getMax(l3vals) << " MIN l3val: " << getMin(l3vals) <<  endl;
  //cout << "MAX l4val: " << getMax(l4vals) << " MIN l4val: " << getMin(l4vals) <<  endl;
  //cout << "MAX l5val: " << getMax(l5vals) << " MIN l5val: " << getMin(l5vals) <<  endl;
}
