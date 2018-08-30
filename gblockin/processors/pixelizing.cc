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
static TH2D* _bxzPos;
static TH2D* _byzPos;
static TH3D* _bxyzPos;
static TH1I* _blayers;
static TH1D* _radVals;
static int _nEvt = 0;

static vector<double> bposxVals;
static vector<double> bposyVals;
static vector<double> bposzVals;
static vector<double> bposVals;
static vector<pair<double, int>> data;
static vector<int> blayers;
static vector<double> xyradius;

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
  _bzPos = new TH1D("bzposhits", "bzposhits", 57, -100, 100);
  _bxyzPos = new TH3D("bxyzpos", "bxyzpos", 200, -62.0, 62.0, 200, -62.0, 62.0, 200, -65.0, 65.0);
  _bpos = new TH1D("bpos", "bpos", 100, -100, 100);
  _radVals = new TH1D("radVals", "radVals", 200, 12, 27);
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
      double bposy = hit->getPosition()[1];
      double bposz = hit->getPosition()[2];
      double pos = sqrt( (bposx*bposx) + (bposy*bposy) + (bposz*bposz) );
      double xyrad = sqrt( (bposx*bposx) + (bposy*bposy) );
      int layer = idDec( hit) [ILDCellID0::layer];
      Entry entry = std::make_pair(xyrad, layer);

      bposxVals.push_back(bposx);
      bposyVals.push_back(bposy);
      bposzVals.push_back(bposz);
      xyradius.push_back(xyrad);
      blayers.push_back(layer);
      data.push_back(entry);
      
    
      _bzPos->Fill(bposz);
      _bpos->Fill(pos);

      _bxyzPos->Fill(bposx, bposy, bposz);
      _radVals->Fill(xyrad);
/*      std::sort(
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
  /*cout << "MAX bPosX: " << getMax(bposxVals) << " MIN bPosX: " << getMin(bposxVals) <<  endl;
  cout << "MAX bPosY: " << getMax(bposyVals) << " MIN bPosY: " << getMin(bposyVals) <<  endl;
  cout << "MAX bPosZ: " << getMax(bposzVals) << " MIN bPosZ: " << getMin(bposzVals) <<  endl;
  cout << _nEvt << endl;
  _bxyzPos->GetXaxis()->SetTitle("X (mm)");
  _bxyzPos->GetYaxis()->SetTitle("Y (mm)");
  _bxyzPos->GetZaxis()->SetTitle("Z (mm)");
  */
  _rootfile->Write();
  _rootfile->Close();
}
