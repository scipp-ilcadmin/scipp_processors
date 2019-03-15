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
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <UTIL/ILDConf.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph2D.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;

TrackerOccupancyAnalysis TrackerOccupancyAnalysis;

typedef vector<vector<int>> PixelGrid;
typedef vector<PixelGrid> Layers;
typedef vector<string> PixIDs;
typedef vector<PixIDs> layerpixIDs;

static TH2D* totes;
static TFile* _rootfile;
static int _nEvt = 0;

//static vector<int> layers;
static vector<string> pixid;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static vector<string> uniquepix;
static layerpixIDs layerpixids(4, vector<string>());
static layerpixIDs layeruniqueids(4, vector<string>());
static Layers layers(4, PixelGrid(160*100, vector<int>(160*100,0)));
static vector<TH2D*> graphs;
static vector<TH1D*> angles;
static int hitcount = 0;
static vector<double> radii;
static TH3D* threedim;

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
  threedim = new TH3D("threedim", "3-D model", 100, -100, 100, 100, -100, 100, 100, -300, 300);
  totes = new TH2D("totes", "totesmagotes", 100, -100, 100, 100, -100, 100);
  for (int i =0; i < 4; i++)
    {
      graphs.push_back(new TH2D(Form("layer%d ", i), "layers", 1000, -80, 80, 1000, -80, 80));
      angles.push_back(new TH1D(Form("angles%d", i), "angles", 100, -10, 370));
    }
  
  _nEvt = 0;

}

void TrackerOccupancyAnalysis::processRunHeader( LCRunHeader* run)
{

}

void TrackerOccupancyAnalysis::processEvent( LCEvent * evt)
{
  cout << "DO STUFF    " << _nEvt << endl;
  _nEvt++;
  //LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHits");
  static const double xmin = 80;
  static const double ymin = 80;
  static const int step = 1;
  for (int i = 0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( endcapHits );
      hitcount++;
      int layer = idDec( hit )[ILDCellID0::layer];
      //int side = idDec ( hit )[ILDCellID0::side];
      //int module = idDec (hit)[ILDCellID0::side];
      int posx = (hit->getPosition()[0] + xmin)*100;
      int posy = (hit->getPosition()[1] + ymin)*100;
      threedim->Fill(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
      totes->Fill(hit->getPosition()[0], hit->getPosition()[1]);
      
      [&] ()
	{
	  int index = layer-1;
	  double theta = (atan2(posy, posx) + M_PI) * 180/M_PI;
	  string id = to_string(posx/step) + to_string(posy/step);
	  //pixid.push_back(id);
	  radii.push_back(sqrt((hit->getPosition()[0]*hit->getPosition()[0])+(hit->getPosition()[1]*hit->getPosition()[1])));	  
	  if(std::find(layeruniqueids[index].begin(), layeruniqueids[index].end(), id) == layeruniqueids[index].end())
	    {
	      layeruniqueids[index].push_back(id);
	    }
	  layerpixids[index].push_back(id);
	  angles[index]->Fill(theta);
	  graphs[index]->Fill(hit->getPosition()[0], hit->getPosition()[1]);
	  if ((posx < 160000 && posy < 160000) && (posx >=0 && posy >= 0))
	    {
	      layers[index][posx/step][posy/step]++;
	    }
	    else
	      {
		cout << "posx: " << posx << ", posy: " << posy << "ERROR" << endl;
	      }

	}();
    }
}

void TrackerOccupancyAnalysis::check( LCEvent * evt)
{

}

void TrackerOccupancyAnalysis::end()
{

  cout << "analysis finished" << endl;
  //cout << " max rad: " << getMax(radii) << " min rad: " << getMin(radii) << endl;
  //cout << " max y: " << getMax(posyVals) << " min y: " << getMin(posyVals) << endl;
  //cout << " max x: " << getMax(posxVals) << " min x: " << getMin(posxVals) << endl;
  int matters = (2*M_PI*73*73) - (2*M_PI*15*15);
  for (int i = 0; i < 4; i++)
    {
      cout << "total hits in layer " << i+1 << ": " << layerpixids[i].size() << endl;
      cout << "unique pixels hit in layer " << i+1 << ": " << layeruniqueids[i].size() << endl;
      double hitperc = static_cast<double>(layeruniqueids[i].size()) / static_cast<double>(matters) * 100;
      cout << "Percent of pixels hit in layer " << i+1 << ": " << hitperc << endl;
      cout << "Average percent of unique pixels hits per event: " << hitperc/_nEvt << endl;
      cout << endl << endl << endl;
    }
  //cout << layeruniqueids[0].size() << endl << endl << endl;
  cout << "hitcount: " << hitcount << endl;
  //cout << "number of events: " << _nEvt << endl;
  //double hitperc = static_cast<double>(uniquepix.size()) / static_cast<double>(matters) * 100;
  //cout << "Percentage of pixels that were hit in layer1: " << hitperc << endl;
  _rootfile->Write();
}
