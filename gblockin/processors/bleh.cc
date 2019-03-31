#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * bleh.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
 */

#include "bleh.h"
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

bleh bleh;

typedef vector<vector<int>> PixelGrid;
typedef vector<PixelGrid> Layers;
typedef vector<string> PixIDs;
typedef vector<PixIDs> layerpixIDs;
typedef vector<vector<double>> VecDubs;

static TH2D* totes;
static TFile* _rootfile;
static int _nEvt = 0;

//static vector<int> layers;
static vector<string> pixid;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static vector<string> uniquepix;
static layerpixIDs layerpixids(5, vector<string>());
static layerpixIDs layeruniqueids(4, vector<string>());
static Layers layers(4, PixelGrid(160*100, vector<int>(160*100,0)));
static VecDubs posxvals(4, vector<double>());
static VecDubs posyvals(4, vector<double>());
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

bleh::bleh() : Processor("bleh") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void bleh::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("yes.root", "RECREATE");
  threedim = new TH3D("threedim", "3-D model", 100, -100, 100, 100, -100, 100, 100, -300, 300);
  totes = new TH2D("totes", "totesmagotes", 1000, -100, 100, 1000, -100, 100);
  for (int i =0; i < 4; i++)
    {
      graphs.push_back(new TH2D(Form("layer%d ", i), "layers", 1000, -80, 80, 1000, -80, 80));
      angles.push_back(new TH1D(Form("angles%d", i), "angles", 100, -10, 370));
    }  
  _nEvt = 0;

}

void bleh::processRunHeader( LCRunHeader* run)
{

}

void bleh::processEvent( LCEvent * evt)
{
  cout << "DO STUFF    " << _nEvt << endl;
  _nEvt++;
  LCCollection* barrelHits = evt->getCollection("SiVertexEndcapHits");
  static const int step = 1;
  static const double xmin = 80;
  static const double ymin = 80;
  for (int i = 0; i < barrelHits->getNumberOfElements(); i++)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      hitcount++;
      CellIDDecoder<SimTrackerHit> idDec( barrelHits );
      int layer = idDec( hit )[ILDCellID0::layer];
      int posx = (hit->getPosition()[0] + xmin)*100;
      int posy = (hit->getPosition()[1] + ymin)*100;

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
	  posxvals[index].push_back(hit->getPosition()[0]);
	  posyvals[index].push_back(hit->getPosition()[1]);
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

void bleh::check( LCEvent * evt)
{

}

void bleh::end()
{

  //int matters = SOME FUCKING MATH BOi
  for (int i = 0; i < 4; i++)
    {
      cout << "total hits in layer " << i+1 << ": " << layerpixids[i].size() << endl;
      cout << "max X and Y in this layer: " << getMax(posxvals[i]) << ", " << getMax(posyvals[i]) << endl;
      cout << "min X and Y in this layer: " << getMin(posxvals[i]) << ", " << getMin(posyvals[i]);
      //cout << "unique pixels hit in layer " << i+1 << ": " << layeruniqueids[i].size() << endl;
      //double hitperc = static_cast<double>(layeruniqueids[i].size()) / static_cast<double>(matters) * 100;
      //cout << "Percent of pixels hit in layer " << i+1 << ": " << hitperc << endl;
      //cout << "Average percent of unique pixels hits per event: " << hitperc/_nEvt << endl;
      cout << endl << endl << endl;
    }

  cout << "hitcount: " << hitcount << endl;
  _rootfile->Write();
}
