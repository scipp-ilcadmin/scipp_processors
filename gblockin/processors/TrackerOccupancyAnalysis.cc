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
static vector<TH1D*> posangles;
static vector<TH1D*> negangles;
static vector<TH1D*> posrads;
static vector<TH1D*> negrads;
static int hitcount = 0;
static vector<double> radii;
static TH3D* threedim;
static vector<TH2D*> posxy;
static vector<TH2D*> negxy;
static TH1D* posside;
static TH1D* negside;
static int poshits =0;
static int neghits =0;
static TH1D* totposang;
static TH1D* totnegang;
static TH1D* totposrad;
static TH1D* totnegrad;


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
  totes = new TH2D("totes", "totesmagotes", 1000, -100, 100, 1000, -100, 100);
  posside = new TH1D("posside", "pos z vals", 100, -10, 250);
  negside = new TH1D("negside", "neg z vals", 100, -250, 10);
  totposang = new TH1D("totposang", "angles", 1000, -4, 4);
  totnegang = new TH1D("totnegang", "angles", 1000, -4, 4);
  totposrad = new TH1D("totposrad", "rads", 250, 0, 100);
  totnegrad = new TH1D("totnegrad", "rads", 250, 0, 100);
  for (int i =0; i < 4; i++)
    {
      posxy.push_back(new TH2D(Form("posxy%d ", i+1), "posxy", 2000, -80, 80, 2000, -80, 80));
      negxy.push_back(new TH2D(Form("negxy%d ", i+1), "negxy", 2000, -80, 80, 2000, -80, 80));
      graphs.push_back(new TH2D(Form("layer%d ", i+1), "layers", 1000, -80, 80, 1000, -80, 80));
      posangles.push_back(new TH1D(Form("posangles%d", i+1), "angles", 1000, -4, 4));
      negangles.push_back(new TH1D(Form("negangles%d", i+1), "angles", 1000, -4, 4));
      posrads.push_back(new TH1D(Form("posrads%d", i+1), "rad vals", 250, 0, 100));
      negrads.push_back(new TH1D(Form("negrads%d", i+1), "rad vals", 250, 0, 100));
      
    }  
  _nEvt = 0;

}

void TrackerOccupancyAnalysis::processRunHeader( LCRunHeader* run)
{

}

void TrackerOccupancyAnalysis::processEvent( LCEvent * evt)
{
  //cout << "DO STUFF    " << _nEvt << endl;
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
      double radval = sqrt((hit->getPosition()[0])*(hit->getPosition()[0])+(hit->getPosition()[1])*(hit->getPosition()[1]));
      double theta = atan2(hit->getPosition()[1], hit->getPosition()[0]);
      if (radval > 15.5)
	{
	  hitcount++;
	  totes->Fill(hit->getPosition()[0], hit->getPosition()[1]);
	  threedim->Fill(hit->getPosition()[0],
			 hit->getPosition()[1],
			 hit->getPosition()[2]);
	}
      int layer = idDec( hit )[ILDCellID0::layer];
      //int side = idDec ( hit )[ILDCellID0::side];
      //int module = idDec (hit)[ILDCellID0::side];
      double posz = hit->getPosition()[2];
      if (posz > 0 && radval > 15.5)
	{
	  posxy[layer-1]->Fill(hit->getPosition()[0], hit->getPosition()[1]);
	  posside->Fill(hit->getPosition()[2]);
	  posangles[layer-1]->Fill(theta);
	  totposang->Fill(theta);
	  posrads[layer-1]->Fill(radval);
	  totposrad->Fill(radval);
	  ++poshits;
	}
      if (posz < 0 && radval > 15.5)
	{
	  negxy[layer -1]->Fill(hit->getPosition()[0], hit->getPosition()[1]);
	  negside->Fill(hit->getPosition()[2]);
	  negangles[layer-1]->Fill(theta);
	  totnegang->Fill(theta);
	  negrads[layer-1]->Fill(radval);
	  totnegrad->Fill(radval);
	  ++neghits;
	}
      //int posx = (hit->getPosition()[0] + xmin)*100;
      //int posy = (hit->getPosition()[1] + ymin)*100;
      //threedim->Fill(hit->getPosition()[0], hit->getPosition()[1], hit->getPosition()[2]);
      //totes->Fill(hit->getPosition()[0], hit->getPosition()[1]);
      //int index = layer-1;
      //double theta = (atan2(posy, posx) + M_PI) * 180/M_PI;
      //string id = to_string(posx/step) + to_string(posy/step);
      //pixid.push_back(id);
      //if(std::find(layeruniqueids[index].begin(), layeruniqueids[index].end(), id) == layeruniqueids[index].end())
      //{
      //layeruniqueids[index].push_back(id);
      //}
      //layerpixids[index].push_back(id);
      //angles[index]->Fill(theta);
      //graphs[index]->Fill(hit->getPosition()[0], hit->getPosition()[1]);
      /*if ((posx < 160000 && posy < 160000) && (posx >=0 && posy >= 0))
	{
	  layers[index][posx/step][posy/step]++;
	}
      else
	{
	  cout << "posx: " << posx << ", posy: " << posy << "ERROR" << endl;
	  }*/
    }
}

void TrackerOccupancyAnalysis::check( LCEvent * evt)
{

}

void TrackerOccupancyAnalysis::end()
{

  cout << "analysis finished" << endl << endl << endl << endl;
  //cout << " max rad: " << getMax(radii) << " min rad: " << getMin(radii) << endl << endl << endl;
  //cout << " max y: " << getMax(posyVals) << " min y: " << getMin(posyVals) << endl;
  //cout << " max x: " << getMax(posxVals) << " min x: " << getMin(posxVals) << endl;
  //double matters [4];
  //matters[0] = 100 * ((2*M_PI*75*75) - (2*M_PI*16*16));
  //matters[1] = 100 * ((2*M_PI*76*76) - (2*M_PI*17*17));
  //matters[2] = 100 * ((2*M_PI*77*77) - (2*M_PI*18*18));
  //matters[3] = 100 * ((2*M_PI*78*78) - (2*M_PI*19*19));

  //for (int i = 0; i < 4; i++)
  //{
  //cout << "total hits in layer " << i+1 << ": " << layerpixids[i].size() << endl;
  //cout << "Pixels in layer " << i+1 << ": " << matters[i] << endl;
  //cout << "unique pixels hit in layer " << i+1 << ": " << layeruniqueids[i].size() << endl;
  //double hitperc = static_cast<double>(layeruniqueids[i].size()) / static_cast<double>(matters[i]) * 100;
  //cout << "Percent of pixels hit in layer " << i+1 << ": " << hitperc << endl;
  //cout << "Average percent of unique pixels hits per event: " << hitperc/matters[i] << endl;
  //cout << endl << endl << endl;
  //}
  //cout << layeruniqueids[0].size() << endl << endl << endl;
  cout << "total hitcount: " << hitcount << endl;
  cout << "pos hits: " << poshits << endl << "neg hits: " << neghits << endl;
  cout << "pos and neg added: " << poshits+neghits << endl;
  //cout << "number of events: " << _nEvt << endl;
  //double hitperc = static_cast<double>(uniquepix.size()) / static_cast<double>(matters) * 100;
  //cout << "Percentage of pixels that were hit in layer1: " << hitperc << endl;
  _rootfile->Write();
}
