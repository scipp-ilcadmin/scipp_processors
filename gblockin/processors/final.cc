#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * final.cc
 * @author Gregory Blockinger
 * October 4th, 2018
 *
 */

#include "final.h"
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

final final;
static TFile* _rootfile;

static TH1D* _zPos;
static TH2D* _xyPos;
static TH2D* _yzPos;
static TH2D* _yPos;
static TH3D* _xyzPos;
static TH2D* _zandr;
static TH3D* _l1xyzPos;
static TH1I* _layers;
static TH1D* _l1radVals;
static TH1D* _l1thetas;
static TH1D* _zpix;
static TH2D* _l1m0;
static int _nEvt = 0;

static vector<TH2D> files;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static vector<pair<double, int>> data;
static vector<int> layers;
static vector<double> xyradius;
static vector<double> l1vals;
static vector<double> thetavals;
static vector<double> modvals;

const double pixsize= 0.005; // pixel size in microns
const double modlen = 9.546;
const double zerox = 13.5182;
const double zeroy = 8.56747;
const int pixnum = static_cast<int>(modlen/pixsize);
static int arr[pixnum];

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

final::final() : Processor("final") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));
}

string str(int i)
{
  return std::to_string(i);
}

void final::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("final.root", "RECREATE");
  //_xyPos = new TH2D("xypos", "(X,Y) Distribution (All layers)", 500, -80, 80, 500, -80, 80);
  //_xyzPos = new TH3D("xyzpos", "(X,Y,Z) Distribution", 124, -80.0, 80.0, 124, -80.0, 80.0, 124, -190.0, 190.0);
  //_l1radVals = new TH1D("l1radVals", "Radial (X,Y) Distribution (All layers)", 100, -100, 100);
  //_l1thetas = new TH1D("l1thetas", "l1thetas", 50, -20, 380);
  //_zandr = new TH2D("zandr", "zandr", 131, -70, 70, 131, -20, 20);
  //_zpix = new TH1D("zpix", "zpix", 131, -65, 65);
  _l1m0 = new TH2D("l1m0", "l1m0", 10, -100, 100, 10, -100, 100);
  /*  for (int i = 0; i < 12; ++i)
    {
      TH2D *test = new TH2D(Form("l1test%d ", i), "test", 500, -80, 80, 500, -100, 100);
      files.push_back(*test);
    }
  */
  _nEvt = 0;
 }

void final::processRunHeader( LCRunHeader* run)
{

}

void final::processEvent( LCEvent * evt)
{
  LCCollection* hits = evt->getCollection("SiVertexBarrelHits");
  int size = hits->getNumberOfElements();
  _nEvt++;
  
  for (int i = 0; i < hits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(hits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( hits );
      double posx = hit->getPosition()[0]; // indecies for x,y,z components;
      double posy = hit->getPosition()[1];
      double posz = hit->getPosition()[2];  //positions are in mm
      //double xyrad = sqrt( (posx*posx) + (posy*posy) );
      int layer = idDec(hit) [ILDCellID0::layer];
      int module = idDec(hit) [ILDCellID0::module];
      //Entry entry = std::make_pair(xyrad, layer);
      switch (layer)
	{
	case 1:
	  {
	   if (module == 0)
	     {
	       //cout <<  "hits in pixel" << arr[dist] << endl;
	       //int dist = static_cast<int>(sqrt((zerox-posx)*(zerox-posx) + (zeroy-posy)*(zeroy-posy))/0.005);
	       //cout << dist << endl;
	       _l1m0->Fill(posx,posy);
	       //arr[dist] += 1;
	     } 
	  }
	case 2:
	  {
	  }
	case 3:
	  {
	  }
	case 4:
	  {
	  }
	case 5:
	  {
	  }
	case 6:
	  {
	  }
	
	}
    }
}

void final::check( LCEvent * evt)
{

}

void final::end()
{
  //cout << "MAX bPosZ: " << getMax(poszVals) << " MIN bPosZ: " << getMin(poszVals) <<  endl;
  //cout << "MAX bPosY: " << getMax(posyVals) << " MIN bPosY: " << getMin(posyVals) <<  endl;
  //cout << "MAX bPosX: " << getMax(posxVals) << " MIN bPosX: " << getMin(posxVals) <<  endl;
  //cout << getMax(modvals) << "    " << getMin(modvals) << endl;
  //double maxx = getMax(posxVals);
  //double maxy = getMax(posyVals);
  //double minx = getMin(posxVals);
  //double miny = getMin(posyVals);
 
  //cout << "MODULE LENGTH: " << endl;
  //cout << sqrt((maxx-minx)*(maxx-minx) + (maxy-miny)*(maxy-miny)) << endl;
  //_xyPos->GetXaxis()->SetTitle("X (mm)");
  //_xyPos->GetYaxis()->SetTitle("Y (mm)");
  //_xyzPos->GetZaxis()->SetTitle("Z (mm)");

  //_l1radVals->SetFillColor(kRed);
  //_l1thetas->SetFillColor(kRed);
  //copy (arr, arr + sizeof(arr) / sizeof(arr[0]), ostream_iterator<short>(cout, "\n"));
  _rootfile->Write();

}
