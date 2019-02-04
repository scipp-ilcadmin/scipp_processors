#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * barrel.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
 */

#include "barrel.h"
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

barrel barrel;

typedef vector<vector<int>> PixelGrid;
typedef vector<PixelGrid> Layers;

static TFile* _rootfile;
static int _nEvt = 0;

//static vector<int> layers;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static Layers layers(5, PixelGrid(150*100, vector<int>(150*100, 0)));  // vector size = (2*shift)*ma 
static vector<TH2D*> graphs;
static vector<TH1D*> angles;
static vector<vector<int>> l1inner(150*100, vector<int>(150*100, 0));
static vector<vector<int>> l1outer(150*100, vector<int>(150*100, 0));
static int l1m1hitcount = 0;
static int l1m2hitcount = 0;

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

barrel::barrel() : Processor("barrel") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void barrel::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("barrel.root", "RECREATE");
  _nEvt = 0;
  for (int i =0; i < 6; i++)
    {
      graphs.push_back(new TH2D(Form("layer%d ", i), "test", 1000, -80, 80, 1000, -80, 80)); 
      angles.push_back(new TH1D(Form("angles%d", i), "fuckyou", 100, -10, 370));
    }
}

void barrel::processRunHeader( LCRunHeader* run)
{

}

void barrel::processEvent( LCEvent * evt)
{
  _nEvt++;
  //LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexBarrelHits");
  static const double ymin = 75;
  static const double xmin = 75;
  static const int step = 1;

  for (int i = 0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( endcapHits );
      int layer = idDec ( hit )[ILDCellID0::layer];
      int side = idDec ( hit )[ILDCellID0::side];
      int module = idDec ( hit )[ILDCellID0::module];
      int posx = (hit->getPosition()[0] + xmin)*100;
      int posy = (hit->getPosition()[1] + ymin)*100;

      [&] () 
	{
	  int index = layer -1;
	  double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees                                                                                  
	  angles[index]->Fill(theta);
	  graphs[index]->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	  if ((posx < 160000 && posy < 160000) && (posx >= 0 && posy >= 0))
	    {
	      layers[index][posx/step][posy/step]++;
	    }
	  else
	    {
	      cout << posx << ", " << posy << ":::::::::::Error in 1" << endl;
	    }
	}();

      /*switch (layer)
	{
	case(1):
	  if (module % 2 == 0)
	    {
	      l1m1hitcount++;
	      double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
	      _angles->Fill(theta);
	      _xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	      if ((posx < 160 && posy < 160) && (posx >= 0 && posy >= 0))
		{ 
		  l1inner[posx/step][posy/step]++;
		}
	      else 
		{
		  cout << posx << ", " << posy << ":::::::::::Error in 1" << endl;
		}
	    }
	  if (module % 2 !=0)
	    {
	      l1m2hitcount++;
              double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
	      _angles2->Fill(theta);
              _xyPos2->Fill(hit->getPosition()[0],hit->getPosition()[1]);
              if ((posx < 160 && posy < 160) && (posx >= 0 && posy >= 0))
                {
                  l1outer[posx/step][posy/step]++;
                }
              else
                {
                  cout << posx << ", " << posy << ":::::::::::Error in 1" << endl;
                }
      
	    }

    	}
      */}
}

void barrel::check( LCEvent * evt)
{

}

void barrel::end()
{

  //cout << " max energy: " << getMax(eDepVals) << " MinEnergy: " << getMin(eDepVals) << endl;
  //cout << " max x: " << getMax(posxVals) << " min x: " << getMin(posxVals) << endl;
  //cout << " max y: " << getMax(posyVals) << " min y: " << getMin(posyVals) << endl;
  //cout << " max z: " << getMax(poszVals) << " min z: " << getMin(poszVals) << endl;
  /*  for (auto vec : l1inner)
    {
      for (auto hit: vec)
	{
	  cout << setw(3) << hit << " ";
	}
      cout << endl;
    }
  */
  cout << endl << endl << endl << endl;
  cout << "l1m1hitcount: " << l1m1hitcount << endl;
  cout << _nEvt << endl;
  _rootfile->Write();
}
