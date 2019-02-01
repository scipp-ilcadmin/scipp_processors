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


static TFile* _rootfile;
static TH2D* _l1xyPos;
static TH2D* _xyPos;
static TH1D* _angles;
static TH2D* _xyPos2;
static TH1D* _angles2;
static int _nEvt = 0;

static vector<int> layers;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static vector<double> angles;
static vector<vector<int>> layer1(25, vector<int>(25, 0));
static vector<vector<int>> layer2(25, vector<int>(25, 0));
static int hitcount = 0;

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
  _xyPos = new TH2D(  "xypos",   "xyPos",   100, -30.0, 30.0, 100, -30.0, 30.0);
  _xyPos2 = new TH2D("xypos2", "xypos2", 100, -30.0, 30.0, 100, -30.0, 30.0);
  _angles = new TH1D("angles", "angles", 100, -10, 370);
  _angles2 = new TH1D("angels2", "angles2", 100, -10, 370);
  _nEvt = 0;

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
  static const int step = 10;
  
  for (int i = 0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( endcapHits );
      int layer = idDec ( hit )[ILDCellID0::layer];
      int side = idDec ( hit )[ILDCellID0::side];
      int module = idDec ( hit )[ILDCellID0::module];
      int posx = hit->getPosition()[0] + xmin;
      int posy = hit->getPosition()[1] + ymin;
      switch (layer)
	{
	case(2):
	  if (module % 2 == 0)
	    {
	      hitcount++;
	      double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
	      _angles->Fill(theta);
	      _xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	      if ((posx < 320 && posy < 320) && (posx >= 0 && posy >= 0))
		{ 
		  layer1[posx/step][posy/step]++;
		}
	      else 
		{
		  cout << posx << ", " << posy << ":::::::::::Error in 1" << endl;
		}
	    }
	  if (module % 2 != 0)
	    {
	      hitcount++;
	      double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
	      _angles2->Fill(theta);
	      _xyPos2->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	      if ((posx < 320 && posy < 320) && (posx >= 0 && posy >= 0))
                {
                  layer2[posx/step][posy/step]++;
                }
              else
                {
                  cout << posx << ", " << posy << ":::::::::::Error in 1" << endl;
                }

	    }
	}
    }
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
  for (auto vec : layer1)
    {
      for (auto hit: vec)
	{
	  cout << setw(3) << hit << " ";
	}
      cout << endl;
    }
  cout << endl << endl << endl << endl;
  cout << "hitcount: " << hitcount << endl;
  //  _rootfile->WriteObject(gr,"gr");
  _rootfile->Write();
}
