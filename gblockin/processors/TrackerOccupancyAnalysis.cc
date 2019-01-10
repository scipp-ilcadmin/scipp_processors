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


static TFile* _rootfile;
static TH2D* _l1xyPos;
static TH2D* _l2xyPos;
static TH2D* _l3xyPos;
static TH2D* _l4xyPos;
static TH2D* _xyPos;
static TH1D* _angles;
static int _nEvt = 0;

static vector<int> layers;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static vector<double> angles;
static vector<vector<int>> layer1(16, vector<int>(16, 0));
static vector<vector<int>> layer2(16, vector<int>(16, 0));
static vector<vector<int>> layer3(16, vector<int>(16, 0));
static vector<vector<int>> layer4(16, vector<int>(16, 0));
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
  _l1xyPos = new TH2D("l1xypos", "l1xypos", 100, -100.0, 100.0, 100, -100.0, 100.0);
  _l2xyPos = new TH2D("l2xypos", "l2xypos", 100, -100.0, 100.0, 100, -100.0, 100.0);
  _l3xyPos = new TH2D("l3xypos", "l3xypos", 100, -100.0, 100.0, 100, -100.0, 100.0);
  _l4xyPos = new TH2D("l4xypos", "l4xypos", 100, -100.0, 100.0, 100, -100.0, 100.0);
  _xyPos = new TH2D(  "xypos",   "xyPos",   100, -100.0, 100.0, 100, -100.0, 100.0);
  _angles = new TH1D("angles", "angles", 100, -10, 370);
  _nEvt = 0;

}

void TrackerOccupancyAnalysis::processRunHeader( LCRunHeader* run)
{

}

void TrackerOccupancyAnalysis::processEvent( LCEvent * evt)
{
  _nEvt++;
  //LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHits");
  static const double xmin = 75;
  static const double ymin = 75;
  static const int step = 10;
  
  for (int i = 0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( endcapHits );
      int layer = idDec( hit )[ILDCellID0::layer];
      int side = idDec ( hit )[ILDCellID0::side];
      layers.push_back(layer);
      switch (layer)
	{
	case(1):
	  if (side == 1)
	    {
	      hitcount++;
	      int posx = hit->getPosition()[0] + xmin; // indecies for x,y,z components;
	      int posy = hit->getPosition()[1] + ymin;
	      posxVals.push_back(posx);
	      double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
	      _angles->Fill(theta);
	      _l1xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	      _xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	      if ((posx < 160 && posy < 160) && (posx >= 0 && posy >= 0))
		{ 
		  layer1[posx/step][posy/step]++;
		}
	      else 
		{
		  cout << posx << ", " << posy << ":::::::::::Error in 1" << endl;
		}
	    }
	case (2):
	  if (side == 1)
	    {
	      hitcount++;
	      int posx = hit->getPosition()[0] + xmin; // indecies for x,y,z components;
	      int posy = hit->getPosition()[1] + ymin;
	      double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
              _angles->Fill(theta);
	      _l2xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
       	      _xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	      if ((posx < 160 && posy < 160) && (posx >= 0 && posy >= 0))
		{
		  layer2[posx/step][posy/step]++;
		}
	      else
		{
		  cout << posx << ", " << posy << ":::::::::::Error in 2" << endl;
		} 
	    }
	case (3):
	  if (side == 1)
	    {
              hitcount++;
              int posx = hit->getPosition()[0] + xmin; // indecies for x,y,z components;
	      int posy = hit->getPosition()[1] + ymin;
	      double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
              _angles->Fill(theta);
              _l3xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
       	      _xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
              if ((posx < 160 && posy < 160) && (posx >= 0 && posy >= 0))
                {
                  layer3[posx/step][posy/step]++;
                }
              else
                {
                  cout << posx << ", " << posy << ":::::::::::Error in 3" << endl;
                }
	    }
	case (4):
	  if (side == 1)
            {
	      hitcount++;
              int posx = hit->getPosition()[0] + xmin; // indecies for x,y,z components;
	      int posy = hit->getPosition()[1] + ymin;
	      double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
              _angles->Fill(theta);
              _l4xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
       	      _xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
              if ((posx < 160 && posy < 160) && (posx >= 0 && posy >= 0))
                {
		  layer4[posx/step][posy/step]++;
                }
	      else
                {
                  cout << posx << ", " << posy << ":::::::::::Error in 4" << endl;
                }
            }
	}
    }
}

void TrackerOccupancyAnalysis::check( LCEvent * evt)
{

}

void TrackerOccupancyAnalysis::end()
{

  //cout << " max energy: " << getMax(eDepVals) << " MinEnergy: " << getMin(eDepVals) << endl;
  //cout << " max x: " << getMax(posxVals) << " min x: " << getMin(posxVals) << endl;
  //cout << " max y: " << getMax(posyVals) << " min y: " << getMin(posyVals) << endl;
  //cout << " max z: " << getMax(poszVals) << " min z: " << getMin(poszVals) << endl;
  /*  for (auto vec : layer1)
    {
      for (auto hit: vec)
	{
	  cout << setw(2) << hit << " ";
	}
      cout << endl;
    }
  cout << endl << endl << endl << endl;
    for (auto vec : layer2)
    {
      for (auto hit : vec)
	{
	  cout << setw(2) << hit << " ";
	}
      cout << endl;
    }
  cout << endl << endl << endl << endl;
  for (auto vec : layer3)
    {
      for (auto hit : vec)
        {
          cout << setw(2) << hit << " ";
        }
      cout << endl;
    }
  cout << endl << endl << endl << endl;
  for (auto vec : layer4)
    {
      for (auto hit : vec)
        {
          cout << setw(2) << hit << " ";
        }
      cout << endl;
    }
  */
  cout << "hitcount: " << hitcount << endl;
  //  _rootfile->WriteObject(gr,"gr");
  _rootfile->Write();
}
