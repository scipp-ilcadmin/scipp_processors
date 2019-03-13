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

static TH2D* totes;
static TFile* _rootfile;
static int _nEvt = 0;

//static vector<int> layers;
static vector<string> pixid;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static vector<string> uniquepix;
static Layers layers(4, PixelGrid(160*100, vector<int>(160*100,0)));
static vector<TH2D*> graphs;
static vector<TH1D*> angles;
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
  totes = new TH2D("totes", "totesmagotes", 100, -100, 100, 100, -100, 100);
  for (int i =0; i < 5; i++)
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
      int layer = idDec( hit )[ILDCellID0::layer];
      int side = idDec ( hit )[ILDCellID0::side];
      int module = idDec (hit)[ILDCellID0::side];
      int posx = (hit->getPosition()[0] + xmin)*100;
      int posy = (hit->getPosition()[1] + ymin)*100;
      totes->Fill(hit->getPosition()[0], hit->getPosition()[1]);
      
      [&] ()
	{
	  int index = layer-1;
	  hitcount++;
	  double theta = (atan2(posy, posx) + M_PI) * 180/M_PI;
	  string id = to_string(posx/step) + to_string(posy/step);
	  pixid.push_back(id);
	  if(std::find(uniquepix.begin(), uniquepix.end(), id) == uniquepix.end())
	    {
	      uniquepix.push_back(id);
	    }
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
      /*switch (layer)
	{
	case(1):
	  if (true)
	    {
	      hitcount++;
	      int posx = (hit->getPosition()[0] + xmin)*100; // indecies for x,y,z components;
	      int posy = (hit->getPosition()[1] + ymin)*100;
	      double theta = (atan2(hit->getPosition()[1], hit->getPosition()[0]) + M_PI) * 180/M_PI; // angles in degrees
	      _angles->Fill(theta);
	      _l1xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	      _xyPos->Fill(hit->getPosition()[0],hit->getPosition()[1]);
	      string id = to_string(posx/step) + to_string(posy/step);
	      if ((posx < 16000 && posy < 16000) && (posx >= 0 && posy >= 0))
		{ 
		  layer1[posx/step][posy/step]++;
		  pixid.push_back(id);
		  if(std::find(uniquepix.begin(), uniquepix.end(), id) == uniquepix.end())
		    {
		      uniquepix.push_back(id);
		    }
		}
	      else 
		{
		  cout << posx << ", " << posy << ":::::::::::Error in 1" << endl;
		}
	    }
	  	case (2):
	  if (side == 2)
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
	  if (side == 2)
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
	  if (side == 2)
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
	  */
    }
}

void TrackerOccupancyAnalysis::check( LCEvent * evt)
{

}

void TrackerOccupancyAnalysis::end()
{

  //cout << " max rad: " << getMax(radii) << " min rad: " << getMin(radii) << endl;
  //cout << " max y: " << getMax(posyVals) << " min y: " << getMin(posyVals) << endl;
  //cout << " max x: " << getMax(posxVals) << " min x: " << getMin(posxVals) << endl;
  /*for (auto str : uniquepix)
    {
      cout << str << endl;
    }
  
  cout << endl << endl << endl << endl;
  cout << pixid.size() << endl;
  for (auto vec : layer1)
    {
      for (auto hit: vec)
	{
	  cout << setw(2) << hit << " ";
	}
      cout << endl;
    }
  cout << endl << endl << endl << endl;
  
  int pixhit=0;
  for(int i = 0; i < layer1[0].size(); i++)
    {
      for(int j = 0; j < layer1[i].size(); j++)
	{
	  if ( layer1[i][j] !=0 )
	    pixhit++;
	}
    }
  */
  cout << "hitcount: " << hitcount << endl;
  int matters = (2*M_PI*73*73) - (2*M_PI*15*15);
  cout << "number of pixels in layer1: " << matters << endl;
  cout << "amount unique pixels hit: " << uniquepix.size() << endl;
  double hitperc = static_cast<double>(uniquepix.size()) / static_cast<double>(matters) * 100;
  cout << "Percentage of pixels that were hit in layer1: " << hitperc << endl;
  //  cout << "pixels hit: " << pixhit << endl;
  //  _rootfile->WriteObject(gr,"gr");
  _rootfile->Write();
}
