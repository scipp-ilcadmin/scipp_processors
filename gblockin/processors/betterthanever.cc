#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * Ok, so I like C++11. Unfortunately,
 * Marlin is uilt with ansi C, so the processor
 * constructor freaks out about the string that is
 * passed to it as an argument. The above two lines
 * fix that issue, allowing our code to be compatible
 * with ansi C class declarations.
 * Big thanks to Daniel Bittman for helping me fix this.
 */

/*
 * author Gregory Blockinger
 * May 11th, 2019
 */

#include "betterthanever.h"
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


betterthanever betterthanever;

typedef vector<vector<int>> PixelGrid;
typedef vector<PixelGrid> Layers; // creating structures to store the pixels and their hits in
typedef vector<Layers> Barrel;
typedef vector<string> PixIDs;
typedef vector<PixIDs> layerpixIDs;

static TFile* _rootfile;
static int _nEvt; 
static int nhit;
static Barrel barrel{};
static double rsuba[5][30][2] = {{{0.0}}};
static double thetas[5][30][1];
static TH2D* newxylays[5];
static int modhits[5][30] = {{0}};
static TH1D* posmomzvals;
static TH1D* negmomzvals;
static TH1D* phigone[5];
static TH1D* phigoneneg[5];
static TH1D* phigonemid[5];
static TH1D* phigonepos[5];
static vector<double> vecids;
static TH2D* newmods[5][30];
static TH2D* ogxyplane;
static TH2D* trimmedogplane;
static TH2D* newxyplane;
static vector<int> particleids;
static TH1D* zgone[5];
static double avgradvals[5] = {0};
static int nhitlayer[5] = {0};
static int hithit;
static vector<string> pixidsbreh;
static TH1D* zgonetot;
static TH1D* phigonetot;
static TH1D* posphigone;
static TH1D* negphigone;
static int unqhits = 0;
static TH2D* layer1;

template<typename T>
static T getMax(vector<T> &vec) //this gets the greatest value in a vector
{
  return *max_element(vec.begin(), vec.end());   
}
template<typename T>
static T getMin(vector<T> &vec)  //this one gets the smallest value
{
  return *min_element(vec.begin(), vec.end());
}

betterthanever::betterthanever() : Processor("betterthanever") 
{
  // modify processor description
    _description = "Protype Processor" ;
}


void betterthanever::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    _rootfile = new TFile("first.root","RECREATE");  //Any histograms would be created in this block
    barrel.emplace_back(12, PixelGrid(9.8*10, vector<int>(125*10,0)));
    //barrel.emplace_back(12, PixelGrid(14, vector<int>(126,0)));
    //barrel.emplace_back(18, PixelGrid(14, vector<int>(126,0)));
    //barrel.emplace_back(24, PixelGrid(14, vector<int>(126,0)));
    //barrel.emplace_back(30, PixelGrid(14, vector<int>(126,0)));    
    ogxyplane = new TH2D("ogxyplane", "Original Plane", 500, -100, 100, 500, -100, 100);
    trimmedogplane = new TH2D("BarrelHits", "Module View of Barrel Array; X-component (cm); Y-component (cm)", 500, -70, 70, 500, -70, 70);
    newxyplane = new TH2D("newxyplane", "Rotated Plane", 500, -12, 12, 500, -80, 80);
    zgonetot = new TH1D("zgonetot", "Collapsed in z", 628, -3.5, 3.5);
    phigonetot = new TH1D("phigonetot", "Collapsed in Phi", 130, -65, 65);
    posphigone = new TH1D("posphigone", "Collapsed in Phi", 200, -2, 65);
    negphigone = new TH1D("negphigone", "Collapsed in Phi", 200, -65, 2);
    layer1 = new TH2D("layer1","layer1",500, -20, 20, 500, -20, 20);
    posmomzvals = new TH1D("posmomzvals", "Z momentum values; momentum in GeV", 50, .18, .5);
    negmomzvals = new TH1D("negmomzvals", "Z momentum values; momentum in GeV", 50, -.5, -.18);
    for (int i=0; i < 5; ++i)
      {
	phigone[i] = new TH1D(Form("phigone%d", i+1), "Collapsed in Phi; Z value of hit", 126, -64, 64);
	phigoneneg[i] = new TH1D(Form("phigoneneg%d", i+1), "Collapsed in Phi, Lesser Vals", 42, -65, -21);
	phigonemid[i] = new TH1D(Form("phigonemid%d", i+1), "Collapsed in Phi, Middle Vals", 42,  -21, 21);
	phigonepos[i] = new TH1D(Form("phigonepos%d", i+1), "Collapsed in Phi, Higher Vals;", 42,   21 ,65);
	zgone[i] = new TH1D(Form("zgone%d", i+1), "Collapsed in Z; Aziumthal Value of hit", 628, -3.5, 3.5); 
	newxylays[i] = new TH2D(Form("newxylays%d", i+1), "breh breh", 500, -12, 12, 500, -80, 80);
	//for(int j=0; j <30; ++j)
	//{
	//newmods[i][j] = (new TH2D(Form("newmod%d%d",i,j), "module", 500, -70, 70, 500, -130, 130));
	//}
      }
    _nEvt = 0;
    nhit = 0;
}

void betterthanever::processRunHeader( LCRunHeader* run)
{ 

} 


void betterthanever::processEvent( LCEvent * evt )
{ 
  LCCollection* barrelhits = evt->getCollection( "SiVertexBarrelHits" ); //Whatever we put in quotes is the type of collection we will be looking at
  _nEvt++;                                                               //e.g. SiVertexBarrelHits, SiVertexEndcaphits, MCParticle, etc...
  int step = 1;
    for(int i=0; i < barrelhits->getNumberOfElements(); ++i)  //loop for each hit in the collection
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec( hit )[ILDCellID0::layer];   //info about hits stored in these types
	//int side = idDec ( hit )[ILDCellID0::side];
	int module = idDec (hit)[ILDCellID0::module];
	//int sensor = idDec (hit)[ILDCellID0::sensor];
	//int subdet = idDec (hit)[ILDCellID0::subdet];
	double posx = hit->getPosition()[0]; //position stored in array, indexed for each coordinate
	double posy = hit->getPosition()[1];
	double posz = hit->getPosition()[2];
	if (posz > 0)
	  {
	    posphigone->Fill(posz); 
	  }
	if (posz < 0)
	  {
	    negphigone->Fill(posz);
	  }
	MCParticle* particle=hit->getMCParticle();
	int particleid = particle->getPDG();
	double momz = particle->getMomentum()[2];
	ogxyplane->Fill(posx, posy);
	double radval = sqrt((posx)*(posx)+(posy)*(posy));
	//double azimuth = atan2(posy, posx);
	//avgradvals[layer-1]+=radval;
	rsuba[layer-1][module][0]+=posx;
	rsuba[layer-1][module][1]+=posy;
	//++nhitlayer[layer-1];
	double azimuth = atan2(posy, posx);
	zgonetot->Fill(azimuth);
	//zgone[layer-1]->Fill(azimuth);
	phigonetot->Fill(posz);
	++modhits[layer-1][module];
	if (radval > 12.0 && radval < 16.0)
	  {
	    layer1->Fill(posx,posy);
	    trimmedogplane->Fill(posx,posy);
	    rsuba[layer-1][module][0]+=posx;
	    rsuba[layer-1][module][1]+=posy;
	    nhitlayer[layer-1]++;
	    //double azimuth = atan2(posy, posx);
	    zgone[layer-1]->Fill(azimuth);
	    phigone[layer-1]->Fill(posz);
	    modhits[layer-1][module]++;
	  }
	if (radval > 20.9 && radval < 24.99)
	  {
	    trimmedogplane->Fill(posx,posy);
	    rsuba[layer-1][module][0]+=posx;
	    rsuba[layer-1][module][1]+=posy;
	    nhitlayer[layer-1]++;
	    //double azimuth = atan2(posy, posx);
	    zgone[layer-1]->Fill(azimuth);
	    phigone[layer-1]->Fill(posz);
	    modhits[layer-1][module]++;
	  }
	if (radval > 33.0 && radval < 37.75)
	  {
	    trimmedogplane->Fill(posx,posy);
	    rsuba[layer-1][module][0]+=posx;
	    rsuba[layer-1][module][1]+=posy;
	    nhitlayer[layer-1]++;
	    //double azimuth = atan2(posy, posx);
	    zgone[layer-1]->Fill(azimuth);
	    phigone[layer-1]->Fill(posz);
	    modhits[layer-1][module]++;
	  }
	if (radval > 45.0 && radval < 50.3)
	  {
	    trimmedogplane->Fill(posx,posy);
	    rsuba[layer-1][module][0]+=posx;
	    rsuba[layer-1][module][1]+=posy;
	    nhitlayer[layer-1]++;
	    //double azimuth = atan2(posy, posx);
	    zgone[layer-1]->Fill(azimuth);
	    phigone[layer-1]->Fill(posz);
	    modhits[layer-1][module]++;
	  }
	if (radval > 57.6 && radval < 63.7)
	  {
	    trimmedogplane->Fill(posx,posy);
	    rsuba[layer-1][module][0]+=posx;
	    rsuba[layer-1][module][1]+=posy;
	    nhitlayer[layer-1]++;
	    //double azimuth = atan2(posy, posx);
	    zgone[layer-1]->Fill(azimuth);
	    phigone[layer-1]->Fill(posz);
	    modhits[layer-1][module]++;
	    }	
	if(std::find(particleids.begin(), particleids.end(), particleid) == particleids.end())                
	  {
	    particleids.push_back(particleid);
	    //cout << "particle id: " << particleid << endl;
	  }
	if (momz < -.18)
	  {
	    negmomzvals->Fill(momz);
	  }
	if (momz > 0.18)
	  {
	    posmomzvals->Fill(momz);
	  }
	phigone[layer-1]->Fill(posz);
	if (posz <= -21)
	  {
	    phigoneneg[layer-1]->Fill(posz);
	  }
	if (posz > -21 && posz <= 21)
	  {
	    phigonemid[layer-1]->Fill(posz);
	  }
	if (posz > 21)
	  {
	    phigonepos[layer-1]->Fill(posz);
	  }
      }
    for (int i=0; i < 5; ++i)
      {                                           
	for (int j = 0; j < 30; ++j)  
	  { 
	    //rsuba[i][j][0]/=modhits[i][j];
	    //rsuba[i][j][1]/=modhits[i][j];
	    thetas[i][j][0] = atan(rsuba[i][j][0]/rsuba[i][j][1]);
	  }                                     
      }
    /*for(int i=0; i < barrelhits->getNumberOfElements(); ++i)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec(hit)[ILDCellID0::layer];
	int module = idDec(hit)[ILDCellID0::module];
	double x = hit->getPosition()[0];
	double y = hit->getPosition()[1];
	double z = hit->getPosition()[2];
	double radval = sqrt((x)*(x)+(y)*(y));
	double newx = ((x * cos(thetas[layer-1][module][0])) - (y * sin(thetas[layer-1][module][0])));
	double newy = ((y * cos(thetas[layer-1][module][0])) + (x * sin(thetas[layer-1][module][0])));
	//newxyplane->Fill(newx, newy);
	newxylays[layer-1]->Fill(newx, newy);
	if (radval > 12.0 && radval < 16.0)
	  {
	    //double newx = ((x * cos(thetas[layer-1][module][0])) - (y * sin(thetas[layer-1][module][0])));
	    //double newy = ((y * cos(thetas[layer-1][module][0])) + (x * sin(thetas[layer-1][module][0])));
	    //newmods[layer-1][module]->Fill(newx, z);
	    if (newx > 0 && newx < 10)
	      {
		newx += (9.8/2);
		z += (125/2);
		newxyplane->Fill(newx, newy);
		string id = to_string(static_cast<int>((newx*100))) + to_string(static_cast<int>((z*100)));
		if(std::find(pixidsbreh.begin(), pixidsbreh.end(), id) == pixidsbreh.end())
		  {
		    pixidsbreh.push_back(id);
		  }
	    //cout << id << endl;
		++hithit;
		if (((newx*10 <= 98) && (z*10 <= 1260)) && ((newx >=0) && (z >=0)))
		  {
		    barrel[layer-1][module][(newx*10)/step][(z*10)/step]++;
		  }
	      }
	  }
	if (radval > 20.9 && radval < 24.99)
          {
            double newx = ((x * cos(thetas[layer-1][module][0])) - (y * sin(thetas[layer-1][module][0])));
            double newy = ((y * cos(thetas[layer-1][module][0])) + (x * sin(thetas[layer-1][module][0])));
            newmods[layer-1][module]->Fill(newx, z);
            newxyplane->Fill(newx, newy);
          }

	if (radval > 33.0 && radval < 37.75)
          {
            double newx = ((x * cos(thetas[layer-1][module][0])) - (y * sin(thetas[layer-1][module][0])));
            double newy = ((y * cos(thetas[layer-1][module][0])) + (x * sin(thetas[layer-1][module][0])));
            newmods[layer-1][module]->Fill(newx, z);
            newxyplane->Fill(newx, newy);
          }

	if (radval > 45.0 && radval < 50.3)
          {
            double newx = ((x * cos(thetas[layer-1][module][0])) - (y * sin(thetas[layer-1][module][0])));
            double newy = ((y * cos(thetas[layer-1][module][0])) + (x * sin(thetas[layer-1][module][0])));
            newmods[layer-1][module]->Fill(newx, z);
            newxyplane->Fill(newx, newy);
          }

	if (radval > 57.6 && radval < 63.7)
          {
            double newx = ((x * cos(thetas[layer-1][module][0])) - (y * sin(thetas[layer-1][module][0])));
            double newy = ((y * cos(thetas[layer-1][module][0])) + (x * sin(thetas[layer-1][module][0])));
            newmods[layer-1][module]->Fill(newx, z);
	    newxyplane->Fill(newx, newy);
	    }
	}*/
}


void betterthanever::check( LCEvent * evt )
{ 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void betterthanever::end()
{
  /*cout << "Amount of pixels in a layer 1 module of pixel size 10x10 micrometers: " << 9.8*126*1000000/(10*10) << endl;
  cout << "Amount of pixels in a layer 1 module of pixel size 50x50 micrometers: " << 9.8*126*1000000/(50*50) << endl;
  cout << "Amount of pixels in a layer 1 module of pixel size 20x20 micrometers: " << 9.8*126*1000000/(20*20) << endl;
  cout << "Amount of pixels in a layer 1 module of pixel size 100x100 micrometers: " << 9.8*126*1000000/(100*100) << endl;
  cout << "Amount of pixels in a layer 2-5 module of pixel size 10x10 micrometers: " << 14*126*1000000/(10*10) << endl;
  cout << "Amount of pixels in a layer 2-5 module of pixel size 50x50 micrometers: " << 14*126*1000000/(50*50) << endl;
  cout << "Amount of pixels in a layer 2-5 module of pixel size 20x20 micrometers: " << 14*126*1000000/(20*20) << endl;
  cout << "Amount of pixels in a layer 2-5 module of pixel size 100x100 micrometers: " << 14*126*1000000/(100*100) << endl;
  */

  nhit = nhitlayer[0] + nhitlayer[1] + nhitlayer[2] + nhitlayer[3] + nhitlayer[4];
  int totsbruh=0;
  //for (int i =0; i < 12; ++i)
  //{
  //totsbruh += modhits[0][i];
  //}
  //cout << "layer 1, all mods hits in modhits: " << totsbruh << endl;
  //cout << "layer1 hits in rot block : " << hithit << endl;
  //cout << "uniqueids in pixidsbreh: " << pixidsbreh.size() << endl;
  cout << "nhit: " << nhit << endl;
  cout << "_nEvt: " << _nEvt << endl;

  for (int i=0; i < 5; ++i)
    {
      cout << "hits in layer " << i+1 << ": " << nhitlayer[i] << endl; 
      for (int j = 0; j < 30; ++j)
	{
	  // cout << "layer: " << i+1 << " module: " << j+1 << " number of hits: " << modhits[i][j];
	  //cout << endl;
	}
      cout << endl;
    }

  _rootfile->Write();

}
