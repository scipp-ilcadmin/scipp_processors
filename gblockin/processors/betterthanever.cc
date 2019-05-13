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
static TH1D* phigone;
static vector<double> vecids;
static TH2D* newmods[5][30];
static TH2D* ogxyplane;
static TH2D* newxyplane;
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
    _rootfile = new TFile("tits.root","RECREATE");  //Any histograms would be created in this block
    barrel.emplace_back(12, PixelGrid(10, vector<int>(126,0)));
    barrel.emplace_back(12, PixelGrid(14, vector<int>(126,0)));
    barrel.emplace_back(18, PixelGrid(14, vector<int>(126,0)));
    barrel.emplace_back(24, PixelGrid(14, vector<int>(126,0)));
    barrel.emplace_back(30, PixelGrid(14, vector<int>(126,0)));
    phigone = new TH1D("phigone","collapsed in phi",100, -100, 100);
    ogxyplane = new TH2D("ogxyplane", "L", 500, -100, 100, 500, -100, 100);
    newxyplane = new TH2D("newxyplane", "L", 500, -100, 100, 500, -100, 100);
    for (int i=0; i < 5; ++i)
      {
	for(int j=0; j <30; ++j)
	  {
	    newmods[i][j] = (new TH2D(Form("newmod%d%d",i,j), "module", 500, -200, 200, 500, -200, 200));
	  }
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
       	ogxyplane->Fill(posx, posy);
	++nhit;
	rsuba[layer-1][module][0]+=posx;
	rsuba[layer-1][module][1]+=posy;
	//rsuba[layer][module][2]++;
	phigone->Fill(posz);
      }
    for(int i=0; i < barrelhits->getNumberOfElements(); ++i)  //loop for each hit in the collection
      {
        SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
        CellIDDecoder<SimTrackerHit> idDec(barrelhits);
        int layer = idDec( hit )[ILDCellID0::layer];   //info about hits stored in these types 
	int module = idDec (hit)[ILDCellID0::module];
	//rsuba[layer-1][module][0] = rsuba[layer-1][module][0]/nhit;
        //rsuba[layer-1][module][1] = rsuba[layer-1][module][1]/nhit;
        //rsuba[layer][module][2] /= nhit;
      }


    for (int i=0; i < 5; ++i)
      {                                           
	for (int j = 0; j < 30; ++j)  
        { 
          rsuba[i][j][0]/=nhit;   
          rsuba[i][j][1]/=nhit;
          thetas[i][j][0] = atan2(rsuba[i][j][0], rsuba[i][j][1]);       
	}                                     
      }
    for(int i=0; i < barrelhits->getNumberOfElements(); ++i)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec(hit)[ILDCellID0::layer];
	int module = idDec(hit)[ILDCellID0::module];
	double x = hit->getPosition()[0];
	double y = hit->getPosition()[1];
	double z = hit->getPosition()[2];
	double newx = ((x * cos(thetas[layer-1][module][0])) - (y * sin(thetas[layer-1][module][0])));
	double newy = ((y * cos(thetas[layer-1][module][0])) + (x * cos(thetas[layer-1][module][0])));
	newmods[layer-1][module]->Fill(newx, newy);
	newxyplane->Fill(newx, newy);
      }
}


void betterthanever::check( LCEvent * evt )
{ 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void betterthanever::end()
{
  cout << "nhit: " << nhit << endl;
  cout << "_nEvt: " << _nEvt << endl;
  /*for (int i=0; i < 5; ++i)
    {
      for (int j = 0; j < 30; ++j)
	{
	  rsuba[i][j][0]/=nhit;
	  rsuba[i][j][1]/=nhit;
	  thetas[i][j][0] = atan2(rsuba[i][j][0], rsuba[i][j][1]);
	  cout << "layer: " << i << ", module: " << j <<", x-val: "
	       << rsuba[i][j][0] << ", y-val: " << rsuba[i][j][1]
	       << ", theta val:" <<thetas[i][j][0]; 
	  //cout << ", z-val: " << rsuba[i][j][2];
	  cout << endl;
	}
      cout << endl;
  }*/

  _rootfile->Write();

}
