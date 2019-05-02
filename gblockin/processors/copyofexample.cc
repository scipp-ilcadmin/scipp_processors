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
 * April 5, 2019
 */

#include "copyofexample.h"
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


copyofexample copyofexample;

typedef vector<vector<int>> PixelGrid;
typedef vector<PixelGrid> Layers;
typedef vector<string> PixIDs;
typedef vector<PixIDs> layerpixIDs;

static TFile* _rootfile;
static int _nEvt=0;
static int nhit = 0;
static TH2D* totes;
static vector<TH2D*> ogmods;
static vector<TH2D*> newmods;
static layerpixIDs layerpixids(5, vector<string>());
static layerpixIDs layeruniqueids(5, vector<string>());
static Layers layers(5, PixelGrid(160*100, vector<int>(160*100,0)));

static double rsuba [12][3];
static double thetavals [12][1];

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

copyofexample::copyofexample() : Processor("copyofexample") 
{
    // modify processor description
    _description = "Protype Processor" ;
}


void copyofexample::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    _rootfile = new TFile("bread.root","RECREATE");
    totes = new TH2D("totes", "magotes", 1000, -30, 30, 1000, -30, 30);
    for (int i = 1; i < 13; i++)
      {
	ogmods.push_back(new TH2D(Form("ogmod%d ", i), "mods", 1000, -30, 30, 1000, -30, 30));
	newmods.push_back(new TH2D(Form("newmod%d ", i), "mods", 1000, -30, 30, 1000, -30, 30));
      }

    _nEvt = 0 ;
}



void copyofexample::processRunHeader( LCRunHeader* run)
{ 
} 


void copyofexample::processEvent( LCEvent * evt ) { 
    LCCollection* barrelhits = evt->getCollection( "SiVertexBarrelHits" );
    _nEvt++;
    for(int i=0; i < barrelhits->getNumberOfElements(); ++i)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec( hit )[ILDCellID0::layer];
	//int side = idDec ( hit )[ILDCellID0::side];
	int module = idDec (hit)[ILDCellID0::module];
	//int sensor = idDec (hit)[ILDCellID0::sensor];
	//int subdet = idDec (hit)[ILDCellID0::subdet];
	if (layer == 1)
	  {
	    double posx = hit->getPosition()[0];
	    double posy = hit->getPosition()[1];
	    double posz = hit->getPosition()[2];
	    double radval = sqrt(posx*posx + posy*posy);
	    if (radval > 12)
	      {
		totes->Fill(posx, posy);
		++nhit;
		[&] ()
		  {
		    int index = module;
		    rsuba[index][0] += posx;
		    rsuba[index][1] += posy;
		    rsuba[index][2] += posz;
		    ogmods[index]->Fill(posx, posy);
		  }();
	      }
	  }
      }
    for (int j = 0; j < 12; ++j)
      {
	rsuba[j][0] /= nhit;
	rsuba[j][1] /= nhit;
	thetavals[j][0] = atan2(rsuba[j][0], rsuba[j][1]);
	cout << "rsuba x: " << rsuba[j][0] << " rsuba y: " << rsuba[j][1];
	cout << " theta for module " << j << ": " << thetavals[j][0] << endl; 
      }
    for (int i = 0; i < barrelhits->getNumberOfElements(); ++i)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec ( hit ) [ILDCellID0::layer];
	int module = idDec (hit ) [ILDCellID0::module];
	if (layer == 1)
	  {
	    double posx = hit->getPosition()[0];
	    double posy = hit->getPosition()[1];
	    double posz = hit->getPosition()[2];
	    double radval = sqrt(posx*posx + posy*posy);
	    if (radval > 12)
	      {
		double newx = posx * cos(thetavals[module][0]) - posy * sin(thetavals[module][0]);
		double newy = posy * cos(thetavals[module][0]) + posx * sin(thetavals[module][0]);
		[&]()
		  {
		    int index = module;
		    newmods[index]->Fill(newx, newy);
		  }();
	      }
	  }
      }
}


void copyofexample::check( LCEvent * evt )
{ 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void copyofexample::end()
{
  cout << "number of hits total: " << nhit << endl;
 
  _rootfile->Write("Update");

}
