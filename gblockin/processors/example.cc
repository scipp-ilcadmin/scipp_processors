#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * Ok, so I like C++11. Unfortunately,
 * Marlin is built with ansi C, so the processor
 * constructor freaks out about the string that is
 * passed to it as an argument. The above two lines
 * fix that issue, allowing our code to be compatible
 * with ansi C class declarations.
 * Big thanks to Daniel Bittman for helping me fix this.
 */

/*
 * author Christopher Milke
 * April 5, 2016
 */

#include "example.h"
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


example example;

static TFile* _rootfile;
static TH2D* mods;
static TH3D* threedim;
static TH2D* ogmod;
static TH2D* newmod;
static vector<double> ogxvals;
static vector<double> ogyvals;
static vector<double> newxvals;

//static const double normal[1.0, 0, 1.0];
//static vector<string> vec;
static int _nEvt=0;
static int nhit = 0;


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

example::example() : Processor("example") 
{
    // modify processor description
    _description = "Protype Processor" ;
}


void example::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    _rootfile = new TFile("stahp.root","RECREATE");
    mods = new TH2D("mods", "mods", 1000, -100, 100, 1000, -100, 100);
    ogmod =  new TH2D("ogmod" , "ogmod" , 300, -100, 100, 300, -100, 100);
    newmod = new TH2D("newmod", "newmod", 300, -100, 100, 300, -100, 100);
    threedim = new TH3D("threedim", "3-D model", 100, -50, 50, 100, -50, 50, 100, -150, 150);
    _nEvt = 0 ;
}



void example::processRunHeader( LCRunHeader* run)
{ 
} 


void example::processEvent( LCEvent * evt ) { 
    LCCollection* barrelhits = evt->getCollection( "SiVertexBarrelHits" );
    _nEvt++;
    //double rsubi[0.0, 0.0, 0.0];
    double rsuba[3] = {0.0, 0.0, 0.0};
    for(int i=0; i < barrelhits->getNumberOfElements(); ++i)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec( hit )[ILDCellID0::layer];
	//int side = idDec ( hit )[ILDCellID0::side];
	int module = idDec (hit)[ILDCellID0::module];
	//int sensor = idDec (hit)[ILDCellID0::sensor];
	//int subdet = idDec (hit)[ILDCellID0::subdet];
	//string all = to_string(layer) + ", " + to_string(module) + ", " + to_string(sensor) + ", " + to_string(subdet);
	double posx = hit->getPosition()[0];
	double posy = hit->getPosition()[1];
	double posz = hit->getPosition()[2];
	if (layer ==1 && module ==1)
	  {
	    ogxvals.push_back(posx);
	    ogyvals.push_back(posy);
	    ogmod->Fill(posx, posy);
	    nhit++;; 
	    rsuba[0] += posx;
	    rsuba[1] += posy;
	    rsuba[2] += posz;
	  }
      }
    for(int i=0; i < barrelhits->getNumberOfElements(); i++)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int module = idDec( hit )[ILDCellID0::module];
	int layer = idDec(  hit )[ILDCellID0::layer];
	double posx = hit->getPosition()[0];
	double posy = hit->getPosition()[1];
	double posz = hit->getPosition()[2];
	if (layer == 1 && module == 1)
	  {
	    rsuba[0] = rsuba[0]/nhit;
	    rsuba[1] = rsuba[1]/nhit;
	    rsuba[2] = rsuba[2]/nhit;
	    newmod->Fill(posx - ( posx - rsuba[0] + posz - rsuba[2]) / sqrt(posx*posx), 0.0);
	    newxvals.push_back(posx - (posx - rsuba[0] + posz - rsuba[2]) / sqrt(posx*posx));
	  }

      }
}


void example::check( LCEvent * evt )
{ 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void example::end()
{ 
  cout << "number of events: " << _nEvt << endl;
  cout << "max og xval: "  << getMax(ogxvals);
  cout << " min og xval: " << getMin(ogxvals) << endl;
  cout << " max new xval: " << getMax(newxvals);
  cout << " min new xval: " << getMin(newxvals) << endl;
  _rootfile->Write();
  cout << nhit << endl;
}
