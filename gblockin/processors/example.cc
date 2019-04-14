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
static vector<string> vec;
static int _nEvt=0;


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
    threedim = new TH3D("threedim", "3-D model", 100, -50, 50, 100, -50, 50, 100, -150, 150);
    _nEvt = 0 ;
}



void example::processRunHeader( LCRunHeader* run)
{ 
} 


void example::processEvent( LCEvent * evt ) { 
    LCCollection* barrelhits = evt->getCollection( "SiVertexBarrelHits" );
    _nEvt++;
    for(int i=0; i < barrelhits->getNumberOfElements(); ++i)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec( hit )[ILDCellID0::layer];
	int side = idDec ( hit )[ILDCellID0::side];
	int module = idDec (hit)[ILDCellID0::module];
	int sensor = idDec (hit)[ILDCellID0::sensor];
	int subdet = idDec (hit)[ILDCellID0::subdet];
	string all = to_string(layer) + ", " + to_string(module) + ", " + to_string(sensor) + ", " + to_string(subdet);
	double posx = hit->getPosition()[0];
	double posy = hit->getPosition()[1];
	double posz = hit->getPosition()[2];
	mods->Fill(posx,posy);
	if (layer ==1 && module ==1)
	  {
	    threedim->Fill(posx, posy, posz); 
	  }
	//	if(std::find(vec.begin(), vec.end(), subdet) == vec.end())
	//  {
	//    vec.push_back(subdet);
	//  }
      }
}


void example::check( LCEvent * evt )
{ 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void example::end()
{ 
  cout << "number of events: " << _nEvt << endl;
  _rootfile->Write();
}
