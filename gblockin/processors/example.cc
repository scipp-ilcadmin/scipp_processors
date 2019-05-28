#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * example.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
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
static TH3D* threedim;


example::example() : Processor("example") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void example::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("example.root", "RECREATE");
  threedim = new TH3D("threedim", "3-D model of hits in Vertex Detector", 200, -100, 100, 200, -100, 100, 200, -300, 300);
  _nEvt = 0;

}

void example::processRunHeader( LCRunHeader* run)
{

}

void example::processEvent( LCEvent * evt)
{
  ++_nEvt;
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHits");
  for (int i = 0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      double posx = hit->getPosition()[0];
      double posy = hit->getPosition()[1];
      double posz = hit->getPosition()[2];
      threedim->Fill(posx, posy, posz);
    }
  for (int i =0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      double posx = hit->getPosition()[0];
      double posy = hit->getPosition()[1];
      double posz = hit->getPosition()[2];
      threedim->Fill(posx,posy,posz);
    }
}

void example::check( LCEvent * evt)
{

}

void example::end()
{

  cout << "analysis finished" << endl << endl << endl << endl;
  _rootfile->Write();
}
