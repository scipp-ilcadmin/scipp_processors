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
static TH2D* twodimexy;
static TH2D* twodimbxy;
static TH1D* onedimbzvals;
static TH1D* onedimbphivals;
static TH1D* onedimeposrvals;
static TH1D* onedimenegrvals;
static TH1D* onedimepostvals;
static TH1D* onedimenegtvals;

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
  threedim = new TH3D("threedim", "3-D Model of Occupancy in Vertex Detector;;;z-axis / BeamLine Direction", 200, -80, 80, 200, -80, 80, 200, -250, 250);
  twodimbxy = new TH2D("twodimbxy", "2-D View of Back Scatter hits in the Barrel Array;x (cm); y (cm)", 1000, -80, 80, 1000, -80, 80);
  twodimexy = new TH2D("twodimexy", "2-D View of Back Scatter hits in the End-Cap Arrays;x (cm); y (cm)", 200, -80, 80, 200, -80, 80);
  onedimbzvals = new TH1D("onedimbzvals", "Z Component of Barrel Array Hits", 200, -100, 100);
  onedimbphivals = new TH1D("onedimbphivals", "Angular Distribution of Barrel Array Hits", 200, -100, 100);
  onedimeposrvals = new TH1D("onedimeposrvals", "Positive End-Cap Array Radial Distribution of Hits", 200, -10, 100);
  onedimenegrvals = new TH1D("onedimenegrvals", "Negative End-Cap Array Radial Distribution of Hits", 200, -10, 100);
  onedimepostvals = new TH1D("onedimepostvals", "Positive End-Cap Array Angular Distribution of Hits", 200, -100, 100);
  onedimenegtvals = new TH1D("onedimenegtvals", "Negative End-Cap Array Angular Distribution of Hits", 200, -100, 100);
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
      double radval = sqrt((posx)*(posx)+(posy)*(posy));
      double phi = atan2(posy,posx);
      threedim->Fill(posx, posy, posz);
      if (posz > 0)
	{
	  onedimeposrvals->Fill(radval); 
	}
      if (posz < 0)
	{
	  onedimenegrvals->Fill(radval);
	}
      if (radval > 0)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimexy->Fill(posx, posy);
	}
    }
  for (int i =0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      double posx = hit->getPosition()[0];
      double posy = hit->getPosition()[1];
      double posz = hit->getPosition()[2];
      double radval = sqrt((posx)*(posx)+(posy)*(posy));
      threedim->Fill(posx,posy,posz);
      if (radval > 13.0 && radval < 16.0)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	}
      if (radval > 21.4 && radval < 24.99)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	}
      if (radval > 33.5 && radval < 37.75)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	}
      if (radval > 45.5 && radval < 50.3)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	}
      if (radval > 58.0 && radval < 63.7)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	}

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
