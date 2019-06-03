#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * momanalysis.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
 */

#include "momanalysis.h"
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

momanalysis momanalysis;

static TFile* _rootfile;
static TH1D* momzvals;

momanalysis::momanalysis() : Processor("momanalysis") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void momanalysis::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("momanalysis.root", "RECREATE");
  momzvals = new TH1D("momzvals" , "Momentum of Simulated Particles; Momentum (GeV)", 50, -.2, .2);
  _nEvt = 0;
}

void momanalysis::processRunHeader( LCRunHeader* run)
{

}

void momanalysis::processEvent( LCEvent * evt)
{
  ++_nEvt;
  LCCollection* col = evt->getCollection("MCParticle");
  for (int i = 0; i < col->getNumberOfElements(); ++i)
    {
      MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
      double momz = particle->getMomentum()[2];
      momzvals->Fill(momz);
    }
}

void momanalysis::check( LCEvent * evt)
{

}

void momanalysis::end()
{

  cout << "analysis finished" << endl << endl << endl << endl;

  _rootfile->Write();


}
