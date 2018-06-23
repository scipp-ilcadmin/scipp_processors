#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/*
 * OccupancyAnalysis.cc
 *
 * 
 * @author Gregory Blockinger
 *
 *
 */


#include "OccupancyAnalysis.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;

static TFile* _rootfile;
static TH1F* _histo;
static TH2F* _plot;
static int _nEvt=0;


OccupancyAnalysis OccupancyAnalysis;

OccupancyAnalysis::OccupancyAnalysis() : Processor("OccupancyAnalysis")
{
  _description = "Occupancy Analysis Processor";

  registerInputCollection( LCIO::MCPARTICLE, "CollectionName", "Name of MCPArticle collection", _colName, std::string("BeamCalHits"));
    }

void OccupancyAnalysis::check(LCEvent *evt) 
{
  //nothing here could be used to fill checkplots in rconstruction processor
}

void OccupancyAnalysis::processRunHeader(LCRunHeader *run)
{

}

void  OccupancyAnalysis::init()
{

  _rootfile = new TFile("occupancy.root", "RECREATE");

  _histo = new TH1F("TH1D", "z position", 2000, 0, 2000); 
  _plot = new TH2F("TH2D", "X Y Hit Occupancy Over All Layers", 1000, -1750, 1750, 1000, -1750, 1750);

  
}

void OccupancyAnalysis::end()
{
  _rootfile->Write();
}

void OccupancyAnalysis::processEvent(LCEvent* evt)
{
  LCCollection* col = evt->getCollection(_colName);
  _nEvt++;
  for (int i=0; i < col->getNumberOfElements(); i++)
  {
    SimCalorimeterHit* hit=dynamic_cast<SimCalorimeterHit*>(col->getElementAt(i));
    

  }
}
