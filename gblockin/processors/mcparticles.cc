#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/* 
 * Marlin is built with ansi C, so the processor
 * constructor freaks out about the string that is
 * passed to it as an argument. The above two lines
 * fix that issue, allowing our code to be compatible
 * with ansi C class declarations.
 */

/*
 * @author Gregory Blockinger
 * July 4, 2018
 */

#include "mcparticles.h"
#include "scipp_ilc_utilities.h"
#include <iostream>

#include <EVENT/SimTrackerHit.h>
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


mcparticles mcparticles;

static TFile* _rootfile;
static TH2F* _plot;
static TH1F* _momentum;
static TH1F* _pos;
static TH1F* _energy;
static int _nEvt=0;


mcparticles::mcparticles() : Processor("mcparticles") {
  // modify processor description
  _description = "Protype Processor" ;

  //  register steering parameters: name, description, class-variable, default value
  //  registerInputCollection( LCIO::SIMCALORIMETERHIT, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("BeamCalHits") );
  registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void mcparticles::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  //  cout << "Initialized\n\n\n " << endl;
  _rootfile = new TFile("example.root","RECREATE");
  _momentum = new TH1F("Momentum", "Momentum Distribution", 50, -0.00177908, 0.002);
  _energy = new TH1F("EnergyDeposited", "Energy Distribution", 50, 0.0, 11000);
  _pos = new TH1F("Position", "Position", 50, -1200.0, 1200.0);
  _nEvt = 0 ;
}



void mcparticles::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 
void mcparticles::processEvent( LCEvent * evt ) 
{
  vector<float> evals;
  vector<float> pvals;
  vector<int> cell1vals;
  vector<int> cell2vals;
   
  std::vector<std::string> collectionNames = *evt->getCollectionNames();
  LCCollection* bcalcol = evt->getCollection("BeamCalHits");
  for (int i=0; i < bcalcol->getNumberOfElements(); ++i)
  {
    //    cout << "THIS IS FOR SIMCALHITS:::::::" << endl;
    SimCalorimeterHit* hit=dynamic_cast<SimCalorimeterHit*>(bcalcol->getElementAt(i));
    float energy = hit->getEnergy();
    evals.push_back(energy);
    float pos = *hit->getPosition();
    pvals.push_back(pos);
    int cell0 = hit->getCellID0();
    cell1vals.push_back(cell0);
    int cell1 = hit->getCellID1();
    cell2vals.push_back(cell1);
    //    int nmcparts = hit->getNMCParticles(); !!DEPRECATED!!
    int nmcconts = hit->getNMCContributions();
    //    printf("Energy: %d, Position: %f, Cell1: %d, Cell2: %d, Number of MCParticles: %d\n", energy, pos, cell0, cell1, nmcconts);
    }
  _nEvt++;
  cout << endl;
}


void mcparticles::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void mcparticles::end(){ 
  cout << "number of events: " << _nEvt << endl;
  _rootfile->Write();
}
