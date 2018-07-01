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
 * June30, 2018
 */

#include "example.h"
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


example example;

static TFile* _rootfile;
static TH2F* _plot;
static TH1F* _momentum;

static int _nEvt=0;


example::example() : Processor("example") {
  // modify processor description
  _description = "Protype Processor" ;

  // register steering parameters: name, description, class-variable, default value
  //  registerInputCollection( LCIO::SIMCALORIMETERHIT, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("BeamCalHits") );
  //  registerInputCollection( LCIO::SIMTRACKERHIT, "CollectionName", "Name of MCParticle collection", _colName, std::string("SiTrackerHit"));    
  registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void example::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  //  cout << "Initialized\n\n\n " << endl;
  _rootfile = new TFile("example.root","RECREATE");
  _momentum = new TH1F("Momentum", "Momentum Distribution", .01 ,0, 2);
  _nEvt = 0 ;
}



void example::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 
void example::processEvent( LCEvent * evt ) 
{ 
  std::vector<std::string> collectionNames = *evt->getCollectionNames();
  LCCollection* bcalcol = evt->getCollection("BeamCalHits");
  LCCollection* sitrackcol = evt->getCollection("SiTrackerBarrelHits");
  /*for (int i=0; i < bcalcol->getNumberOfElements(); ++i)
  {
    //    cout << "THIS IS FOR SIMCALHITS:::::::" << endl;
    SimCalorimeterHit* hit=dynamic_cast<SimCalorimeterHit*>(bcalcol->getElementAt(i));
    float energy = hit->getEnergy();
    float pos = *hit->getPosition();
    int cell0 = hit->getCellID0();
    int cell1 = hit->getCellID1();
    //    int nmcparts = hit->getNMCParticles(); !!DEPRECATED!!
    int nmcconts = hit->getNMCContributions();
    //    printf("Energy: %d, Position: %f, Cell1: %d, Cell2: %d, Number of MCParticles: %d\n", energy, pos, cell0, cell1, nmcconts);
    }*/
  //  cout << endl << endl << endl << endl << endl;
  for (int i=0; i < 481; ++i)
  {
    //    cout << "THIS IS FOR SIMTRACKERHITS::::::" << endl;
    SimTrackerHit* hit2=dynamic_cast<SimTrackerHit*>(sitrackcol->getElementAt(i));
    int sicell0 = hit2->getCellID0();
    int sicell1 = hit2->getCellID1();
    double sipos = *hit2->getPosition();
    //    float edx = hit2->getEdx();
    float edep = hit2->getEDep();
    float time = hit2->getTime();
    float mom = *hit2->getMomentum();
    float pathlength = hit2->getPathLength();
    _momentum->Fill(mom);
    //    printf("Cell1: %d, Cell2: %d, Position: %f, Energy Deposited: %f, Time: %f, Momentum: %f Path Length: \n", sicell0, sicell1, sipos, edep, time, mom, pathlength );
    }
  _nEvt++;
  cout << endl;
}


void example::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void example::end(){ 
  cout << "number of events: " << _nEvt << endl;

  _momentum->Write();  
  _rootfile->Write();
}
