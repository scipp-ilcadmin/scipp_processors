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

#include "ecalendcaphits.h"
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


ecalendcaphits ecalendcaphits;

static TFile* _rootfile;
static TH2F* _plot;
static TH1F* _momentum;
static TH1F* _pos;
static TH1F* _energy;
static int _nEvt=0;


ecalendcaphits::ecalendcaphits() : Processor("ecalendcaphits") {
  // modify processor description
  _description = "Protype Processor" ;

  //  register steering parameters: name, description, class-variable, default value
  //  registerInputCollection( LCIO::SIMCALORIMETERHIT, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("BeamCalHits") );
  //  registerInputCollection( LCIO::SIMTRACKERHIT, "CollectionName", "Name of MCParticle collection", _colName, std::string("SiTrackerHit"));    
  registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void ecalendcaphits::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  //  cout << "Initialized\n\n\n " << endl;
  _rootfile = new TFile("example.root","RECREATE");
  _momentum = new TH1F("Momentum", "Momentum Distribution", 50, -0.00177908, 0.002);
  _energy = new TH1F("EnergyDeposited", "Energy Distribution", 50, 0.0, 11000);
  _pos = new TH1F("Position", "Position", 50, -1200.0, 1200.0);
  _nEvt = 0 ;
}



void ecalendcaphits::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 
void ecalendcaphits::processEvent( LCEvent * evt ) 
{ 
  vector<float> eDepVals;
  vector<double> posVals;
  std::vector<std::string> collectionNames = *evt->getCollectionNames();
  LCCollection* ecalcol = evt->getCollection("EcalEndcapHits");
  for (int i=0; i < ecalcol->getNumberOfElements(); ++i)
  {
    SimTrackerHit* hit2=dynamic_cast<SimTrackerHit*>(ecalcol->getElementAt(i));
    int sicell0 = hit2->getCellID0();
    int sicell1 = hit2->getCellID1();
    double sipos = *hit2->getPosition();
    posVals.push_back(sipos);
    //    float edx = hit2->getEdx();
    float edep = hit2->getEDep();
    eDepVals.push_back(edep * 1000000000);
    float time = hit2->getTime();
    float mom = *hit2->getMomentum();
    float pathlength = hit2->getPathLength();
    _pos->Fill(sipos);
    _momentum->Fill(mom);
    _energy->Fill(edep * 1000000000);
    //    printf("Cell1: %d, Cell2: %d, Position: %f, Energy Deposited: %f, Time: %f, Momentum: %f Path Length: \n", sicell0, sicell1, sipos, edep, time, mom, pathlength );
    }
  _nEvt++;
  cout << endl;
}


void ecalendcaphits::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ecalendcaphits::end(){ 
  cout << "number of events: " << _nEvt << endl;

  //  _momentum->Write();  
  _rootfile->Write();

    }
