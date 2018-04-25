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

#include "test.h"
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


test test;

static TFile* _rootfile;
static TH2F* _plot;
static TH1F* _histo;

static int _nEvt=0;


test::test() : Processor("test") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void test::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    cout << "Initialized " << endl;

    _rootfile = new TFile("test.root","RECREATE");

    _plot = new TH2F("hh", "Hit-Hit HeatMap", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _histo = new TH1F("energy", "Energy Distribution", 300, 0, 550);

    _nEvt = 0 ;
}



void test::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 
void test::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    _nEvt++;
    for(int i=0; i < col->getNumberOfElements(); ++i){
      MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
      int pid=particle->getPDG();
      cout << pid << endl;
      double energy = particle->getEnergy();
      double mass = particle->getMass();
      double charge = particle->getCharge();
      //      double momentum = particle->getMomentum();
      cout << "energy: ";
      cout << energy << endl;
      cout << "mass: ";
      cout << mass << endl;
      cout << "charge: ";
      cout << charge << endl;
      //      cout << "momentum: ";
      //      cout << momentum << endl;
      _histo->Fill(energy);
    }
    cout << endl;
}


void test::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void test::end(){ 
  cout << "number of events: " << _nEvt << endl;
  
  _rootfile->Write();
}
