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

#include "hcalbarrelhits.h"
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


hcalbarrelhits hcalbarrelhits;

static TFile* _rootfile;
static TH2F* _plot;
static TH1F* _momentum;
static TH1F* _pos;
static TH1F* _energy;
static int _nEvt=0;


hcalbarrelhits::hcalbarrelhits() : Processor("hcalbarrelhits") {
  // modify processor description
  _description = "Protype Processor" ;

  //  register steering parameters: name, description, class-variable, default value
  //  registerInputCollection( LCIO::SIMCALORIMETERHIT, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("BeamCalHits") );
  //  registerInputCollection( LCIO::SIMTRACKERHIT, "CollectionName", "Name of MCParticle collection", _colName, std::string("SiTrackerHit"));    
  registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void hcalbarrelhits::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  //  cout << "Initialized\n\n\n " << endl;
  _rootfile = new TFile("example.root","RECREATE");
  _momentum = new TH1F("Momentum", "Momentum Distribution", 50, -0.00177908, 0.002);
  _energy = new TH1F("EnergyDeposited", "Energy Distribution", 50, 0.0, 11000);
  _pos = new TH1F("Position", "Position", 50, -1200.0, 1200.0);
  _nEvt = 0 ;
}



void hcalbarrelhits::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 
void hcalbarrelhits::processEvent( LCEvent * evt ) 
{ 
  vector<float> posVals;
  vector<int> cell1Vals;
  vector<int> cell2Vals;
  vector<int> nmcVals;
  vector<float> eVals;
  std::vector<std::string> collectionNames = *evt->getCollectionNames();
  LCCollection* hcalcol = evt->getCollection("HCalBarrelHits");
  for (int i=0; i < hcalcol->getNumberOfElements(); ++i)
    {
      SimCalorimeterHit* hit2=dynamic_cast<SimCalorimeterHit*>(hcalcol->getElementAt(i));
      float pos = *hit2->getPosition();
      posVals.push_back(pos);
      int cell1 = hit2->getCellID0();
      cell1Vals.push_back(cell1);
      int cell2 = hit2->getCellID1();
      cell2Vals.push_back(cell2);
      int nmcconts = hit2->getNMCContributions();
      nmcVals.push_back(nmcconts);
      float energy = hit2->getEnergy();
      eVals.push_back(energy);
    }
  _nEvt++;
  cout << endl;
}


void hcalbarrelhits::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void hcalbarrelhits::end(){ 
  cout << "number of events: " << _nEvt << endl;

  //  _momentum->Write();  
  _rootfile->Write();

    }
