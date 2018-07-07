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

#include "ecalbarrelhits.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <algorithm>

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


ecalbarrelhits ecalbarrelhits;

static TFile* _rootfile;
static TH2F* _plot;
static TH1F* _momentum;
static TH1F* _pos;
static TH1F* _energy;
static int _nEvt=0;


ecalbarrelhits::ecalbarrelhits() : Processor("ecalbarrelhits") {
  // modify processor description
  _description = "Protype Processor" ;

  //  register steering parameters: name, description, class-variable, default value
  //  registerInputCollection( LCIO::SIMCALORIMETERHIT, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("BeamCalHits") );
  //  registerInputCollection( LCIO::SIMTRACKERHIT, "CollectionName", "Name of MCParticle collection", _colName, std::string("SiTrackerHit"));    
  registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void ecalbarrelhits::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;
  //  cout << "Initialized\n\n\n " << endl;
  _rootfile = new TFile("example.root","RECREATE");
  _nEvt = 0 ;
}



void ecalbarrelhits::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 
void ecalbarrelhits::processEvent( LCEvent * evt ) 
{ 
  vector<float> eVals;
  vector<double> posVals;
  vector<int> cell1Vals;
  vector<int> cell2Vals;
  vector<int> nmcVals;
  std::vector<std::string> collectionNames = *evt->getCollectionNames();
  LCCollection* ecalcol = evt->getCollection("EcalBarrelHits");
  for (int i=0; i < ecalcol->getNumberOfElements(); ++i)
  {
    SimCalorimeterHit* hit=dynamic_cast<SimCalorimeterHit*>(ecalcol->getElementAt(i));
    float energy = hit->getEnergy();
    eVals.push_back(energy);
    float pos = *hit->getPosition();
    posVals.push_back(pos);
    int cell0 = hit->getCellID0();
    cell1Vals.push_back(cell0);
    int cell1 = hit->getCellID1();
    cell2Vals.push_back(cell1);
    int nmcconts = hit->getNMCContributions();
    nmcVals.push_back(nmcconts);  
  }
  float maxEnergy = *max_element(posVals.begin(), posVals.end());
  cout << maxEnergy << endl;
  _nEvt++;
  cout << endl;
}


void ecalbarrelhits::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void ecalbarrelhits::end()
{ 
  cout << "number of events: " << _nEvt << endl;
  _rootfile->Write();
}
