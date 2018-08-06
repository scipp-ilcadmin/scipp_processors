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
#include <algorithm>

#include <EVENT/SimTrackerHit.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;


hcalbarrelhits hcalbarrelhits;

static TFile* _rootfile;
static TH1F* _momentum;
static TH1F* _posx;
static TH1F* _posy;
static TH1F* _posz;
static TH3D* _posxyz;
static TH1F* _energy;
static int _nEvt=0;

static vector<float> posxvals;
static vector<float> posyvals;
static vector<float> poszvals;

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
  _rootfile = new TFile("hcal.root","RECREATE");
  _posxyz = new TH3D("posxyz", "posxyz", 50, -2500.0, 2500.0, 50, -2500.0, 2500.0, 50, -2500.0, 2500.0);
  _nEvt = 0 ;
}



void hcalbarrelhits::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 
void hcalbarrelhits::processEvent( LCEvent * evt ) 
{ 
  std::vector<std::string> collectionNames = *evt->getCollectionNames();
  LCCollection* hcalcol = evt->getCollection("HcalBarrelHits");
  for (int i=0; i < hcalcol->getNumberOfElements(); ++i)
    {
      SimCalorimeterHit* hit2=dynamic_cast<SimCalorimeterHit*>(hcalcol->getElementAt(i));
      float posx = hit2->getPosition()[0];
      float posy = hit2->getPosition()[1];
      float posz = hit2->getPosition()[2];
      posxvals.push_back(posx);
      posyvals.push_back(posy);
      poszvals.push_back(posz);
      _posxyz->Fill(posx, posy, posz);
    }
  _nEvt++;
  cout << endl;
}


void hcalbarrelhits::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void hcalbarrelhits::end(){ 
  cout << "number of events: " << _nEvt << endl;

  cout << "xmax: " <<  *max_element(posxvals.begin(), posxvals.end());
  cout << "  xmin:  " <<  *min_element(posxvals.begin(), posxvals.end()) << endl;
  cout << "ymax: " <<  *max_element(posyvals.begin(), posyvals.end());
  cout << "  ymin:  " <<  *min_element(posyvals.begin(), posyvals.end()) << endl;
  cout << "zmax: " <<  *max_element(poszvals.begin(), poszvals.end());
  cout << "  zmin:  " <<  *min_element(poszvals.begin(), poszvals.end()) << endl;

  //  _momentum->Write();  
  _rootfile->Write();
  _rootfile->Close();
    }
