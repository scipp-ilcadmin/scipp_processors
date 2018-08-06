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


mcparticles mcparticles;

static TFile* _rootfile;
static TH3D* _posxyz;
static int _nEvt=0;
static vector<float> posxvals;
static vector<float> posyvals;
static vector<float> poszvals;


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
  _rootfile = new TFile("mcparts.root","RECREATE");
  _posxyz = new TH3D("posxyz", "posxyz", 250, -1600.0, 1600.0, 250, -1600.0, 1600.0, 250, -1600.0, 1600.0);
  _nEvt = 0 ;
}



void mcparticles::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 
void mcparticles::processEvent( LCEvent * evt ) 
{
  std::vector<std::string> collectionNames = *evt->getCollectionNames();
  LCCollection* mcparts = evt->getCollection("MCParticle");
  for (int i=0; i < mcparts->getNumberOfElements(); ++i)
  {
    MCParticle* hit=dynamic_cast<MCParticle*>(mcparts->getElementAt(i));
    float posx = hit->getEndpoint()[0];
    posxvals.push_back(posx);
    float posy = hit->getEndpoint()[1];
    posyvals.push_back(posy);
    float posz = hit->getEndpoint()[2];
    poszvals.push_back(posz);
    _posxyz->Fill(posx, posy, posz);
    
    }
  _nEvt++;
}


void mcparticles::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void mcparticles::end(){ 
  cout << "number of events: " << _nEvt << endl;
  cout << "xmax: " <<  *max_element(posxvals.begin(), posxvals.end());
  cout << "  xmin:  " <<  *min_element(posxvals.begin(), posxvals.end()) << endl;
  cout << "ymax: " <<  *max_element(posyvals.begin(), posyvals.end());
  cout << "  ymin:  " <<  *min_element(posyvals.begin(), posyvals.end()) << endl;
  cout << "zmax: " <<  *max_element(poszvals.begin(), poszvals.end());
  cout << "  zmin:  " <<  *min_element(poszvals.begin(), poszvals.end()) << endl;
  _rootfile->Write();
  _rootfile->Close();
}
