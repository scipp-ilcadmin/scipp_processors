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

#include "EventAnalysis_L.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <sstream>

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


EventAnalysis_L EventAnalysis_L;

stringstream _hInfo;
TFile* _rootfile2;
TH2F* _hitmap2;


EventAnalysis_L::EventAnalysis_L() : Processor("EventAnalysis_L") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
}



void EventAnalysis_L::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile2 = new TFile("hitmap.root","RECREATE");
    _hitmap2 = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);

    // usually a good idea to
    printParameters() ;


    _hInfo << "Get some" << " strings " << endl;
    _nRun = 0 ;
    _nEvt = 0 ;

}



void EventAnalysis_L::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void EventAnalysis_L::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );

           const float* pos = hit->getPosition();
           _hitmap2->Fill(pos[0],pos[1]);
        } 
	
	/*for(  LCIterator it( evt, "Tracks" ) ;  SimCalorimeterHit* hit  = it.next()  ; ) {
	  std::cout << trk->getTrackState( TrackState::AtIP ) << std::endl  ;
	  }*/

    }
    //_hInfo << "canyouprintinProcess?";
    //_hInfo << _nEvt;
    _nEvt ++ ;
    //std::cout << _nEvt << std::endl;
    //std::cout << _hInfo.str() << std::endl;
}



void EventAnalysis_L::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void EventAnalysis_L::end(){ 
  _rootfile2->Write();
  //std::cout << _nEvt << std::endl;
  //_hInfo << "END STUFF" << std::endl;
  //std::cout << _hInfo.str();	
  
  
}
