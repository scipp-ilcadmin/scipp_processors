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

#include "NewProc.h"
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


NewProc NewProc;

static TFile* _rootfile;
static TH2F* _hitmap;


NewProc::NewProc() : Processor("NewProc") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void NewProc::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile(_root_file_name.c_str(),"RECREATE");
    _hitmap = new TH2F("hitmap","Hit Distribution",300.0,-150.0,150.0,300.0,-150.0,150.0);

    // usually a good idea to
    //printParameters() ;

    _nRun = 0 ;
    _nEvt = 0 ;

}



void NewProc::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void NewProc::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    LCCollection* col = evt->getCollection( _colName ) ;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
	int Nelectrons = 0;
	int Npositrons = 0;
	int NOthers = 0;

        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           SimCalorimeterHit* hit = dynamic_cast<SimCalorimeterHit*>( col->getElementAt(hitIndex) );
	   
           const float* pos = hit->getPosition();
	   const int id = 11;

	   // cout << "This is type " << id << endl;
	   switch (id){
	   case 11: Nelectrons++;
	     break;
	   case -11: Npositrons++;
	     break;
	   default: NOthers++;
	   }

           _hitmap->Fill(pos[0],pos[1]);
	   // cout << "The Z position is " << pos[2] << endl;
        }
	cout << "The event number is " << _nEvt << ", and there are " << nElements << " elements." << endl;
	cout << Nelectrons << " of these are electrons, and " << Npositrons << " are positrons." <<endl;
    }

    _nEvt ++ ;
}



void NewProc::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void NewProc::end(){ 
    _rootfile->Write();
}
