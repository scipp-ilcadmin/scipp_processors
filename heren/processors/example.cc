

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

#include "example.h"
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


example example;

static TFile* _rootfile;
static TH2F* _plot;
static TH1F* _histo;
static TH1F* _momentum;
static TH1F* _scalar;
static TH1F* _vector;
static TH1F* _mass;

static int _nEvt=0;


example::example() : Processor("example") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}


void example::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    cout << "Initialized " << endl;

    _rootfile = new TFile("example.root","RECREATE");

    _plot = new TH2F("hh", "Hit-Hit HeatMap", 300.0, -150.0, 150.0, 300.0, -150.0, 150.0);
    _histo = new TH1F("energy", "Energy Distribution", 300, 0, 550);
    _momentum = new TH1F("momentum","Momentum Distribution",300,0,500);
    _scalar = new TH1F("scalar","Scalar Distribution",300,0,500);
    _vector = new TH1F("vector","Vector Distribution",300,0,50);
    _mass = new TH1F("mass","Mass Distribution",300,0,5000);
    _nEvt = 0 ;
}



void example::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 
void example::processEvent( LCEvent * evt ) { 
    LCCollection* col = evt->getCollection( _colName );
    _nEvt++;
    double totalMag=0;
    double totalMomx=0;
    double totalMomy=0;
    double totalMomz=0;
    double totalEnergy=0;

    /* double max_e=0.0, max_p=0.0;
    for(int i=0; i < col->getNumberOfElements(); ++i){
       MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
       int id=particle->getPDG();
       double energy=particle->getPDG();
       if(id== 11 && energy>max_e) max_e=energy;
       if(id==-11 && energy>max_p) max_p=energy;
       }*/


    for(int i=0; i < col->getNumberOfElements(); ++i){
      MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
      double energy=particle->getEnergy();
      int pid=particle->getPDG();      
      if(particle->getGeneratorStatus()!=1) continue;
      //      if(energy == max_e || energy == max_p) continue;
      //      if(pid==1000022 || pid==1000023 || pid==1000025 || pid==1000035 || pid==1000015 || pid==2000015 /* pid==12 || pid==14 || pid==16 || pid==18*/) continue; //exclude four main types of neutralinos, neutrinos and stauons.
      _histo->Fill(energy);
      double mom[3];
      mom[0]=particle->getMomentum()[0];
      mom[1]=particle->getMomentum()[1];
      mom[2]=particle->getMomentum()[2];
      double newX=0.0;
      double newE=0.0;
      //      scipp_ilc::transform_to_lab(mom[0],energy,newX,newE);
      // if(scipp_ilc::get_hitStatus(newX,mom[1],mom[2])>2)continue;
      double xyz=sqrt(pow(mom[0],2)+pow(mom[1],2)+pow(mom[2],2));//norm of the totoal momentum
      double cth=abs(mom[2]/xyz);//cosine of the angle between momentum vector and z-direction momentum
      if(cth>0.9999) continue;
      // if(cth>0.9999 || cth<-0.9999) continue;
      double mag=sqrt(pow(mom[0],2)+pow(mom[1],2));//exclude z momentum
      _momentum->Fill(mag);
      totalMag +=mag;
      totalMomx +=mom[0];
      totalMomy +=mom[1];
      totalMomz +=mom[2];
      totalEnergy +=energy;
    }
    _scalar->Fill(totalMag);
    double totalMom=sqrt(pow(totalMomx,2)+pow(totalMomy,2));//exclude z momentum
    double totalMass=sqrt(pow(totalEnergy,2)-pow(totalMomx,2)-pow(totalMomy,2)-pow(totalMomz,2));
    _vector->Fill(totalMom);
    _mass->Fill(totalMass);
    
}


void example::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void example::end(){ 
  cout << "number of events: " << _nEvt << endl;
  
  _rootfile->Write();
}
