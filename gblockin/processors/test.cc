#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/*
 * test.cc
 * @author Gregory Blockinger
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

static TFile* _rootfile;
static TH1F* _histo;
static TH2F* _plot;
static int _nEvt=0;


test test;

test::test() : Processor("test")
{
  _description = "Occupancy Analysis Processor";

  //  registerInputCollection( LCIO::MCPARTICLE, "CollectionName", "Name of MCPArticle collection", _mccolName, std::string("NoCLUE_MCPART"));
  registerInputCollection(LCIO::SIMCALORIMETERHIT, "SimCalorimeterHits", "BeamCalHits", _bcolName, std::string("NoCLUE_BCALLHOPEFULLY"));
}

void test::check(LCEvent *evt) 
{
  //nothing here could be used to fill checkplots in rconstruction processor
}

void test::processRunHeader(LCRunHeader *run)
{

}

void  test::init()
{
  
}

void test::end()
{

}

void test::processEvent(LCEvent* evt)
{
  LCCollection* col = evt->getCollection( _colName );
  _nEvt++;
  double mom[3];

  for(int i=0; i < col->getNumberOfElements(); ++i){
    MCParticle* particle=dynamic_cast<MCParticle*>(col->getElementAt(i));
    int pid=particle->getPDG();
    mom[0]=particle->getMomentum()[0];
    mom[1]=particle->getMomentum()[1];
    mom[2]=particle->getMomentum()[2];
    cout << "Momentum(X,Y,Z): (" << mom[0] << ", " << mom[1] << ", "<< mom[2] << ")" << endl;
    double mag = sqrt(pow(mom[0],2) + pow(mom[1],2) + pow(mom[2],2));
    cout << "Momentum Magnitude: " << mag << endl;
    double energy = particle->getEnergy();
    cout << "Energy: " << energy << endl;
     }
}
