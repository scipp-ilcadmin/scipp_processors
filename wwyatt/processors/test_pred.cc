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

#include "test_pred.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <iomanip> 

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <MyParticle.h>
#include <Will.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;
using namespace Will;
using namespace TwoPhoton;


test_pred test_pred;


static TFile* _rootfile;
static TH2F* _prediction;
static TH2F* _observe;
static TH1F* _vector;
static TH1F* _p_theta;
static TH1F* _e_theta;
static TH1F* zmom;
static TH1F* tmom;
static TH1F* amom;
static TH1F* bmom;
static TH1F* cmom;
static TH1F* dmom;

static vector<Result> positron_results;
static vector<Result> electron_results;
static vector<fourvec> actual;
static vector<fourvec> predicted;
static vector<double> spread_e;
static vector<double> spread_p;
static vector<pair<double,double>> spread;

static int p_scatter=0;
static int e_scatter=0;
static int ph_th = 0;
static int ph_tm = 0;
static int pm_th = 0;
static int pm_tm = 0;
static int total = 0;

static int photons=0;
static int nphotons=0;

static int n_events=0;
static int i_events=0;

test_pred::test_pred() : Processor("test_pred") {
  // modify processor description
  _description = "Protype Processor" ;

  // register steering parameters: name, description, class-variable, default value
  registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

  registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void test_pred::init() { 
  streamlog_out(DEBUG) << "   init called  " << std::endl ;

  _rootfile = new TFile("ppredict.root","RECREATE");
  //usually a good idea to
  //printParameters() ;
  _prediction = new TH2F("predict", "Predicted Angle of Scatter, Correct vs Incorrect Kinematics", 1000, 0.0, 0.01, 1000, 0.0, 0.01);
  _observe = new TH2F("energy", "B vs System Energy", 1000, 0.0, 1.5, 1000, 0.0, 525);
  _p_theta = new TH1F("p_theta", "Theta between positron and hadronic system", 360, 0, .1);
  _e_theta = new TH1F("e_theta", "Theta between positron and hadronic system", 360, 0, .1);
  _vector = new TH1F("vector", "Vector", 200, 0.0, 0.05);
  zmom=new TH1F("zmom", "System energy", 500, 450, 550);
  tmom=new TH1F("tmom", "Eletron system energy", 500, 450, 550);
  amom=new TH1F("amom", "Distribution of Electron Theta", 500, 0,.1);
  bmom=new TH1F("bmom", "Distribution of Positron Theta", 500, 0,.1);
  cmom=new TH1F("cmom", "Distribution of Electron Theta", 500, 0,.1);
  dmom=new TH1F("dmom", "Distribution of Positron Theta", 500, 0,.1);
  _nEvt = 0 ;
}



void test_pred::processRunHeader( LCRunHeader* run) { 
  //    _nRun++ ;
} 

void test_pred::processEvent( LCEvent * evt ) { 
  return;
}

void test_pred::check( LCEvent * evt ) { 
  // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void test_pred::end(){ 
  //Run a test for each function in TwoPhoton
  //transform_to_lab(fourvec)
  
  //get_hitStatus(fourvec)
  //prediction(bundle)
  //getHadronicSystem(LCCollection* col);
  //getBeamcalPosition(fourvec);
  //getHMGrid(vector<fourvec> predicted, vector<fourvec> actual);
  //getHMGrid(vector<Result> input, double energy_cut);
  //recordHMValue(hmgrid &output, fourvec predicted, fourvec actual){
  _rootfile->Write();
}
