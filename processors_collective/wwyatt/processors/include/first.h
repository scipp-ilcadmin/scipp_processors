#ifndef first_h
#define first_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include <string>
#include <vector>
#include "scipp_ilc_utilities.h"
#include <TFile.h>
#include <TH2D.h>

using namespace std;
using namespace lcio ;
using namespace marlin ;

//Normal Processor Header
class first : public Processor {

 public:

  struct AngleStore{
    vector<double> hh_colinearity;
    vector<double> mh_colinearity;
    vector<double> hm_colinearity;
    vector<double> mm_colinearity;
  };
  struct HitMissStore{
    int hit = 0;
    int miss = 0;
    vector<double*> hh;
    vector<double*> hm;
    vector<double*> mh;
    vector<double*> mm;
    vector<double> hh_theta;
    vector<double> hm_theta;
    vector<double> mh_theta;
    vector<double> mm_theta;
  };
  //Verbose is mostly unused.
  const bool VERBOSE = false;
  //Setting PDG constants to particle names
  static  const int PHOTON = 22;
  static  const int POSITRON = 11;
  static  const int ELECTRON = -11;

  //Setting errors static so I can spell check them easily. Mostly unused, I though Bhabhas were going do be a larger project.
  const string ERROR_NOT_BAHBAH_PARTICLE = "Not a Bahbah particle User-Error: Trying to add a particle that is not accounted for.";
  const string ERROR_MAX_PHOTONS = "Max Photon User-Error: Trying to add photon but it is already full.";
  const string ERROR_ALREADY_PARTICLE = "Particle Already There User-Error: Trying to add a paticle but the space is not NULL.";

  bool hasInitialized = false;

  int info = 0;

  MCParticle* photonA = NULL;
  MCParticle* photonB = NULL;
  MCParticle* positron = NULL;
  MCParticle* electron = NULL;
  MCParticle* photons[2] = {NULL, NULL};

  HitMissStore positronStore;
  HitMissStore electronStore;
  AngleStore colinearityStore;
  //Do a check to see if it is ready, then orchestrates other functions.
  void initialize();

  //Calculate and process hh, hm, mh, mm.
  void init_hitmap(bool=true);

  //Used to calculate and plot the colinearity
  double getColinearity(bool=true);

  //Same thing as getColinearity but intead of using the particles from the parent bundle class it takes in two double arrays representing momentums.
  double getOpenAngle(double*, double*);

  //Python-like printing function, because I love python.
  void print(string _input);

  //This function is written to plot hit statuses
  void graphHitStatus(const double*, int id, bool=true);

  //Returns the total number of particles added. MAX is 4 right now.
  int getCount();

  //Checks to see if particle is already there, then adds it for processing. Once full it runs initialize.
  void addParticle(MCParticle* _input);

  //Preparing for updating class for n photons. Currently max is 2 photons.
  void addPhoton(MCParticle* _input);

  //Returns associated particle.
  MCParticle* getElectron();
  MCParticle* getPositron();

  //-- Physics Section --\\
  //  void plotHitMiss(TH2F*, vector<double*>);
  void plotHitMiss(TH2F*, initializer_list<vector<double*>>);

  void plotHistogram(TH1F*, initializer_list<vector<double>>);

  //returns the angle from the x-y axis
  double getPhi(MCParticle* _input, bool=true);

  //returns the angle from the z axis
  double getTheta(MCParticle*, bool=true);

  //returns the cosTheta between the two vectors.
  double getCosTheta(MCParticle*, bool=true);

  //Gets the norm of the vector and finds it's manitude.
  double getMagnitude(MCParticle* _input, bool=true);
  double getMagnitude(double*);

  //Retuns the magnitude of the dot-product.
  double getDotProduct(double*, double*);
  double getDotProduct(MCParticle*, MCParticle*, bool=true);

  //returns the lab-momentum of the particle as an array => {0: x-axis, 1: y-axis, 2: z-axis}
  double* getMomentum(MCParticle* _input, bool=true);


  virtual Processor*  newProcessor() { return new first ; }


  first() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 


  virtual void check( LCEvent * evt ) ; 


  /** Called after data processing for clean up.
   */
  virtual void end() ;


 protected:

  /** Input collection name.
   */
  std::string _colName ;
  std::string _root_file_name;

  int _nRun ;
  int _nEvt ;
};

#endif



