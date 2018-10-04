#ifndef final_h
#define final_h 1

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
class final : public Processor {

 public:

  virtual Processor*  newProcessor() { return new final ; }


  final() ;

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
