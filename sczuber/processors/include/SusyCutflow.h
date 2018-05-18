#ifndef SusyCutflow_h
#define SusyCutflow_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>


using namespace lcio ;
using namespace marlin ;

class SusyCutflow : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new SusyCutflow ; }
  
  
  SusyCutflow() ;
  
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

    int _nRun ;
    int _nEvt ;


    double _cuts_sep[3][3][11];
    double cuts[3][5];
} ;

#endif



