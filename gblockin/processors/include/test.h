#ifndef test_h
#define test_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

using namespace lcio ;
using namespace marlin ;

class test : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new test ;  }
  
  test() ;

  virtual void init() ;
  
  virtual void processRunHeader( LCRunHeader* run  ) ;

  virtual void processEvent( LCEvent * evt  ) ;

  virtual void bcalLoop ( LCEvent * evt ) ;
  
  virtual void check( LCEvent * evt  ) ; 
  
  virtual void end() ;
  
  // declare your new methods here
 protected:

  std::string _colName;
  std::string _bcolName;
  // tell you which event is running 
  int _nRun ;
  int _nEvt ;

  // declare your new variables here
} ;

#endif


