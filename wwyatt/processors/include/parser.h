#ifndef parser_h
#define parser_h 1

#include "marlin/Processor.h"

#include "lcio.h"
#include <string>
#include <vector>
#include <algorithm>
#include "scipp_ilc_utilities.h"
#include <TFile.h>
#include <TH2D.h>

using namespace std;
using namespace lcio ;
using namespace marlin ;


//Normal Processor Header
class parser : public Processor {

 public:
  typedef vector<vector<MCParticle*>*> Community;
  typedef vector<MCParticle*> Family;
  static void printParticle(MCParticle *p);
  static void printAllEvents(LCCollection* col);


  virtual Processor*  newProcessor() { return new parser ; }


  parser() ;

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

  bool sameTree(const Family*,const Family*);
  bool sameTree(const Family&,const Family&);
  bool sameMomentum(const double*,const double*);
  bool sameParticle(const MCParticle*, const MCParticle*, bool=true);
  Community* removeDuplicateTrees(Community*);
  Family* removeDuplicateParticles(Family*);
  enum class Direction {children, parents};

  Family* removeTraversed(Family*, Family*);
  Family* traverse(MCParticle*, Family* =NULL);
  bool inFamily(MCParticle*, Family*);
  Family biggerTree(Family,Family);

  void addToTree(Family*, Community* ,Family*);
  Community* getTrees(LCEvent*, bool=false);
  Family* getAssociates(MCParticle*);
  unsigned int countBhabhas(Community*);
 protected:

  /** Input collection name.
   */
  std::string _colName ;
  std::string _root_file_name;

  int _nRun ;
  int _nEvt ;
};

#endif



