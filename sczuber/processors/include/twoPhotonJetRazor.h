#ifndef twoPhotonJetRazor_h
#define twoPhotonJetRazor_h 1
#include <vector>
#include "marlin/Processor.h"
#include "lcio.h"
#include <iostream>
#include <string>
#include <IMPL/ReconstructedParticleImpl.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Random/RanluxEngine.h>
#include <fastjet/ClusterSequence.hh>

namespace CLHEP{}    // declare namespace CLHEP for backward compatibility
using namespace CLHEP ;
using namespace lcio ;
using namespace marlin ;
using namespace fastjet;

/**  Example processor for marlin.
 * 
 *  If compiled with MARLIN_USE_AIDA 
 *  it creates a histogram (cloud) of the MCParticle energies.
 * 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of MCParticles.
 *
 *  <h4>Output</h4> 
 *  A histogram.
 * 
 * @param CollectionName Name of the MCParticle collection
 * 
 * @author F. Gaede, DESY
 * @version $Id: twoPhotonJetRazor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class twoPhotonJetRazor : public Processor {

    public:

        virtual Processor*  newProcessor() { return new twoPhotonJetRazor ; }


        twoPhotonJetRazor() ;

        /** Called at the begin of the job before anything is read.
         * Use to initialize the processor, e.g. book histograms.
         */
        virtual void init() ;

        /** Called for every run.
        */
        virtual void processRunHeader( LCRunHeader* run ) ;
        virtual void modifyRunHeader( LCRunHeader* run ) {}
        /** Called for every event - the working horse.
        */
        virtual void processEvent( LCEvent * evt ) ; 

        virtual void check( LCEvent * evt ) ; 


        /** Called after data processing for clean up.
        */
        virtual void end() ;
        std::map<int,std::string> pname;

    protected:

        /** Input collection name.
        */
        std::string _colName ;
        std::string _root_file_name;
        std::string j0eventsCheck;
        std::string j1eventsCheck; 

        std::string betaCheck; 
        //bool partMom_0;
        //bool partMom_1; 

      /** Input collection name.
       */
        int _jetDetectability; 
        double _JetRParameter = 1.5;
      
        LCCollection* _inParVec;
        std::vector<PseudoJet> _parp;
        std::string filename;
        RanluxEngine myrnd;

        int _nRun ;
        int _nEvt ;
};

#endif



