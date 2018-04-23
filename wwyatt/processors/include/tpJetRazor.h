#ifndef tpJetRazor_h
#define tpJetRazor_h 1
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
using namespace std; 
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
 * @version $Id: tpJetRazor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class tpJetRazor : public Processor {

    public:

        virtual Processor*  newProcessor() { return new tpJetRazor ; }


        tpJetRazor() ;

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
        vector<vector<PseudoJet>> getMegajets(vector<PseudoJet> jets); 
        vector<vector<double>> boostMegajets(vector<double> j1, vector<double> j2);  
        std::map<int,std::string> pname;

    protected:

        /** Input collection name.
        */
        std::string _colName ;
        std::string _root_file_name;
        std::string j0eventsCheck;
        std::string j1eventsCheck; 

        int MR0check; 
        
        std::string betaCheck;
        std::string Rcheck;   

        int _jetDetectability; 
        double _JetRParameter = 1.5;
     
        int _boost;  
       
        int _cuts[3]; 
        LCCollection* _inParVec;
        std::vector<PseudoJet> _parp;
        std::string filename;
        RanluxEngine myrnd;

        int _nRun ;
        int _nEvt ;
};

#endif



