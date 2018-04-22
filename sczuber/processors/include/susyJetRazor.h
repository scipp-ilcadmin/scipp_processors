#ifndef susyJetRazor_h
#define susyJetRazor_h 1
#include <vector>
#include "marlin/Processor.h"
#include "lcio.h"
#include <iostream>
#include <string>
#include <IMPL/ReconstructedParticleImpl.h>
#include "fastjet/ClusterSequence.hh"
#include <map>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Random/RanluxEngine.h>

namespace CLHEP{}    // declare namespace CLHEP for backward compatibility
using namespace std; 
using namespace CLHEP ;
using namespace lcio ;
using namespace marlin ;
using namespace fastjet ;


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
 * @version $Id: susyJetRazor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class susyJetRazor : public Processor {

    public:

        virtual Processor*  newProcessor() { return new susyJetRazor ; }

        susyJetRazor() ;

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

        //virtual std::vector<PseudoJet> getMegajets(std::vector<PseudoJet> jets);
        //std::vector<PseudoJet> getMegajets(std::vector<PseudoJet> jets);
        //vector<PseudoJet> getMegajets(vector<PseudoJet> jets);
        vector<vector<PseudoJet>> getMegajets(vector<PseudoJet> jets);
        vector<vector<double>> boostMegajets(vector<double> j1, vector<double> j2); 
        //vector<int> doCuts(double MR, double R2); 
        
        std::map<int,std::string> pname;

    protected:

        /** Input collection name.
        */
        std::string _colName ;
        std::string _root_file_name; 

      /** Input collection name.
       */
        std::string j0eventsCheck;
        std::string j1eventsCheck;
        int numj0eventsCheck; 
        int numj1eventsCheck;
        std::string parpCheck;
        std::string betaCheck;
        std::string Rcheck; 
 
        int _jetDetectability; 
        // needed to put jet parameter at 1.5 because twoPhotonJetRazor was crashing 
        double _JetRParameter=1.5; // don't edit // lol 
        
        int _boost; // the type of R-frame transformation to do   
 
        int _cuts[3]; 
  
        LCCollection* _inParVec; 
        std::vector<PseudoJet> _parp;
        std::string filename;
        RanluxEngine myrnd;

        int _nRun ;
        int _nEvt ;
};

#endif



