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
#include <TH2D.h>


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
        virtual void BinLogX(TH2 * h);     
        //virtual void BinLogY(TH2 * h);     
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
        int mr0check; 
        std::string parpCheck;
        std::string betaCheck;
        int totalUnph; 
        std::string Rcheck; 
        std::string Rvals; 
        int totalJets; 
 
        int _jetDetectability; 
        // needed to put jet parameter at 1.5 because twoPhotonJetRazor was crashing 
        double _JetRParameter=0.5; // don't edit // lol 
        
        int _boost; // the type of R-frame transformation to do   
 
        int _cuts[6]; 
 
        int _concut;  
        int _concut_green;  
        int _concut_blue;  
        int _concut_yellow;  
        LCCollection* _inParVec; 
        std::vector<PseudoJet> _parp;
        std::string filename;
        RanluxEngine myrnd;

        int _nRun ;
        int _nEvt ;
};

#endif



