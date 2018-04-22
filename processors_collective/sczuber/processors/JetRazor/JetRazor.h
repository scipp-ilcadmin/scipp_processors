#ifndef JetRazor_h
#define JetRazor_h 1
#include <vector>
#include "marlin/Processor.h"
#include "lcio.h"
#include <iostream>
#include <string>
#include <IMPL/ReconstructedParticleImpl.h>

#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Random/RanluxEngine.h>

namespace CLHEP{}    // declare namespace CLHEP for backward compatibility
using namespace CLHEP ;
using namespace lcio ;
using namespace marlin ;


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
 * @version $Id: JetRazor.h,v 1.4 2005-10-11 12:57:39 gaede Exp $ 
 */

class JetRazor : public Processor {

    public:

        virtual Processor*  newProcessor() { return new JetRazor ; }


        JetRazor() ;

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


    protected:

        /** Input collection name.
        */
        std::string _colName ;
        std::string _root_file_name;
        int TassoJetRazor();
        int JetsetJetRazor();
        double sign(double a,double b);
        double min(double a,double b);

      /** Input collection name.
       */
        std::string partMomCheck;
        std::string betaCheck;

        bool partMom1;
        bool partMom0; 

        int _typeOfJetRazorFinder;
        int _thrustDetectability;
        float _principleJetRazorValue;
        float _majorJetRazorValue;
        float _minorJetRazorValue;
        Hep3Vector _principleJetRazorAxis;
        Hep3Vector _majorJetRazorAxis;
        Hep3Vector _minorJetRazorAxis;
        float _min,_max;
        LCCollection* _inParVec;
        std::vector<Hep3Vector> _partMom;
        std::string filename;
        RanluxEngine myrnd;

        int _nRun ;
        int _nEvt ;
};

#endif



