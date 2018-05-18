#ifndef TwoPhotonCutflow_h
#define TwoPhotonCutflow_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>

using namespace lcio ;
using namespace marlin ;

class TwoPhotonCutflow : public Processor {

    public:

        virtual Processor*  newProcessor() { return new TwoPhotonCutflow ; }


        TwoPhotonCutflow() ;

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

      /** Input collection name.
       */ 
        LCCollection* _inParVec;   

        int _nRun ;
        int _nEvt ;

        double _cuts_sep[3][3][11];
        double cuts[3][5];
};

#endif



