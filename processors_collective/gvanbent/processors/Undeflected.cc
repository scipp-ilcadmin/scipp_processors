#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * C. Milke: 
 * Ok, so I like C++11. Unfortunately,
 * Marlin is built with ansi C, so the processor
 * constructor freaks out about the string that is
 * passed to it as an argument. The above two lines
 * fix that issue, allowing our code to be compatible
 * with ansi C class declarations.
 * Big thanks to Daniel Bittman for helping me fix this.
 */

/*
 * author Yevgeniya Shtalenkova
 * June 6, 2017
 */

#include "Undeflected.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <cmath>

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

#define _USE_MATH_DEFINES


using namespace lcio;
using namespace marlin;
using namespace std;
using namespace scipp_ilc;

Undeflected Undeflected;

static TFile* _rootfile;
static TH2D* _pos;

Undeflected::Undeflected() : Processor("Undeflected") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}



void Undeflected::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("BB_pos.root","RECREATE");
    _pos = new TH2D("pos", "Undeflected Beam Particle Position", 600, -150.0, 150.0, 600, -150.0, 150.0);    
    
    // usually a good idea to
    //printParameters() ;
    _nEvt = 0 ;

    _tot = 0;
}



void Undeflected::processRunHeader( LCRunHeader* run) { 
//    _nRun++ ;
} 



void Undeflected::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...
    double tot[]={0, 0};
    LCCollection* col = evt->getCollection( _colName ) ;

    int stat, id = 0;
    double theta_e, theta_p;
    double pT_e, pT_p;
    double E_e, E_p = 0;
    double mom_e[4], mom_p[4];
    double pred_e[3], pred_p[3];
    double pos_e[3], pos_p[3];

    double hadronic[4];
    double electronic[4];
     
    double diff;
     
    int hit = 0;
    //total scalar magnitude of hadronic system
    double S = 0;
    //total vector magnitude of system
    double V;

    // this will only be entered if the collection is available
    if( col != NULL ){
        int nElements = col->getNumberOfElements()  ;
    cout << endl;
    cout << endl;
    cout << "***************************EVENT: " << _nEvt << "****************************" << endl;

    vector<MCParticle*> system;
        //create final state subsystem
        //determine beam particle energies for identification
        for(int hitIndex = 0; hitIndex < nElements ; hitIndex++){
           MCParticle* hit = dynamic_cast<MCParticle*>( col->getElementAt(hitIndex) );
        
            id = hit->getPDG();
            stat = hit->getGeneratorStatus();
            if(stat==1){
                //fill final state system 
                system.push_back(hit);
                //identify highest e/p energies
                if(id==11){E_e = (E_e < hit->getEnergy()) ? hit->getEnergy() : E_e;}
                if(id==-11){E_p = (E_p < hit->getEnergy()) ? hit->getEnergy() : E_p;}
                //cout << "Particle " << hitIndex << " with ID: " << id;
                //cout << " status: " << stat;

                double mom[4];
                mom[0] = hit->getMomentum()[0]; 
                mom[1] = hit->getMomentum()[1]; 
                mom[2] = hit->getMomentum()[2];
                mom[3] = hit->getEnergy();

                //x and y mom sums for V
                tot[0]+=mom[0];
                tot[1]+=mom[1];    
                //cout << " momentum [" << mom[0] << ", " << mom[1] << ", " << mom[2] << "] with energy: " << mom[3] << endl;

            }//end final state
        }//end for
        
        //calculate V
        V=sqrt(pow(tot[0], 2)+pow(tot[1], 2));
        //create pseudo particle to balance momentum
        double pseudo_x = -tot[0];
        double pseudo_y = -tot[1];

        //Now only for final state system
        for(MCParticle* particle : system){
            id = particle->getPDG();
            
            //find beam electron
            if(id==11&&particle->getEnergy()==E_e){
                    const double* mom = particle->getMomentum();
                    mom_e[0]=mom[0];
                    mom_e[1]=mom[1];
                    mom_e[2]=mom[2];
                    mom_e[3]=particle->getEnergy();
                    //if electron is deflected, calculate angle 
                    if(abs(mom_e[0])!=0||abs(mom_e[1])!=0){
                        hit++;
                        electronic[0]=mom_e[0];
                        electronic[1]=mom_e[1];
                        electronic[2]=mom_e[2];
                        electronic[3]=mom_e[3];
                        pT_e = sqrt(pow(mom_e[0], 2)+pow(mom_e[1], 2));
                        double mag = sqrt(pow(mom_e[0], 2)+pow(mom_e[1], 2)+pow(mom_e[2], 2));
                        theta_e = asin(pT_e/mag);
                    }//end deflection
                    else{
                        scipp_ilc::transform_to_lab(mom_e[0], mom_e[3], mom_e[0], mom_e[3]);    
                        double pos[3];
                        pos[2] = 3265;
                        pos[0] = mom_e[0]*pos[2]/mom_e[2];
                        pos[1] = mom_e[1]*pos[2]/mom_e[2];
                        pos[0] = pos[0] - pos[2]*0.007;
                        _pos->Fill(pos[0], pos[1]);                       
                        if(abs(pos[0])>0 || abs(pos[1])>0){_tot++; cout << "ELECTRON" << endl; cout << "X: " << pos[0] << endl; cout << "Y: " << pos[1] << endl;}
                    }
            }//end beam electron
              
            //find beam positron
            else if(id==-11&&particle->getEnergy()==E_p){
                    const double* mom = particle->getMomentum();
                    mom_p[0] = mom[0]; 
                    mom_p[1] = mom[1]; 
                    mom_p[2] = mom[2];
                    mom_p[3] = particle->getEnergy();
                    //if electron is deflected, calculate angle
                    if(abs(mom_p[0])!=0||abs(mom_p[1])!=0){
                        hit++;
                        electronic[0]=mom_p[0];
                        electronic[1]=mom_p[1];
                        electronic[2]=mom_p[2];
                        electronic[3]=mom_p[3];
                        pT_p = sqrt(pow(mom_p[0], 2)+pow(mom_p[1], 2));
                        double mag = sqrt(pow(mom_p[0], 2)+pow(mom_p[1], 2)+pow(mom_p[2], 2));
                        theta_p = asin(pT_p/mag);
                    }//end delfection
                    else{
                        scipp_ilc::transform_to_lab(mom_p[0], mom_p[3], mom_p[0], mom_p[3]);    
                        double pos[3];
                        pos[2] = -3265;
                        pos[0] = mom_p[0]*pos[2]/mom_p[2];
                        pos[1] = mom_p[1]*pos[2]/mom_p[2];
                        pos[0] = pos[0] + pos[2]*0.007;
                        _pos->Fill(pos[0], pos[1]);                       
                        if(abs(pos[0])>0 || abs(pos[1])>0){_tot++; cout << "POSITRON" << endl; cout << "X: " << pos[0] << endl; cout << "Y: " << pos[1] << endl;}
                    
                    }
            }//end beam positron    
            
            else{
                double mom[4];
                mom[0] = particle->getMomentum()[0]; 
                mom[1] = particle->getMomentum()[1]; 
                mom[2] = particle->getMomentum()[2];
                mom[3] = particle->getEnergy();

                double hyp = sqrt(pow(mom[0], 2)+pow(mom[1], 2)+pow(mom[2], 2));
                double cos = mom[2]/hyp;

                //if(cos<0.9){
                    hadronic[0]+=mom[0];
                    hadronic[1]+=mom[1];
                    hadronic[2]+=mom[2];
                    hadronic[3]+=mom[3];
                //}
                double mag = sqrt(pow(mom[0], 2)+pow(mom[1], 2));
                S+=mag;    
            }
        }//end for


        
        //add psuedo particle
        hadronic[0]+=pseudo_x;
        hadronic[1]+=pseudo_y;


    }//end collection



    _nEvt ++ ;
}



void Undeflected::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Undeflected::end(){ 
    cout << _tot << endl;
    _rootfile->Write();
}
