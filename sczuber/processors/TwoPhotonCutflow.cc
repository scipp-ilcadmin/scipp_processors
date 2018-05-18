#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * Ok, so I like C++11. Unfortunately,
 * Marlin is built with ansi C, so the processor
 * constructor freaks out about the string that is
 * passed to it as an argument. The above two lines
 * fix that issue, allowing our code to be compatible
 * with ansi C class declarations.
 * Big thanks to Daniel Bittman for helping me fix this.
 */

/*
 * author Christopher Milke
 * April 5, 2016
 *
 * author Summer Zuber 
 * March 4, 2017 
 *
 * In order to analyze the kinematic event oversables 
 * 'S', 'V', and 'M' of the hadronic system of two-photon events 
 * (The hadronic system being all final state particles less the initial 
 * scattered electron and positron)
 * S = sum of the magnitudes of the transerve momentum
 * of each particle 
 * V = magnitude of the vector sum of teh transverse mometum 
 * of each particle 
 * M = sqrt((sum(energies))^(2)-((sum(px))^(2)+(sum(py))^(2)+(sum(pz))^(2)) )
 * (mass of event) 
 *
 * We analyze these at three different levels of "detectability"
 * True: entire hadronic system
 * Detectable: True less neutrinos 
 * Detected: Detectable less forward particle with |cos(theta)| > 0.9 
 *
 */

#include "TwoPhotonCutflow.h"
#include "scipp_ilc_utilities.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <iostream>
#include <fstream>
#include <vector>

// #include <CLHEP/Vector/ThreeVector.h>
// #include <CLHEP/Random/RanluxEngine.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCIO.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCTOOLS.h>


using namespace lcio;
using namespace marlin;
using namespace std;

TwoPhotonCutflow TwoPhotonCutflow;

TwoPhotonCutflow::TwoPhotonCutflow() : Processor("TwoPhotonCutflow") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
}

void TwoPhotonCutflow::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ; 
    _nEvt = 0 ;
}



void TwoPhotonCutflow::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ; 

} 

void TwoPhotonCutflow::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

    _inParVec = evt->getCollection( _colName) ; 


    cout << _colName << endl;
    cout << "event = " << _nEvt << endl;

    if( _colName != NULL ){
        double vec[4][3];
        double scalars[4];
        double energy[4];

        int id, stat;

        int nElements = _inParVec->getNumberOfElements()  ;
        cout << _inParVec->getNumberOfElements() << endl;

        // For each particle in Event ...
        for(int particleIndex = 0; particleIndex < nElements ; particleIndex++){
            MCParticle* particle = dynamic_cast<MCParticle*>( _inParVec->getElementAt(particleIndex));

            try{ id = particle->getPDG();
                //cout << id << endl; 
            }
            catch(...){cout << "could not get particle id" << endl;}

            stat = particle->getGeneratorStatus();
            // If Particle is FINAL-STATE 
            if(stat==1){
                bool isDarkMatter = (id == 1000022);
                if(isDarkMatter) continue ;
                double E = particle->getEnergy();
                const double* P = particle->getMomentum();
                double px = P[0];
                double py = P[1];
                double pz = P[2];
                double Pmag = sqrt(px*px+py*py+pz*pz);
                double cos = pz/Pmag;
                double scalar = sqrt(px*px+py*py);
                bool isNeutrino = (

                        id == 12 || id == -12 ||
                        id == 14 || id == -14 ||
                        id == 16 || id == -16 ||
                        id == 18 || id == -18);
                bool isForward = ( cos > 0.9 || cos < -0.9);
                scalars[0]+=scalar;
                vec[0][0]+=px;
                vec[0][1]+=py;
                vec[0][2]+=pz;
                energy[0]+=E;
                if(!isDarkMatter && !isNeutrino){
                    scalars[2]+=scalar;
                    vec[2][0]+=px;
                    vec[2][1]+=py;
                    vec[2][2]+=pz;
                    energy[2]+=E;
                    if(!isForward){
                        scalars[1]+=scalar;
                        vec[1][0]+=px;
                        vec[1][1]+=py;
                        vec[1][2]+=pz;
                    }
                }

            }//end final state
        }//end for

        //all
        double total_true_scalar = scalars[0];
        double total_detected_scalar = scalars[1];
        double total_detectable_scalar = scalars[2];

        double total_true_vector = sqrt(vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]);
        double total_detected_vector = sqrt(vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]);
        double total_detectable_vector = sqrt(vec[2][0]*vec[2][0]+vec[2][1]*vec[2][1]);


        double total_true_mass_squared = energy[0]*energy[0]-
            (vec[0][0]*vec[0][0]+vec[0][1]*vec[0][1]+
             vec[0][2]*vec[0][2]);
        double total_true_mass = sqrt(total_true_mass_squared);
        double total_detected_mass_squared = energy[1]*energy[1]-
            (vec[1][0]*vec[1][0]+vec[1][1]*vec[1][1]+
             vec[1][2]*vec[1][2]);
        double total_detected_mass = sqrt(total_detected_mass_squared);
        double total_detectable_mass_squared = energy[2]*energy[2]-
            (vec[2][0]*vec[2][0]+vec[2][1]*vec[2][1]+
             vec[2][2]*vec[2][2]);
        double total_detectable_mass = sqrt(total_detectable_mass_squared);
        cuts[0][0]+=1;
        cuts[1][0]+=1;
        cuts[2][0]+=1;
        if(total_true_scalar > 0.5){
            cuts[0][1]+=1;
            if(total_true_mass > 0.5){
                cuts[0][2]+=1;
                if(total_true_scalar > 1){
                    cuts[0][3]+=1;
                    if(total_true_mass > 1){
                        cuts[0][4]+=1;
                    }
                }
            }
        }

        if(total_detected_scalar > 0.5){
            cuts[1][1]+=1;
            if(total_detected_mass > 0.5){
                cuts[1][2]+=1;
                if(total_detected_scalar > 1){
                    cuts[1][3]+=1;
                    if(total_detected_mass > 1){
                        cuts[1][4]+=1;
                    }
                }
            }
        }

        if(total_detectable_scalar > 0.5){
            cuts[2][1]+=1;
            if(total_detectable_mass > 0.5){
                cuts[2][2]+=1;
                if(total_detectable_scalar > 1){
                    cuts[2][3]+=1;
                    if(total_detectable_mass > 1){
                        cuts[2][4]+=1;
                    }
                }
            }
        }


        /*cout << "TRUE" << endl;
          cout<< "cut_0 " << cuts[0][0] << endl;
          cout<< "cut_1 " << cuts[0][1] << endl;
          cout<< "cut_2 " << cuts[0][2] << endl;
          cout<< "cut_3 " << cuts[0][3] << endl;
          cout<< "cut_4 " << cuts[0][4] << endl; 
          cout << "DETECTABLE" << endl;
          cout<< "cut_0 " << cuts[2][0] << endl;
          cout<< "cut_1 " << cuts[2][1] << endl;
          cout<< "cut_2 " << cuts[2][2] << endl;
          cout<< "cut_3 " << cuts[2][3] << endl;
          cout<< "cut_4 " << cuts[2][4] << endl;
          cout << "DETECTED" << endl;
          cout<< "cut_0 " << cuts[1][0] << endl;
          cout<< "cut_1 " << cuts[1][1] << endl;
          cout<< "cut_2 " << cuts[1][2] << endl;
          cout<< "cut_3 " << cuts[1][3] << endl;
          cout<< "cut_4 " << cuts[1][4] << endl;*/
        double cut = 0;
        int i = 0;
        while(cut < 2.1){
            cout << cut << endl;
            if(total_true_scalar>cut){
                _cuts_sep[0][0][i]++;
            }
            if(total_true_vector>cut){
                _cuts_sep[0][1][i]++;
            }
            if(total_true_mass>cut){
                _cuts_sep[0][2][i]++;
            }
            if(total_detectable_scalar>cut){
                _cuts_sep[1][0][i]++;
            }
            if(total_detectable_vector>cut){
                _cuts_sep[1][1][i]++;
            }
            if(total_detectable_mass>cut){
                _cuts_sep[1][2][i]++;
            }
            if(total_detected_scalar>cut){
                _cuts_sep[2][0][i]++;
            }
            if(total_detected_vector>cut){
                _cuts_sep[2][1][i]++;
            }
            if(total_detected_mass>cut){
                _cuts_sep[2][2][i]++;
            }
            cut += 0.2;
            i +=1;

        }
    }
    cout << "end Event "<<_nEvt <<endl; 
    _nEvt ++ ; // different from original-moved out of for loop - summer                     
}  // end process Event 


void TwoPhotonCutflow::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}
void TwoPhotonCutflow::end(){  
    for (int i = 0; i<3; i++){
        cout << "i = "<< i <<"TRUE/DAB/DED"<< endl;
        for(int j = 0; j < 3; j++){
            cout <<" j = " << j << " S/V/M"<< endl;
            for (int c = 0; c <11; c++){

                cout << _cuts_sep[i][j][c] << endl;
            }
        }
    }
}

