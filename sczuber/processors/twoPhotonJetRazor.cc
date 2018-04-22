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
 * January 18, 2017 
 *
 * Using fastjet jet reconstruction algorithm  
 */

#include "twoPhotonJetRazor.h"
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
using namespace fastjet;

static TFile* _rootfile;


static TH1F* _beta_T;
static TH1F* _beta_DAB;
static TH1F* _beta_DED;

static TH1F* _MR_T;
static TH1F* _MR_DAB;
static TH1F* _MR_DED;

static TH1F* _MRT_T;
static TH1F* _MRT_DAB;
static TH1F* _MRT_DED;

static TH1F* _R_T;
static TH1F* _R_DAB;
static TH1F* _R_DED;

twoPhotonJetRazor twoPhotonJetRazor;

twoPhotonJetRazor::twoPhotonJetRazor() : Processor("twoPhotonJetRazor") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "jetDetectability" ,
            "Detectability Level particles used in the Jet reconstruction:\n#\t0 : True\n#\t1 : Detectable\n#\t2 : Detected" ,
            _jetDetectability, 2  );
}

void twoPhotonJetRazor::init() { 
    //streamlog_out(DEBUG)  << "   init called  " << endl;


    if(_jetDetectability==0){_rootfile = new TFile("twoPhotonJetRazor_eW.pW.I39212._T.root","RECREATE");
        _R_T = new TH1F("R_T", "R=MTR/MR",250,0,50);
        _MR_T = new TH1F("MR_T", "MR",250,0,50);
        _MRT_T = new TH1F("MRT_T", "MRT",250,0,50);
        _beta_T = new TH1F("beta_T", "beta",250,0,50);
        
        freopen("twoPhotonJetRazor_eW.pW.I39212._T.log", "w",stdout);    
    }
    if(_jetDetectability==1){_rootfile = new TFile("twoPhotonJetRazor_.eW.pW.I39212._DAB.root","RECREATE");
        _R_DAB = new TH1F("R_DAB", "R=MTR/MR",250,0,50);
        _MR_DAB = new TH1F("MR_DAB", "MR",250,0,50);
        _MRT_DAB = new TH1F("MRT_DAB", "MRT",250,0,50);
        _beta_DAB = new TH1F("beta_DAB", "beta",130,-3,10);
        freopen("twoPhotonJetRazor_.eW.pW.I39212._DAB.log", "w",stdout);   
    }
    if(_jetDetectability==2){_rootfile = new TFile("twoPhotonJetRazor_eW.pW.I39212._DED.root","RECREATE");
        _R_DED = new TH1F("R_DED", "R=MTR/MR",250,0,50); 
        _MR_DED = new TH1F("MR_DED", "MR",250,0,50); 
        _MRT_DED = new TH1F("MRT_DED", "MRT",250,0,50); 
        _beta_DED = new TH1F("beta_DED", "beta",130,-3,10); 
        freopen("twoPhotonJetRazor_eW.pW.I39212._DED.log", "w",stdout);    
    }

    //irameters() ;

    // config ranlux 
    filename = "Ranlux.coonf";
    ifstream rndcfgfile( filename.c_str() );

    if (!rndcfgfile)
    {
        long int ss=1234;
        myrnd.setSeeds(&ss,4);
        myrnd.showStatus();
    }
    else
    {
        rndcfgfile.close();
        myrnd.restoreStatus(filename.c_str());
        myrnd.showStatus();
    } // if file not existusually a good idea to
    //printParameters() ;
    _nEvt = 0 ;
    pname[11  ]="               electron |";
    pname[-11  ]="              positron |";
    pname[12  ]="      electron neutrino |";
    pname[-12 ]="anti electron neutrino |";
    pname[13  ]="                   muon |";
    pname[-13 ]="             anti muon |";
    pname[14  ]="          muon neutrino |";
    pname[-14 ]="    anti muon neutrino |";
    pname[15  ]="                   tau |";
    pname[16  ]="           tau neutrino |";
    pname[-16 ]="     anti tau neutrino |";
    pname[22  ]="                 photon |";
    pname[1000022]="               X10 |";
    pname[211 ]="                   pi+ |";
    pname[-211]="                  pi- |";
    pname[321]="                   K+ |";
    pname[-321]="                  K- |";
    pname[130]="                 K_L^0 |";
    pname[2112]="                    n |";
    pname[-2112]="                   -n |";
    pname[2212]="                   p |";
    pname[-2212]="                 -p |";
}



void twoPhotonJetRazor::processRunHeader( LCRunHeader* run) {  
    //    _nRun++ ; 

} 



void twoPhotonJetRazor::processEvent( LCEvent * evt ) { 
    
    cout << "EVENT: "<< _nEvt<< endl;   

    _inParVec = evt->getCollection( _colName) ;
    //cout << " get collection " << endl;
    //cout << _inParVec->getNumberOfElements() << endl;
    if (!_parp.empty()) _parp.clear(); 

    MCParticle* high_e;
    MCParticle* high_p; 

    int id, stat; 

    // setting high_e and high_p to a final state electron and positron
 
    for(int n=0; n<_inParVec->getNumberOfElements(); n++)
    {
        MCParticle* hit = dynamic_cast<MCParticle*>(_inParVec->getElementAt(n));
        id = hit->getPDG();
        stat = hit->getGeneratorStatus();
        if(stat==1){
            if(id==11){
                high_e = hit;
            }
            if(id==-11){
                high_p = hit;
            }
        }// end final state 
    } // end for loop

    for (int n = 0; n<_inParVec->getNumberOfElements(); n++)
    {
        MCParticle* hit = dynamic_cast<MCParticle*>(_inParVec->getElementAt(n));
        id = hit->getPDG();
        stat = hit->getGeneratorStatus();
        if(stat==1){
            if(id==11){
                if(hit->getEnergy()>high_e->getEnergy()){
                    high_e = hit;
                }
            }
            if(id == -11){
                if(hit->getEnergy()>high_p->getEnergy()){
                    high_p = hit;
                }
            }
        }// end final state 
    }//end particle for loop 

    for (int n=0;n<_inParVec->getNumberOfElements() ;n++){
        MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );
        try{
            id = aPart->getPDG();
            stat = aPart->getGeneratorStatus();
        }
        catch(const std::exception& e){
            cout << "exception caught with message " << e.what() << "\n";
        }
        if(stat==1){
            const double* parp = aPart->getMomentum();
            double penergy= aPart->getEnergy();
            double pMag = sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]);
            double cos = parp[2]/pMag;
            bool isForward = (cos > 0.9 || cos < -0.9);
            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable && !isForward);
            // exclude highest energy electron and postron:
            if(aPart != high_e && aPart != high_p){ 
                if(_jetDetectability == 0){ 
                    _parp.push_back( PseudoJet(parp[0], parp[1], parp[2], penergy) ); 
                }
                if(_jetDetectability == 1){
                    if(isDetectable){
                        _parp.push_back(PseudoJet(parp[0], parp[1], parp[2], penergy));
                    }
                }
                if(_jetDetectability == 2){
                    if(isDetected){
                        _parp.push_back(PseudoJet(parp[0], parp[1], parp[2], penergy));
                    }
                }
            } // end particle not highest e+e-
        } //stat==1
    } // for particle 



    for (int n=0;n<_inParVec->getNumberOfElements() ;n++){
        MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );
        try{
            id = aPart->getPDG();
            stat = aPart->getGeneratorStatus();
        }
        catch(const std::exception& e){
            cout << "exception caught with message " << e.what() << "\n";
        }
        if(stat==1){
            // remove highest energy electron/posistron 
            if(aPart != high_e && aPart != high_p){
                const double* parp = aPart->getMomentum();
                double par4p[4] = {aPart->getEnergy(), parp[0], parp[1], parp[2]} ; 

                bool isDarkMatter = (id == 1000022);
                bool isNeutrino = (
                        id == 12 || id == -12 ||
                        id == 14 || id == -14 ||
                        id == 16 || id == -16 ||
                        id == 18 || id == -18);
                double cos = parp[2]/(sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]));
                bool isForward = ( cos > 0.9 || cos < - 0.9);
                bool isDetectable = (!isDarkMatter && !isNeutrino);
                bool isDetected = (isDetectable && !isForward);

                if(_jetDetectability == 0){
                    if(!isDarkMatter){
                        cout << "id: " << id << " "<< pname[id] << " "<< parp[0]<< " "<<parp[1]<< " "<< parp[2]<< endl;
                    }
                }
                if(_jetDetectability == 1){
                    if(isDetectable){
                        cout << "id: " << id << " "<< pname[id] << " "<< parp[0]<< " "<<parp[1]<< " "<< parp[2]<< endl;
                    }
                }
                if(_jetDetectability == 2){
                    if(isDetected){
                        cout << "id: " << id << " "<< pname[id] << " "<< parp[0]<< " "<<parp[1]<< " "<< parp[2]<< endl;
                    }
                }
            }// end particle not original 
        } //stat == 1 
    } // for particle

    // identify the jets using the fastjet clustering algorithm:
    JetDefinition jet_def(antikt_algorithm, _JetRParameter );
    // run the clustering, extract the jets
    ClusterSequence cs(_parp, jet_def);

    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());
    //print the number of jets in event:
    cout << "NUMBER OF JETS: "<< jets.size()<< endl;
    //check if 0 jets:
    if(jets.size()==0){
        j0eventsCheck +=" ";
        j0eventsCheck +=std::to_string(_nEvt);
        j0eventsCheck +=" ";

    }
    //check if 1 jet: 
    if(jets.size()==1){
        j1eventsCheck +=" ";
        j1eventsCheck +=std::to_string(_nEvt);
        j1eventsCheck +=" ";
    }

    // RAZOR MEGAJET ALGORITHM ------------------------------------------------
    int num_partitions = pow(2, jets.size());  
    double m2_values[num_partitions]; 

    int p = 0;
    while(p<num_partitions){
        m2_values[p]=0;
        p++;
    }
    for(int i = 0; i<num_partitions; i++){
        vector<PseudoJet> partition[2]; 
        for(int j = 0; j<jets.size(); j++){
            int bit = pow(2,j);
            if((i & bit) !=0){
                partition[0].push_back(jets[j]);
            }
        }
        for (int j = 0; j<jets.size(); j++){
            if(std::find(partition[0].begin(),partition[0].end(), jets[j]) != partition[0].end()){
                // partition[0] contains jets[j]
                //cout << "no" << endl; 
            }
            else{
                partition[1].push_back(jets[j]);
            }
        }
        double subset_p[2][4] = {{0,0,0,0},{0,0,0,0}};
        // for both subsets add up 4 momenta of each jet in subset 
        for(int p = 0; p<2; p++){
            for(int j = 0; j<partition[p].size(); j++){
                subset_p[p][0]=subset_p[p][0]+partition[p][j].e();
                subset_p[p][1]=subset_p[p][1]+partition[p][j].px();
                subset_p[p][2]=subset_p[p][2]+partition[p][j].py();
                subset_p[p][3]=subset_p[p][3]+partition[p][j].pz();
            }
        }
        double subset_m2[2]= {0,0};
        subset_m2[0]=pow(subset_p[0][0],2)-(pow(subset_p[0][1],2)+pow(subset_p[0][2],2)+pow(subset_p[0][3],2));
        subset_m2[1]=pow(subset_p[1][0],2)-(pow(subset_p[1][1],2)+pow(subset_p[1][2],2)+pow(subset_p[1][3],2));
        double partition_m2 = subset_m2[0]+subset_m2[1]; // sum of the two megajets 
        // print the m2 value for this partition: 
        m2_values[i]= partition_m2;
        //cout <<"m2: "<< m2_values[i] << endl ;
    } // num_partitions
    // set the min m2 value to the last valid 
    double min_value = m2_values[num_partitions-2]; // MIGHT BE PROBLEM IF *NO* JETS OR *ONE* JET!!
    int min_index = num_partitions-2;
    cout <<"initial min_value: "<< min_value<< endl;
    for(int i=0; i<num_partitions; i++){
        // if subset is not null set or set itself:
        if(i!=0 && i!=(num_partitions-1) ){
            if(m2_values[i]<min_value){
                min_index = i;
                min_value = m2_values[i];
            }
        }
    }
    cout <<"~~minimum index: "<< min_index<<" ~~minimum value: "<< min_value<< endl;
    // re-create the partition corresponding to the min_index: 
    vector<PseudoJet> megajet[2];
    for(unsigned int j=0; j<jets.size(); j++){
        unsigned int bit = pow(2,j);
        if ((min_index & bit) != 0){
            megajet[0].push_back(jets[j]);
        }
    }
    for(unsigned int j=0;j<jets.size(); j++){
        if(std::find(megajet[0].begin(),megajet[0].end(), jets[j]) !=megajet[0].end()){
            //subset contains jets[j]
        }
        else {
            //subset does not contain jets[j]
            megajet[1].push_back(jets[j]);
        }
    }
    // megajet is now the correct megajet to use

    // ----------------------------------------------------------------------------------------------
    //
    //     // CALCULATE THE MR, MRT, AND R VARIABLED:
    //         // the 4-momenta of the two megajets: 
    double j1[4]={0,0,0,0};
    double j2[4]={0,0,0,0};
    for(int i = 0; i<megajet[0].size(); i++){
        j1[0] = j1[0]+megajet[0][i].e();
        j1[1] = j1[1]+megajet[0][i].px();
        j1[2] = j1[2]+megajet[0][i].py();
        j1[3] = j1[3]+megajet[0][i].pz();
    }
    for(int i = 0; i<megajet[1].size(); i++){
        j2[0] = j2[0]+megajet[1][i].e();
        j2[1] = j2[1]+megajet[1][i].px();
        j2[2] = j2[2]+megajet[1][i].py();
        j2[3] = j2[3]+megajet[1][i].pz();
    }
    // the squared-invariant-mass of teh two megajets:
    double m2j1 = pow(j1[0],2)-pow(j1[1],2)-pow(j1[2],2)-pow(j1[3],2);
    double m2j2 = pow(j2[0],2)-pow(j2[1],2)-pow(j2[2],2)-pow(j2[3],2);
    // sum of the squared-invariant-mass values: 
    double m2 = m2j1+m2j2;

    cout <<"Megajet m2 value: " << endl;
    cout << "m2: "<< m2 << endl;

    // boost into R-frame:  ------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------
    double beta = (j1[0]-j2[0])/(j1[3]-j2[3]);
    cout << "beta: "<< beta<< endl;
    double gamma = pow((1-pow(beta,2)), -0.5);
    cout << "gamma: "<< gamma<< endl;
    double j1R[4]={gamma*j1[0]-gamma*beta*j1[3], j1[1],j1[2],-gamma*beta*j1[0]+gamma*j1[3] } ;
    double j2R[4]={gamma*j2[0]-gamma*beta*j2[3], j2[1],j2[2],-gamma*beta*j2[0]+gamma*j2[3] } ;

    double MRj1 = 2*((pow(j1R[1],2), pow(j1R[2],2), pow(j1R[3],2)));
    double MRj2 = 2*((pow(j2R[1],2), pow(j2R[2],2), pow(j2R[3],2)));
    cout << "MRj1: "<< MRj1 << endl;
    cout << "MRj2: "<< MRj2 << endl;

    if(beta>1){

        betaCheck += " ";
        betaCheck += std::to_string(_nEvt);
        betaCheck += " ";
    }
    //----------------------------------------------------------------------------------------------

    // The Razor Variables *without* using the R-frame: -------------------------------------------- 
    double vpj1[3] = {j1[1],j1[2],j1[3]};  // 3-momentum jet 1
    double vpj2[3] = {j2[1],j2[2],j2[3]};  // 3-momentum jet 2

    double pj1 = pow(pow(j1[1],2)+pow(j1[2],2)+pow(j1[3],2) , 0.5   );      // momentum jet 1
    double pj2 = pow(pow(j2[1],2)+pow(j2[2],2)+pow(j2[3],2) , 0.5   );      // momentum jet 2

    double ptj1 = pow(pow(j1[1],2)+pow(j1[2],2) , 0.5   );      // transverse momentum jet 1
    double ptj2 = pow(pow(j2[1],2)+pow(j2[2],2) , 0.5   );      // transverse momentum jet 2

    double MR2 = pow(pj1+pj2,2)-pow(j1[3]+j2[3],2);
    double MR = pow(MR2, 0.5);
    cout << "MR: "<< MR<< endl; 

    double vptm[2] = {-j1[1]-j2[1],-j1[2]-j2[2]};
    double ptm = pow(pow(vptm[0],2)+pow(vptm[1],2) ,0.5);

    double MRT2 = (ptm*(ptj1+ptj2)-(vptm[0]*(j1[1]+j2[1])+vptm[1]*(j1[2]+j2[2])))/2 ;
    double MRT = pow(MRT2,0.5);
    cout << "MRT: "<< MRT<< endl; 
    
    double R;
    
    // WHAT AM I DOING HERE????????
    try{ 
    if(MR == 0){
        cout << "entered if MR==0 if" <<endl;
        R = 1; 
    } 
    else{
        cout << " else: MR not 0 "<< endl; 
        R=MRT/MR;
    }
    }
    catch(const char* error)
    {
        cerr << "ERROR: " << error << " ";
    }
    //double R = MRT/MR;
    cout << "R: "<< R<< endl; 

    // fill the razor variable plots: 
    if(_jetDetectability == 0){
        _MR_T->Fill(MR);
        _MRT_T->Fill(MRT);
        _R_T->Fill(R);
        //_NJ_T->Fill(jets.size());
        _beta_T->Fill(beta);
    }
    if(_jetDetectability == 1){
        _MR_DAB->Fill(MR);
        _MRT_DAB->Fill(MRT);
        _R_DAB->Fill(R);
        //_NJ_DAB->Fill(jets.size());
        _beta_DAB->Fill(beta);
    }
    if(_jetDetectability == 2){
        _MR_DED->Fill(MR);
        _MRT_DED->Fill(MRT);
        _R_DED->Fill(R);
        //_NJ_DED->Fill(jets.size());
        _beta_DED->Fill(beta);
    }
    cout << "End EVENT "<< _nEvt<< endl;


    _nEvt ++ ; 
}



void twoPhotonJetRazor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void twoPhotonJetRazor::end(){ 
    _rootfile->Write();
}
