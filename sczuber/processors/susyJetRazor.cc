#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * author Summer Zuber 
 * January 18, 2017 
 *
 * Including fastjet/ClusterSequence.hh to use fastjet algorithms to identify jets. Uses those jets to make Razor Variables.  
 */

#include "susyJetRazor.h"
#include "scipp_ilc_utilities.h"

#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <TMath.h> 

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio> // trying to creat log file
//#include <map> 

#include <CLHEP/Vector/ThreeVector.h>
// #include <CLHEP/Random/RanluxEngine.h>

#include <EVENT/LCCollection.h>
#include <EVENT/LCIO.h>
#include <EVENT/MCParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCTOOLS.h>

using namespace lcio;
using namespace marlin;
using namespace fastjet;
using namespace std;

static TFile* _rootfile;


//static TH2F* _njbeta;

//static TH1F* _beta_T;
//static TH1F* _beta_DAB;
//static TH1F* _beta_DED;

static TH1F* _R_TRU;
static TH1F* _R_DAB;
static TH1F* _R_DED;

static TH1F* _MR_TRU;
static TH1F* _MR_DAB;
static TH1F* _MR_DED; 

static TH1F* _MRT_TRU;
static TH1F* _MRT_DAB;
static TH1F* _MRT_DED; 

static TH2F* _MRR_TRU;
static TH2F* _MRR_DAB;
static TH2F* _MRR_DED; 

static TH2F* _MRR2_TRU;
static TH2F* _MRR2_DAB;
static TH2F* _MRR2_DED; 

static TH2F* _MRTMR_TRU;
static TH2F* _MRTMR_DAB;
static TH2F* _MRTMR_DED; 

static TH1F* _NJ_TRU;   // histogram of the number of jets using TRUE jets
static TH1F* _NJ_DAB; // histogram of the number of jets using DAB jets
static TH1F* _NJ_DED; // histogram of the number of jets using DED jets

static TH2F* _NJR_TRU;   // histogram of the number of jets using TRUE jets
static TH2F* _NJR_DAB; // histogram of the number of jets using DAB jets
static TH2F* _NJR_DED; // histogram of the number of jets using DED jets
susyJetRazor susyJetRazor;

susyJetRazor::susyJetRazor() : Processor("susyJetRazor") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );
    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") ); 
    registerProcessorParameter( "jetDetectability",
            "Detectability of the Thrust Axis/Value to be used:\n#\t0 : True \n#t1 : Detectable \n#t2 : Detected" ,
            _jetDetectability, 2);
    registerProcessorParameter("boost", 
            "Which R-frame transformation to do:\n#\t0 : None \n#t1 : Original (equalizes magnitude of 3-momenta) \n#t2 : Modified (equalizes the z-momenta) \n#t3 : New (using beta_{L}^{R}*, should always be physical) ",
             _boost, 1 ); 
}


void susyJetRazor::init() { 
    streamlog_out(DEBUG)  << "   init called  " << std::endl ;
    
    int nxbins = 400; 
    int nybins = 400;  
    Double_t xEdges[nxbins+1];
    Double_t yEdges[nybins+1]; 
    double p = -2; 
    size_t size = 0;
    while( size < nxbins+1){
         xEdges[size++] = pow(10,p );
         p += 0.01; 
    }
    size = 0; 
    p = -2;
    while(size < nybins +1){
        yEdges[size++] = pow(10,p);
        p += 0.006; 
    }
    for(int i = 0 ; i<nxbins+1; i++){
        cout << xEdges[i] << endl; 
    }
    if(_jetDetectability==0){_rootfile = new TFile("susyJetRazor_.39117._TRU0.5.root","RECREATE");
        _R_TRU = new TH1F("R_TRU", "R =MTR/MR",1000,0,10); // the razor variable 
        _MR_TRU = new TH1F("MR_TRU","MR", 500, 0.0 ,100); // the M_{R} variable = 2|pR|
        _MRT_TRU = new TH1F("MRT_TRU","MRT", 500, 0 ,100); // the M_{T}^{R} variable 
        _MRR_TRU = new TH2F("MRR_TRU","MRR", 500, 0 ,100, 1000, 0,10); // the M_{T}^{R} variable 
        //_MRR2_TRU = new TH2F("MRR2_TRU","MRR2", 500, 0.1 ,100, 1000, 0.001, 1.4); // the M_{T}^{R} variable 
        _MRR2_TRU = new TH2F("MRR2_TRU","MRR2", nxbins, xEdges, nybins, yEdges); // the M_{T}^{R} variable  
        _MRTMR_TRU = new TH2F("MRTMR_TRU","MRTMR", 1000, 0 ,1000, 1000, 0,1000); // the M_{T}^{R} variable 
        _NJ_TRU = new TH1F("NJ_TRU","NJ",20, 0, 20); 
        _NJR_TRU = new TH2F("NJR_TRU","NJR",20, 0, 20, 1000,0,2); 
       // _beta_T = new TH1F("beta_T","beta",80,-20,20);
       // _njbeta = new TH2F("njbeta","njbeta",40,-10,20,40,-20,20);
        
        freopen( "susyJetRazor_.39133._DED0.5.log", "w", stdout ); 
    }
    if(_jetDetectability==1){_rootfile = new TFile("susyJetRazor_.39133._DAB0.5.root","RECREATE");
        _R_DAB = new TH1F("R_DAB", "R =MTR/MR",1000,0,10);
        _MR_DAB = new TH1F("MR_DAB","MR", 500, 0 ,100); 
        _MRT_DAB = new TH1F("MRT_DAB","MRT", 250, 0 ,50); 
        _MRR_DAB = new TH2F("MRR_DAB","MRR", 500, 0 ,100, 1000,0,10); 
        _MRR2_DAB = new TH2F("MRR2_DAB","MRR2", nxbins, xEdges, nybins,yEdges); 
        _MRTMR_DAB = new TH2F("MRTMR_DAB","MRTMR", 1000, 0 ,1000, 1000,0,1000); 
        _NJ_DAB = new TH1F("NJ_DAB","NJ",20,0,20);
        _NJR_DAB = new TH2F("NJR_DAB","NJR",20, 0, 20, 1000,0,2); 
       // _beta_DAB = new TH1F("beta_DAB","beta",80,-20,20);
        //_njbeta = new TH2F("njbeta","njbeta",20,0,20,200,-20,20);
        
        freopen( "susyJetRazor_.39133._DAB0.5.log", "w", stdout ); 
    }
    if(_jetDetectability==2){_rootfile = new TFile("susyJetRazor_.39133._DED0.5.root","RECREATE");
        _R_DED = new TH1F("R_DED", "R =MTR/MR",1000,0,10);
        _MR_DED = new TH1F("MR_DED","MR", 500, 0 ,100); 
        _MRT_DED = new TH1F("MRT_DED","MRT", 500, 0 ,100); 
        _MRR_DED = new TH2F("MRR_DED","MRR", 500, 0 ,100, 1000,0,10); 
        _MRR2_DED = new TH2F("MRR2_DED","MRR2", nxbins, xEdges, nybins,yEdges); 
        _MRTMR_DED = new TH2F("MRTMR_DED","MRTMR", 1000, 0, 1000, 1000,0,1000); 
        _NJR_DED = new TH2F("NJR_DED","NJR",20, 0, 20, 1000,0,2); 
       // _beta_DED = new TH1F("beta_DED","beta",40,-10,20);
        
        freopen( "susyJetRazor_.39133._DED0.5.log", "w", stdout ); 
    }
    // irameters() ;

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

}

void susyJetRazor::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ;
} 

void susyJetRazor::processEvent( LCEvent * evt ) { 

    cout << "EVENT: " << _nEvt << endl; 
    _inParVec = evt->getCollection( _colName) ;
    //cout << "num of elements " << _inParVec->getNumberOfElements() << endl;
    if (!_parp.empty()) _parp.clear();
    int id, stat;

    cout << "Final State Particles: "<< endl;
    // loop through particles and add to _parp: 
    for (int n=0;n<_inParVec->getNumberOfElements() ;n++){
        MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );
        try{
            id = aPart->getPDG();
            stat = aPart->getGeneratorStatus();
        }
        catch(const std::exception& e){
            cout << "exception caught with message " << e.what() << "\n";
        }
        const double* parp = aPart->getMomentum();
        double parp_mag = sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]);
        
        if(stat==1){
            cout << "id: " << id << " " << pname[id] <<" " << parp[0]<<" "<< parp[1]<<" "<< parp[2] <<endl; 
            double penergy = aPart->getEnergy(); 
            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);

            double cos = parp[2]/parp_mag;
            bool isForward = ( cos > 0.9 || cos < - 0.9);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable &&  !isForward  );
            
           
            if(_jetDetectability == 0){
                if(!isDarkMatter){ 
                    _parp.push_back( PseudoJet(parp[0], parp[1], parp[2], penergy) );
                   
                }
            }
            if(_jetDetectability == 1){
                if(isDetectable){ 
                    _parp.push_back( PseudoJet(parp[0], parp[1], parp[2], penergy) );
                  
                }
            }

            if(_jetDetectability == 2){ 
                if(isDetected){  
                    _parp.push_back( PseudoJet(parp[0], parp[1], parp[2], penergy) ); 
                 
                }
            }
        } // stat = 1
    } // for particle
    
    // print the included particles:
    cout <<"----------------------------------------------"<< endl;  
    cout<< "Included Particles: "<< endl; 
    //cerr<< "Included Particles: "<< endl; 
    for (int n=0;n<_inParVec->getNumberOfElements() ;n++){
        MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );
        try{
            id = aPart->getPDG();
            stat = aPart->getGeneratorStatus();
        }
        catch(const std::exception& e){
            cout << "exception caught with message " << e.what() << "\n";
        }
        const double* parp = aPart->getMomentum();
        double parp_mag = sqrt(parp[0]*parp[0]+parp[1]*parp[1]+parp[2]*parp[2]);
        
        if(stat==1){
            double penergy = aPart->getEnergy(); 
            bool isDarkMatter = (id == 1000022);
            bool isNeutrino = (
                    id == 12 || id == -12 ||
                    id == 14 || id == -14 ||
                    id == 16 || id == -16 ||
                    id == 18 || id == -18);

            double cos = parp[2]/parp_mag;
            bool isForward = ( cos > 0.9 || cos < - 0.9);
            bool isDetectable = (!isDarkMatter && !isNeutrino);
            bool isDetected = (isDetectable &&  !isForward  );
            
           
            if(_jetDetectability == 0){
                if(!isDarkMatter){  
                    //cerr << "id: " << id << " " << pname[id] <<" " << parp[0]<<" "<< parp[1]<<" "<< parp[2] <<endl;   
                }
            }
            if(_jetDetectability == 1){
                if(isDetectable){    
                    cout << "id: " << id << " " << pname[id] <<" " << parp[0]<<" "<< parp[1]<<" "<< parp[2] <<endl; 
                }
            }
            if(_jetDetectability == 2){ 
                if(isDetected){    
                    //cerr << "id: " << id << " " << pname[id] <<" " << parp[0]<<" "<< parp[1]<<" "<< parp[2] <<endl; 
                }
            }
        } // stat = 1
    } // end particle loop  
    
    // identify the jets using the fastjet cluserting algorithm:
    //JetDefinition jet_def(antikt_algorithm,   _JetRParameter   ); 
    JetDefinition jet_def(ee_genkt_algorithm, _JetRParameter, -1 ); 
    
    // run the clustering, extract the jets
    ClusterSequence cs(_parp, jet_def);
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets()); 
    cerr << "NUMBER OF JETS: "<< jets.size() << endl;
    
    int NJ = jets.size(); 
    //  check if 0 jets in event:
    if(jets.size()==0){
        //cerr << "Event with 0 jets"<< endl; 
        j0eventsCheck += " ";
        j0eventsCheck += std::to_string(_nEvt);
        j0eventsCheck += " "; 
        numj0eventsCheck+=1; 
    }
    // check if 1 jet in event:
    if(jets.size()==1){
        //cerr << "Event with 1 jet" << endl; 
        j1eventsCheck += " ";
        j1eventsCheck += std::to_string(_nEvt);
        j1eventsCheck += " ";  
        numj1eventsCheck +=1; 
    }

    //vector<PseudoJet> megajet = getMegajets(jets);
    vector<vector<PseudoJet>> megajet = getMegajets(jets);
    //cerr << megajet.size() << endl; 

    // ----------------------------------------------------------------------------------------------
    
    // CALCULATE THE MR, MRT, AND R VARIABLED:
    
    // the 4-momenta of the two megajets: 
    vector<double> j1(4);
    vector<double> j2(4);
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
    
    //cerr <<"megajet m2 value:   " << m2 <<endl;  
    vector<vector<double>> jR = boostMegajets(j1, j2);  

    // The Razor Variables using whatever R-frame is selcted: -------------------------------------------- 
    cerr << "4-momenta jet 1 : "<< jR[0][0] <<' '<<jR[0][1]<< ' '<<jR[0][2]<< ' '<<jR[0][3]<< endl;  
    cerr << "4-momenta jet 2 : "<< jR[1][0] <<' '<<jR[1][1]<< ' '<<jR[1][2]<< ' '<<jR[1][3]<< endl;  
    
    double vpj1[3] = {jR[0][1],jR[0][2],jR[0][3]};  // 3-momentum jet 1
    double vpj2[3] = {jR[1][1],jR[1][2],jR[1][3]};  // 3-momentum jet 2

    double pj1 = pow(pow(jR[0][1],2)+pow(jR[0][2],2)+pow(jR[0][3],2) , 0.5   );      // momentum jet 1
    cerr << "momentum jet 1 = "<< pj1 << endl; 
    cout << "momentum jet 1 = "<< pj1 << endl; 
    double pj2 = pow(pow(jR[1][1],2)+pow(jR[1][2],2)+pow(jR[1][3],2) , 0.5   );      // momentum jet 2
    cerr << "momentum jet 2 = "<< pj2 << endl; 
    cout << "momentum jet 2 = "<< pj2 << endl; 
  
    cerr << "pzj1: "<< jR[0][3] << " pzj2: "<< jR[1][3] << endl;  
    cout << "pzj1: "<< jR[0][3] << " pzj2: "<< jR[1][3] << endl;  
    
    double ptj1 = pow(pow(jR[0][1],2)+pow(jR[0][2],2) , 0.5   );      // transverse momentum jet 1
    double ptj2 = pow(pow(jR[1][1],2)+pow(jR[1][2],2) , 0.5   );      // transverse momentum jet 2
   
    cout << "ptj1: "<< ptj1<< " ptj2: "<< ptj2 << endl;  
    //cerr << "ptj1: "<< ptj1<< " ptj2: "<< ptj2 << endl;  
    double MR2 = pow(pj1+pj2 ,2)-pow(jR[0][3]+jR[1][3], 2); 
    double MR = pow(MR2, 0.5);
    
    // vector missing transverse energy  
    double vptm[2] = {-jR[0][1]-jR[1][1],-jR[0][2]-jR[1][2]};
    double ptm = pow(pow(vptm[0],2)+pow(vptm[1],2) ,0.5); 
    
    double MRT2 = (ptm*(ptj1+ptj2)-(vptm[0]*(jR[0][1]+jR[1][1])+vptm[1]*(jR[0][2]+jR[1][2])))/2 ;
    double MRT = pow(MRT2, 0.5); 
  
    cerr << "MR = sqrt( (pj1+pj2)**2 - (pzj1+pzj2)**2 ) = "<< MR<< endl;  
   
    cout << "MR="<< MR<< endl;  
    cout << "MRT="<< MRT<< endl;  
    cerr << "MR="<< MR<< endl;  
    cerr << "MRT="<< MRT<< endl;  
    
    
    double R;


    //_cuts = doCuts(MR,R*R); 
    
    try{
        if (MR == 0){
            cerr << "found event with MR == 0. Setting R = 1. " << endl;
            mr0check++; 
            R = 0;
        }
        else{
            R = MRT/MR;
        }
    }
    catch (const char* error){
        cerr << "ERROR: " << error << " ";    
    }
    double R2 = R*R;
  
    cout << "R="<<R<< endl; 
    cerr << "R="<<R<< endl; 
    if (MR < 0.001){
        R = 0; 
    }
  
   // reset the underflow bin for log plots:
    if(R2 < 0.01){ // they both start at 10^-2
        R2 = 0.01;
    }
    if(MR < 0.01){
        MR = 0.01;
    }
    // ------------------------------------------------
    if(MR > 2){
        _cuts[0]+=1;
        if(R2 > 0.05){
            _cuts[1]+=1;
        }
    }
    if(R2>0.1 && MR>3){
        _cuts[2]+=1;
    }
    if(R2<0.9 && MR>1){
        _cuts[3]+=1;
    }
    if(R2 > 0.05 && MR > 2){
        _cuts[4]+=1;
    }
    if(R > 1){
        
        Rvals += " ";
        Rvals += std::to_string(R);
        Rvals += " ";
    }


    if (R>1.2){
        cerr << "FOUND EVENT WITH R>1.2!!!"<< endl; 
        Rcheck += " ";
        Rcheck += std::to_string(_nEvt);
        Rcheck += " ";
    }
    
    if(R2 > 1.4){
        R2 = 1.4;
    }
    if(MR > 100){
        MR = 100; 
    }
    // fill the razor variable plots: 
    if(_jetDetectability == 0){ 
        _MR_TRU->Fill(MR);
        _MRT_TRU->Fill(MRT);
        _R_TRU->Fill(R);
        _MRR_TRU->Fill(MR,R);
        _MRR2_TRU->Fill(MR,R2);
        _MRTMR_TRU->Fill(MRT,MR);
        _NJR_TRU->Fill(NJ,R);
        //_NJ_T->Fill(jets.size());
        //_beta_T->Fill(beta);
     }
     if(_jetDetectability == 1){
        _MR_DAB->Fill(MR);
        _MRT_DAB->Fill(MRT);
        _R_DAB->Fill(R);
        _MRR_DAB->Fill(MR,R);
        _MRR2_DAB->Fill(MR,R2);
        _MRTMR_DAB->Fill(MRT,MR);
        _NJR_DAB->Fill(NJ,R);
        //_NJ_DAB->Fill(jets.size());
        //_beta_DAB->Fill(beta);
     }
     if(_jetDetectability == 2){
        _MR_DED->Fill(MR);
        _MRT_DED->Fill(MRT);
        _R_DED->Fill(R);
        _MRR_DED->Fill(MR,R);
        _MRR2_DED->Fill(MR,R2);
        _MRTMR_DED->Fill(MRT,MR);
        _NJR_DED->Fill(NJ,R);
        //_NJ_DED->Fill(jets.size());
        //_beta_DED->Fill(beta);
     }
    cout << "End EVENT "<< _nEvt<< endl;
    _nEvt ++ ; 
}

void susyJetRazor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}

void susyJetRazor::end(){ 
    _rootfile->Write();
    //cerr << "Events with 0 jets: "<< j0eventsCheck << endl; 
    //cerr << "Events with 1 jets: "<< j1eventsCheck << endl;
    cerr << "Total # events with 0 jet: "<< numj0eventsCheck << endl; 
    cerr << "Total # events with 1 jet: "<< numj1eventsCheck << endl; 
    cout << "Total # events with 0 jet: "<< numj0eventsCheck << endl; 
    cout << "Total # events with 1 jet: "<< numj1eventsCheck << endl; 
    //cout << parpCheck << endl;
    //cerr << "Events with beta > 1: "<< betaCheck << endl;
    cerr << "Events with R>1.2: " << Rcheck << endl ; //use cerr to pring in terminal rather than in log file   
    cout << "Events with R>1.2: " << Rcheck << endl ; //use cerr to pring in terminal rather than in log file   
    streamlog_out(DEBUG)  << "   Events with R > 1.20:  " <<Rcheck<< std::endl ;
    cout << "CUTS: "<< endl; 
    for(int i = 0; i< 5; i++){
        cout << _cuts[i] << endl; 
    }
    cerr << "CUTS: "<< endl; 
    for(int i = 0; i< 5; i++){
        cerr << _cuts[i] << endl; 
    }
    cerr << "# of events with MR == 0 : "<< mr0check<< endl;
    cerr << "R VALUES OVER 1 : "<< Rvals<< endl; 
    cerr << "total # Unphysical R-frame: "<< totalUnph;
    cerr <<" "<<  endl;  
}

//vector<PseudoJet> susyJetRazor::getMegajets(vector<PseudoJet> jets){
vector<vector<PseudoJet>> susyJetRazor::getMegajets(vector<PseudoJet> jets){
    vector<vector<PseudoJet>> megajet(2); 
    // RAZOR MEGAJET ALGORITHM -----------------------------------------------
    // megajet algorithm: the combination of jets into two disjoint sets (megajets) is the combination that minimized the 
    // sum of the squared-invariant-mass values of the two megajets 
    // the list of m2 values so that the minumum of these may be found 
    // the number of subsets of the set of jets: 
    // (this includes the null set and the set of jets itself, also includes double counting)
    int num_partitions = pow(2,jets.size()); 
    // array of squared-invariant-mass values (for each possible partitioning of jets): 
    double m2_values[num_partitions];
    // clear all m2 values before calculating all m2 values: 
    int p = 0;
    while(p<num_partitions){
        m2_values[p]=0;
        p++;
    }
    // calculate the m2 values for each i 
    for(unsigned int i = 0; i<num_partitions; i++){
        //cout << "*********************************************" << endl;
        //cout << "partition: " << i << endl; 
        vector<PseudoJet> partition[2]; 
        // create partition i (one subset)
        for(unsigned int j=0; j<jets.size(); j++){
            unsigned int bit = pow(2,j);
            if ((i & bit) != 0){ 
                partition[0].push_back(jets[j]);
            }
            else{
                partition[1].push_back(jets[j]);
            }
        }

        ////cout << "size: "<<partition[1].size() << "   " << partition[0].size() << endl;
        double subset_p[2][4]= {{0,0,0,0},{0,0,0,0}}; // subset 1 and 2: e, px, py, pz  
        // add the four momenta of jets in partition vectorally 
        for(int j=0; j<partition[0].size(); j++){
                subset_p[0][0]=subset_p[0][0]+partition[0][j].e();
                subset_p[0][1]=subset_p[0][1]+partition[0][j].px();
                subset_p[0][2]=subset_p[0][2]+partition[0][j].py();
                subset_p[0][3]=subset_p[0][3]+partition[0][j].pz();
        }
        for(int j=0; j<partition[1].size(); j++){
                subset_p[1][0]=subset_p[1][0]+partition[1][j].e();
                subset_p[1][1]=subset_p[1][1]+partition[1][j].px();
                subset_p[1][2]=subset_p[1][2]+partition[1][j].py();
                subset_p[1][3]=subset_p[1][3]+partition[1][j].pz();
        }
        // add squared-invariant-mass values together
        double subset_m2[2]= {0,0};
        subset_m2[0]=pow(subset_p[0][0],2)-(pow(subset_p[0][1],2)+pow(subset_p[0][2],2)+pow(subset_p[0][3],2));
        subset_m2[1]=pow(subset_p[1][0],2)-(pow(subset_p[1][1],2)+pow(subset_p[1][2],2)+pow(subset_p[1][3],2));
        double partition_m2 = subset_m2[0]+subset_m2[1]; // sum of the two megajets 
        // print the m2 value for this partition: 
        m2_values[i]= partition_m2;
        //cout <<"m2: "<< m2_values[i] << endl ;

   
    }
    // set the min m2 value to the last valid 
    double min_value = m2_values[num_partitions-2]; // MIGHT BE PROBLEM IF *NO* JETS OR *ONE* JET!!
    int min_index = num_partitions-2; 
    //cout <<"initial min_value: "<< min_value<< endl; 
    for(int i=0; i<num_partitions; i++){
        // if subset is not null set or set itself:
        if(i!=0 && i!=(num_partitions-1) ){
            if(m2_values[i]<min_value){
                min_index = i;
                min_value = m2_values[i];   
            }
        }
    }
    //cout <<"~~minimum index: "<< min_index<<" ~~minimum value: "<< min_value<< endl;
    // re-create the partition corresponding to the min_index: 
    //vector<PseudoJet> megajet[2]; 
    for(unsigned int j=0; j<jets.size(); j++){
        unsigned int bit = pow(2,j);
        if ((min_index & bit) != 0){ 
            megajet[0].push_back(jets[j]);
        } 
        else {
            //subset does not contain jets[j]
            megajet[1].push_back(jets[j]);
        }
    } 
    
    // megajet is now the correct megajet to use
    return megajet;
}

vector<vector<double>> susyJetRazor::boostMegajets(vector<double> j1, vector<double> j2){
    
    vector<vector<double>> jR; 

    if(_boost == 0){
            // this is the no boost case
           jR = {j1,j2}; 
    }
    
    if(_boost == 1){    
        // this is the boost to the original R-frame with beta = (e1-e2)/(pz1-pz2)   
        double beta = (j1[0]-j2[0])/(j1[3]-j2[3]); 
        cout << "beta: "<< beta<< endl; 
        //cerr << "beta: "<< beta<< endl; 
        double gamma = pow((1-pow(beta,2)), -0.5);
        cout << "gamma: "<< gamma<< endl;
        if(beta>1){
            //cerr << "Event with beta > 1 !!!" << endl;     
            betaCheck += " ";
            betaCheck += std::to_string(_nEvt);
            betaCheck += " ";  
            totalUnph +=1; 
            beta = 0.999;
            //cerr << "new beta: "<< beta<< endl; 
            gamma = pow((1-pow(beta,2)), -0.5);
            //cerr << "new gamma: "<< gamma<< endl; 
        }
        if(beta< -1){
            //cerr << "Event with beta > 1 !!!" << endl;     
            betaCheck += " ";
            betaCheck += std::to_string(_nEvt);
            betaCheck += " ";  
            totalUnph +=1;
            beta = -0.999;
            //cerr << "new beta: "<< beta<< endl; 
            gamma = pow((1-pow(beta,2)), -0.5);
            //cerr << "new gamma: "<< gamma<< endl; 
        }
        jR = { 
             {gamma*j1[0]-gamma*beta*j1[3], j1[1],j1[2],-gamma*beta*j1[0]+gamma*j1[3] },
             {gamma*j2[0]-gamma*beta*j2[3], j2[1],j2[2],-gamma*beta*j2[0]+gamma*j2[3] }} ;
  
    }
    
    //double MRj1 = 2*((pow(j1R[1],2), pow(j1R[2],2), pow(j1R[3],2)));
    //double MRj2 = 2*((pow(j2R[1],2), pow(j2R[2],2), pow(j2R[3],2)));
    //cerr << "MRj1: "<< MRj1 << endl; 
    //cerr << "MRj2: "<< MRj2 << endl; 
 
    //_njbeta->Fill(jets.size(),beta);
    if(_boost == 2){
        // this is the new one that is supposed to equalize the z-momenta of the two jets 
        double beta = (-j1[3]+j2[3])/(-j1[0]+j2[0]);  
        
        if(beta>1){
            cerr << "Found event with beta>1 "<<beta<< endl;  
            betaCheck += " ";
            betaCheck += std::to_string(_nEvt);
            betaCheck += " ";  
        }
        double gamma = pow((1-pow(beta,2)), -0.5);
        jR = { 
              {gamma*j1[0]-gamma*beta*j1[3], j1[1],j1[2],-gamma*beta*j1[0]+gamma*j1[3] },
              {gamma*j2[0]-gamma*beta*j2[3], j2[1],j2[2],-gamma*beta*j2[0]+gamma*j2[3] }} ; 
    }
    if(_boost == 3){
        // the one that is supposed to be always physical 
        double beta = (j1[3]+j2[3])/(j1[0]+j2[0]);
        double gamma = pow(1-pow(beta,2) , -0.5);
        cerr << "beta="<<beta<< endl; 
        if(beta>1){
            cerr << "Event with beta>1 !!!"<<beta<< endl;  
            betaCheck += " ";
            betaCheck += std::to_string(_nEvt);
            betaCheck += " ";  
        }
        jR = { 
              {gamma*j1[0]-gamma*beta*j1[3], j1[1],j1[2],-gamma*beta*j1[0]+gamma*j1[3] },
              {gamma*j2[0]-gamma*beta*j2[3], j2[1],j2[2],-gamma*beta*j2[0]+gamma*j2[3] }} ; 
    }
    
    return jR;
    
}
void susyJetRazor::BinLogX(TH2*h) {
    TAxis *axis = h->GetXaxis(); 
    int bins = axis->GetNbins();

    Axis_t from = axis->GetXmin();
    Axis_t to = axis->GetXmax();
    Axis_t width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins + 1];

    for (int i = 0; i <= bins; i++) {
        new_bins[i] = TMath::Power(10, from + i * width);
    } 
    axis->Set(bins, new_bins); 
    delete new_bins;
}
/*void susyJetRazor::BinLogY(TH2*h) {
    TAxis *axis = h->GetYaxis(); 
    int bins = axis->GetNbins();

    Axis_t from = axis->GetYmin();
    Axis_t to = axis->GetYmax();
    Axis_t width = (to - from) / bins;
    Axis_t *new_bins = new Axis_t[bins + 1];

    for (int i = 0; i <= bins; i++) {
        new_bins[i] = TMath::Power(10, from + i * width);
    } 
    axis->Set(bins, new_bins); 
    delete new_bins; 
}*/
