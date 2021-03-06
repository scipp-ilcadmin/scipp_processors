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

#include "tpJetRazor.h"
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


//static TH1F* _beta_T;
//static TH1F* _beta_DAB;
//static TH1F* _beta_DED;

static TH1F* _MR_TRU;
static TH1F* _MR_DAB;
static TH1F* _MR_DED;

static TH1F* _MRT_TRU;
static TH1F* _MRT_DAB;
static TH1F* _MRT_DED;

static TH1F* _R_TRU;
static TH1F* _R_DAB;
static TH1F* _R_DED;

static TH2F* _MRR_TRU;
static TH2F* _MRR_DAB;
static TH2F* _MRR_DED;

static TH2F* _MRR2_TRU;
static TH2F* _MRR2_DAB;
static TH2F* _MRR2_DED;

static TH1F* mult;
static TH2F* multjets;

tpJetRazor tpJetRazor;

tpJetRazor::tpJetRazor() : Processor("tpJetRazor") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "jetDetectability" ,
            "Detectability Level particles used in the Jet reconstruction:\n#\t0 : True\n#\t1 : Detectable\n#\t2 : Detected" ,
            _jetDetectability, 0);
    registerProcessorParameter("boost", 
            "Which R-frame transformation to do:  \n#\t0 : None \n#\t1 : Original (equalizes magnitude of 3-momenta) \n#\t2 : Modified (equalizes the z-momenta) \n#t3 : New (Using beta_{L}^{R}*, should always be physical)",
            _boost, 1);
}

void tpJetRazor::init() { 
    //streamlog_out(DEBUG)  << "   init called  " << endl;
    int nxbins = 400; 
    int nybins = 400; 
    Double_t xEdges[nxbins +1];
    Double_t yEdges[nybins +1];
    double p = -2;
    size_t size = 0;
    while(size < nxbins+1){
        xEdges[size++] = pow(10,p);
        p += 0.01;
    }
    size = 0; 
    p = -2; 
    while(size<nybins+1){
        yEdges[size++] = pow(10,p);
        p +=0.006;
    }



    if(_jetDetectability==0){_rootfile = new TFile("tpJetRazor_eW.pW.I39212._TRU.root","RECREATE");
        _R_TRU = new TH1F("R_TRU", "R=MTR/MR",1000,0,10);
        _MR_TRU = new TH1F("MR_TRU", "MR",500,0,100);
        _MRT_TRU = new TH1F("MRT_TRU", "MRT",250,0,50);
        _MRR_TRU = new TH2F("MRR_TRU", "MRR",500,0,100,1000,0,10 );
        _MRR2_TRU = new TH2F("MRR2_TRU", "MRR2",nxbins, xEdges,nybins, yEdges );
        mult = new TH1F("mult", "multiplicity", 500,0,500); 
        multjets = new TH2F("multjets", "jets v multiplicity", 500,0,500, 500,0,500); 
       // _beta_T = new TH1F("beta_T", "beta",250,0,50);
        
        freopen("tpJetRazor_eW.pW.I39212._TRU.log", "w",stdout);    
    }
    if(_jetDetectability==1){_rootfile = new TFile("tpJetRazor_eW.pW.I39212._DAB.root","RECREATE");
        _R_DAB = new TH1F("R_DAB", "R=MTR/MR",1000,0,10);
        _MR_DAB = new TH1F("MR_DAB", "MR",500,0,100);
        _MRT_DAB = new TH1F("MRT_DAB", "MRT",250,0,50);
        _MRR_DAB = new TH2F("MRR_DAB", "MRR",500,0,100, 1000,0,10);
        _MRR2_DAB = new TH2F("MRR2_DAB", "MRR2", nxbins, xEdges, nybins, yEdges);
        mult = new TH1F("mult", "multiplicity", 500,0,500); 
        multjets = new TH2F("multjets", "jets v multiplicity", 500,0,500, 500,0,500); 
      //  _beta_DAB = new TH1F("beta_DAB", "beta",130,-3,10);
        
        freopen("tpJetRazor_eW.pW.I39212._DAB.log", "w",stdout);   
    }
    if(_jetDetectability==2){_rootfile = new TFile("tpJetRazor_eW.pB.I39213._DED.root","RECREATE");
        _R_DED = new TH1F("R_DED", "R=MTR/MR",1000,0,10); 
        _MR_DED = new TH1F("MR_DED", "MR",500,0,100); 
        _MRT_DED = new TH1F("MRT_DED", "MRT",250,0,50); 
        _MRR_DED = new TH2F("MRR_DED", "MRR",500,0,100, 1000, 0, 10); 
        _MRR2_DED = new TH2F("MRR2_DED", "MRR2", nxbins, xEdges, nybins, yEdges); 
     //   _beta_DED = new TH1F("beta_DED", "beta",130,-3,10); 
        
        freopen("tpJetRazor_eW.pB.I39213._DED.log", "w",stdout);    
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
    pname[321]="                    K+ |";
    pname[-321]="                   K- |";
    pname[130]="                 K_L^0 |";
    pname[2112]="                    n |";
    pname[-2112]="                  -n |";
    pname[2212]="                    p |";
    pname[-2212]="                  -p |";
}



void tpJetRazor::processRunHeader( LCRunHeader* run) {  
    //    _nRun++ ; 

} 



void tpJetRazor::processEvent( LCEvent * evt ) { 
    
    //cerr << "EVENT: "<< _nEvt<< endl;   

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


    //cerr << "Included particles:" << endl; 
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
                        //cerr << "id: " << id << " "<< pname[id] << " "<< parp[0]<< " "<<parp[1]<< " "<< parp[2]<< endl;
                    }
                }
                if(_jetDetectability == 1){
                    if(isDetectable){
                        //cerr << "id: " << id << " "<< pname[id] << " "<< parp[0]<< " "<<parp[1]<< " "<< parp[2]<< endl;
                    }
                }
                if(_jetDetectability == 2){
                    if(isDetected){
                        //cerr << "id: " << id << " "<< pname[id] << " "<< parp[0]<< " "<<parp[1]<< " "<< parp[2]<< endl;
                    }
                }
            }// end particle not original 
        } //stat == 1 
    } // for particle

    Strategy strategy = Best;
    // identify the jets using the fastjet clustering algorithm:
    JetDefinition jet_def(ee_genkt_algorithm, _JetRParameter, -1 );
    //JetDefinition jet_def(antikt_algorithm, _JetRParameter );
    ClusterSequence cs(_parp, jet_def);
    // run the clustering, extract the jets

    //JetDefinition::Plugin * plugin = new EECambridgePlugin(0.5);
    //ClusterSequence cs(_parp, plugin);
    //vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets(0.5)); doesn't include all particles
    //vector<PseudoJet> jets = sorted_by_pt(cs.exclusive_jets()); doesn't include all particles 
    vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets()); //doesn't include all particles 
    //vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

    vector<PseudoJet> unclustered = cs.unclustered_particles();
    for (unsigned i = 0; i < unclustered.size(); i++) {
        //cerr << "jet " << i << ": "<< unclustered[i].pt() << " "<< unclustered[i].rap() << " " << unclustered[i].phi() << endl;
        vector<PseudoJet> constituents = unclustered[i].constituents();
        for (unsigned j = 0; j < unclustered.size(); j++) {
            //cerr << "constituent " << j << "’s pt: "<< constituents[j].pt() << endl;
        }
    }
    vector<PseudoJet> childless = cs.childless_pseudojets(); 
    
   
    for (unsigned i = 0; i < childless.size(); i++) {
        //cerr << "jet " << i << ": "<< childless[i].pt() << " "<< childless[i].rap() << " " << childless[i].phi() << endl;
        vector<PseudoJet> constituents = childless[i].constituents();
        for (unsigned j = 0; j < childless.size(); j++) {
            //cerr << "constituent " << j << "’s pt: "<< constituents[j].pt() << endl;
        }
    }
   
    //print the number of jets in event:
    //cerr << "NUMBER OF JETS: "<< jets.size()<< endl;
    totalJets += jets.size(); 
    //check if 0 jets:
    if(jets.size()==0){
        j0eventsCheck +=" ";
        j0eventsCheck +=std::to_string(_nEvt);
        j0eventsCheck +=" ";
        totalj0Events+=1; 

    }
    //check if 1 jet: 
    if(jets.size()==1){
        j1eventsCheck +=" ";
        j1eventsCheck +=std::to_string(_nEvt);
        j1eventsCheck +=" ";
        totalj1Events+=1; 
    }
    
    //multjets->Fill(_parp.size(), jets.size());
    //mult->Fill(_parp.size());
    vector<vector<PseudoJet>> megajet = getMegajets(jets); 
    
    
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

    cout <<"Megajet m2 value: " << endl;
    cout << "m2: "<< m2 << endl;

    vector<vector<double>> jR = boostMegajets(j1,j2);  

    // The Razor Variables *without* using the R-frame: -------------------------------------------- 
    double vpj1[3] = {jR[0][1],jR[0][2],jR[0][3]};  // 3-momentum jet 1
    double vpj2[3] = {jR[1][1],jR[1][2],jR[1][3]};  // 3-momentum jet 2

    double pj1 = pow(pow(jR[0][1],2)+pow(jR[0][2],2)+pow(jR[0][3],2) , 0.5   ); // momentum magnitude jet 1
    double pj2 = pow(pow(jR[1][1],2)+pow(jR[1][2],2)+pow(jR[1][3],2) , 0.5   );          // momentum jet 2

    double ptj1 = pow(pow(jR[0][1],2)+pow(jR[0][2],2) , 0.5   );      // transverse momentum jet 1
    double ptj2 = pow(pow(jR[1][1],2)+pow(jR[1][2],2) , 0.5   );      // transverse momentum jet 2

    double MR2 = pow(pj1+pj2,2)-pow(jR[0][3]+jR[1][3],2);
    double MR = pow(MR2, 0.5);
    cout << "MR: "<< MR<< endl; 
    //cerr << "MR: "<< MR<< endl; 

    double vptm[2] = {-jR[0][1]-jR[1][1],-jR[0][2]-jR[1][2]};
    double ptm = pow(pow(vptm[0],2)+pow(vptm[1],2) ,0.5);

    double MRT2 = (ptm*(ptj1+ptj2)-(vptm[0]*(jR[0][1]+jR[1][1])+vptm[1]*(jR[0][2]+jR[1][2])))/2.0 ;
    double MRT = pow(MRT2,0.5);
    cout << "MRT: "<< MRT<< endl; 
    //cerr << "MRT: "<< MRT<< endl; 
    
    double R;
     
    try{ 
        if(MR == 0){
            //cerr << "MR == 0" <<endl;
            MR0check++; 
            R = 0; 
        } 
        else{ 
            R=MRT/MR;
        }
    }
    catch(const char* error)
    {
        //cerr << "ERROR: " << error << " ";
    }
    //cerr << "R = "<< R<< endl;  
    cout << "R = "<< R<< endl;  
   
     
    double R2 = R*R;
    //underflow bins FOR THE SAKE OF LOG PLOT HISTOGRAM MRR2 :
    if (R2 < 0.01){
        R2 = 0.01;
    }
    if(MR < 0.01){
        MR = 0.01;
    }
    // reset overflow bins for log plots:  
    if(R2 > 1.4){
        R2 = 1.4;
    }
    if(MR > 100){
        MR = 100; 
    }
    //----------------------------------------
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
    if(R2<0.9 && MR>2){
        _cuts[4]+=1;
    }
    if(R2 > 0.015 && MR>1.5){
        _cuts[5]+=1;
    }
    if(R2 > 1.5*exp(-0.8*MR)){
        _concut += 1; 
    }
    if(R2 > 1.5*exp(-0.8*(MR- 0.5))){
        _concut_green += 1; 
    }
    if(R2 > 1.5*exp(-0.8*(MR- 1))){
        _concut_blue += 1; 
    }
    if(R2 > 1.5*exp(-0.8*(MR- 1.4))){
        _concut_yellow += 1; 
    }
    if(R>1.0){
        Rvals += " ";
        Rvals += std::to_string(R);
        Rvals += " ";
        totalRcheck++;
    }
    // fill the razor variable plots: 
    if(_jetDetectability == 0){
        _MR_TRU->Fill(MR);
        _MRT_TRU->Fill(MRT);
        _R_TRU->Fill(R);
        _MRR_TRU->Fill(MR,R);
        _MRR2_TRU->Fill(MR,R2);
        //_NJ_T->Fill(jets.size());
     //   _beta_T->Fill(beta);
    }
    if(_jetDetectability == 1){
        _MR_DAB->Fill(MR);
        _MRT_DAB->Fill(MRT);
        _R_DAB->Fill(R);
        _MRR_DAB->Fill(MR,R);
        _MRR2_DAB->Fill(MR,R2);
        //_NJ_DAB->Fill(jets.size());
     //   _beta_DAB->Fill(beta);
    }
    if(_jetDetectability == 2){
        _MR_DED->Fill(MR);
        _MRT_DED->Fill(MRT);
        _R_DED->Fill(R);
        _MRR_DED->Fill(MR,R);
        _MRR2_DED->Fill(MR,R2);
        //_NJ_DED->Fill(jets.size());
     //   _beta_DED->Fill(beta);
    }
    //cerr << "End EVENT "<< _nEvt<< endl;


    _nEvt ++ ; 
}



void tpJetRazor::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void tpJetRazor::end(){ 
    _rootfile->Write(); 
    cerr << "CUTS: "<< endl; 
    //cout << "CUTS: "<< endl; 
    for(int i = 0; i <6; i++){
        cerr << _cuts[i] << endl; 
    }
    double aveJets = totalJets/1600000; 
    //cerr << "Average # of Jets: "<< aveJets << endl ; 
    cerr << "CONTOUR CUT: "<< endl;
    cerr << _concut << endl;  
    cerr << "GREEN CONTOUR CUT: "<< endl;
    cerr << _concut_green << endl;  
    cerr << "BLUE CONTOUR CUT: "<< endl;
    cerr << _concut_blue << endl; 
    cerr << "------------------------------------"<< endl;  
    cerr << "YELLOW CONTOUR CUT: "<< endl;
    cerr << _concut_yellow  << endl;  
}

vector<vector<PseudoJet>> tpJetRazor::getMegajets(vector<PseudoJet> jets){
    // RAZOR MEGAJET ALGORITHM ------------------------------------------------
    vector<vector<PseudoJet>> megajet(2);
    int num_partitions = pow(2, jets.size());  
    double min_partition_m2 = numeric_limits<double>::max();
    int min_index = 0;

    for(int i = 0; i<num_partitions; ++i){
        vector<PseudoJet> partition[2];
        // the two partitions are symmetric under swapping so only make 1/2 of possible partions: 
        if((i & pow(2,0) != 0)){
            for(int j = 0; j<jets.size(); ++j){
                int bit = pow(2,j);

                if((i & bit) !=0){
                    partition[0].push_back(jets[j]);
                }
                else{
                    partition[1].push_back(jets[j]);
                }
            } // jets 
        } // if zero'th bit set  
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
      
        if(partition_m2<min_partition_m2 && i !=0 && i != (num_partitions-1) )  {
            min_partition_m2 = partition_m2;
            min_index = i;
        } 
    } // num_partitions

    for(int j=0; j<jets.size(); ++j){
        int bit = pow(2,j);
        if ((min_index & bit) != 0){
            megajet[0].push_back(jets[j]);
        }
        else {
            //subset does not contain jets[j]
            megajet[1].push_back(jets[j]);
        }
    }
    return megajet; 
}
vector<vector<double>> tpJetRazor::boostMegajets(vector<double> j1, vector<double> j2){

    vector<vector<double>> jR; 

    if(_boost == 1){
        double beta = (j1[0]-j2[0])/(j1[3]-j2[3]);
        cout << "beta: "<< beta<< endl;
        double gamma = pow((1-pow(beta,2)), -0.5);
        cout << "gamma: "<< gamma<< endl;

        //double MRj1 = 2*((pow(j1R[1],2), pow(j1R[2],2), pow(j1R[3],2)));
        //double MRj2 = 2*((pow(j2R[1],2), pow(j2R[2],2), pow(j2R[3],2)));
        //cout << "MRj1: "<< MRj1 << endl;
        //cout << "MRj2: "<< MRj2 << endl;

        if(beta>1){
            //cerr << "Event with beta > 1 !!!" << endl; 
            betaCheck += " ";
            betaCheck += std::to_string(_nEvt);
            betaCheck += " ";
            totalUnph +=1;  
            beta = .999;
            gamma = pow((1-pow(beta,2)), -0.5);
        }
        if(beta<-1){
            //cerr << "Event with beta > 1 !!!" << endl; 
            betaCheck += " ";
            betaCheck += std::to_string(_nEvt);
            betaCheck += " ";
            totalUnph +=1;            
            beta = -0.999;
            gamma = pow((1-pow(beta,2)), -0.5);
        }
        jR  = 
        {{gamma*j1[0]-gamma*beta*j1[3], j1[1],j1[2],-gamma*beta*j1[0]+gamma*j1[3] }, 
        {gamma*j2[0]-gamma*beta*j2[3], j2[1],j2[2],-gamma*beta*j2[0]+gamma*j2[3] } } ;
    }
    if(_boost == 3){
        double beta = (j1[3]+j2[3])/(j1[0]+j2[0]);  
         
        if(beta>1){
            //cerr << "Event with beta > 1 !!!" << endl; 
            betaCheck += " ";
            betaCheck += std::to_string(_nEvt);
            betaCheck += " "; 
        }
        double gamma = pow(1-pow(beta,2), -0.5);
        jR  = 
        {{gamma*j1[0]-gamma*beta*j1[3], j1[1],j1[2],-gamma*beta*j1[0]+gamma*j1[3] }, 
        {gamma*j2[0]-gamma*beta*j2[3], j2[1],j2[2],-gamma*beta*j2[0]+gamma*j2[3] } } ;
        
    }
        
        return jR;
}
