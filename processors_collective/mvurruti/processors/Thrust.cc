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
 */

#include "Thrust.h"
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


static TFile* _rootfile;

Thrust Thrust;

Thrust::Thrust() : Processor("Thrust") {
    // modify processor description
    _description = "Protype Processor" ;

    // register steering parameters: name, description, class-variable, default value
    registerInputCollection( LCIO::MCPARTICLE, "CollectionName" , "Name of the MCParticle collection"  , _colName , std::string("MCParticle") );

    registerProcessorParameter( "RootOutputName" , "output file"  , _root_file_name , std::string("output.root") );
    registerProcessorParameter( "typeOfThrustFinder" ,
     "Type of thrust reconstruction algorithm to be used:\n#\t1 : Tasso algorithm\n#\t2 : JetSet algorithm"  ,
      _typeOfThrustFinder , 2 ) ;
}



void Thrust::init() { 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;

    _rootfile = new TFile("thrust_test.root","RECREATE");
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

}



void Thrust::processRunHeader( LCRunHeader* run) { 
    //run->parameters().setValue("thrust",12300321);
    //    _nRun++ ;

} 



void Thrust::processEvent( LCEvent * evt ) { 
    // this gets called for every event 
    // usually the working horse ...

 

   _inParVec = evt->getCollection( _colName) ;

    if (!_partMom.empty()) _partMom.clear();

    for (int n=0;n<_inParVec->getNumberOfElements() ;n++)
    {
      MCParticle* aPart = dynamic_cast<MCParticle*>( _inParVec->getElementAt(n) );

      const double* partMom = aPart->getMomentum();
      _partMom.push_back( Hep3Vector(partMom[0], partMom[1], partMom[2]) ); 

    _nEvt ++ ;
    }

  //reset variables for output   
  _principleThrustValue = -1;
  _majorThrustValue     = -1;
  _minorThrustValue     = -1;
  _principleThrustAxis.set(0,0,0);
  _majorThrustAxis.set(0,0,0);
  _minorThrustAxis.set(0,0,0);

  // Switch to the desired type of thrust finder
  if (_typeOfThrustFinder == 1)
    { 
      TassoThrust();
    }
  else if (_partMom.size()<=1)
    { 
      TassoThrust();
    }
  else if (_typeOfThrustFinder == 2)
    {
      JetsetThrust();
    }
  // ###write
  //    evt->parameters().setValue("thrust",_principleThrustValue);

  FloatVec thrax;
  thrax.clear();
  thrax.push_back(_principleThrustAxis.x());
  thrax.push_back(_principleThrustAxis.y());
  thrax.push_back(_principleThrustAxis.z());

  _inParVec->parameters().setValue("principleThrustValue",_principleThrustValue);
  _inParVec->parameters().setValues("principleThrustAxis",thrax);

  if (_typeOfThrustFinder == 2)
    {
      thrax.clear();
      thrax.push_back(_majorThrustAxis.x());
      thrax.push_back(_majorThrustAxis.y());
      thrax.push_back(_majorThrustAxis.z());

      _inParVec->parameters().setValue("majorThrustValue",_majorThrustValue);
      _inParVec->parameters().setValues("majorThrustAxis",thrax);

      thrax.clear();
      thrax.push_back(_minorThrustAxis.x());
      thrax.push_back(_minorThrustAxis.y());
      thrax.push_back(_minorThrustAxis.z());

      _inParVec->parameters().setValue("minorThrustValue",_minorThrustValue);
      _inParVec->parameters().setValues("minorThrustAxis",thrax);

      float Oblateness;
      Oblateness = _majorThrustValue - _minorThrustValue;
      _inParVec->parameters().setValue("Oblateness",Oblateness);
      if ( (_majorThrustValue < 0) || (_minorThrustValue < 0) )
    {
      _inParVec->parameters().setValue("Oblateness",-1);
    }
    }

    streamlog_out( DEBUG4 ) << " thrust: " << _principleThrustValue << " TV: " << _principleThrustAxis << endl;
    streamlog_out( DEBUG4 ) << "  major: " << _majorThrustValue << " TV: " << _majorThrustAxis << endl;
    streamlog_out( DEBUG4 ) << "  minor: " << _minorThrustValue << " TV: " << _minorThrustAxis << endl;
    cout << "EVENT: " << _nEvt << endl;
    cout << " thrust: " << _principleThrustValue << " TV: " << _principleThrustAxis << endl;
    cout << "  major: " << _majorThrustValue << " TV: " << _majorThrustAxis << endl;
    cout << "  minor: " << _minorThrustValue << " TV: " << _minorThrustAxis << endl;

  if (_principleThrustValue >= _max) _max = _principleThrustValue;
  if (_principleThrustValue <= _min) _min = _principleThrustValue;
}




void Thrust::check( LCEvent * evt ) { 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void Thrust::end(){ 
    _rootfile->Write();
}

int Thrust::JetsetThrust(){
  const int nwork=11,iFastMax = 4,iGood=2;
  const float dConv=0.0001; // 0.0001
  int sgn;
  double theta=0,phi=0;
  double thp,thps,tds,tmax,dOblateness;
  vector<Hep3Vector> TAxes(3),Fast(iFastMax+1),Workv(nwork);
  vector<double> Workf(nwork),dThrust(3);
  Hep3Vector tdi,tpr,mytest;

  tmax = 0;
  for ( unsigned int i=0; i < _partMom.size(); i++)
    tmax += _partMom[i].mag();

  // pass = 0: find thrust axis
  // pass = 1: find major axis
  for ( int pass=0; pass <= 1; pass++ )
    {
      if ( pass == 1 )
    {
      phi   = TAxes[0].phi();
      theta = TAxes[0].theta();
      for ( unsigned  int i = 0;i < _partMom.size(); i++)
        {
          _partMom[i].rotateZ(-phi);
          _partMom[i].rotateY(-theta);
        }
      TAxes[0].set(0,0,1);
    } // if pass == 1

      // Find the ifast highest momentum particles and
      // put the highest in Fast[0], next in Fast[1],....Fast[iFast-1].
      // Fast[iFast] is just a workspace.
    for ( unsigned  int i = 0; i < Fast.size(); i++ )
    Fast[i].set(0,0,0);

      for ( unsigned int i = 0; i < _partMom.size(); i++ )
    {
      for ( int ifast = iFastMax -1; ifast >= 0 ; ifast-- )
        {
          if (_partMom[i].mag2() > Fast[ifast].mag2() )
        {
          Fast[ifast + 1] = Fast[ifast];
          if (ifast == 0) Fast[ifast] = _partMom[i];
        }
          else
        {
          Fast[ifast + 1] = _partMom[i];
          break;
        } // if p>p_fast
        } // for ifast 
    } // for i 

      // Find axis with highest thrust (case 0)/ highest major (case 1).
for ( unsigned int iw = 0; iw < Workv.size(); iw++ )
    {
      Workf[iw] = 0.;
    }
      int p = (int) min( iFastMax,_partMom.size() ) - 1 ;
      int nc = 1 << p;
      for ( int n = 0; n < nc; n++ )
    {
      tdi.set(0,0,0);
      for (int i = 0; i < min(iFastMax,nc) ; i++)
        {
          if ( (1 << (i+1)) * ( (n + (1<<i)) / (1<<(i+1)) ) >= n+1) //i+1 
        { sgn = -1;} else {sgn = 1;}
          tdi += sgn*Fast[i];
          if (pass==1) tdi.setZ(0);
        } // for i 
      tds = tdi.mag2();
      for ( int iw = (int) min(n,9); iw >= 0; iw-- )
        {
          if (tds > Workf[iw])
        {
          Workf[iw+1] = Workf[iw];
          Workv[iw+1] = Workv[iw];
          if (iw == 0)
            { Workv[iw] = tdi; Workf[iw] = tds;}
        }
          else // if tds 
        {
          Workv[iw+1] = tdi;
          Workf[iw+1] = tds;
        } // if tds 
        } // for iw
    } // for n 

      // Iterate direction of axis until stable maximum.
    dThrust[pass] = 0;
      int nagree = 0;
      for ( int iw = 0; iw < min(nc,10) && nagree < iGood; iw++ )
    {
      thp = 0;
      thps = -99999.;
      while ( thp > thps + dConv )
        {
          thps = thp;
          if ( thp <= 1E-10 )
        { tdi = Workv[iw]; } else { tdi=tpr; }
          tpr.set(0,0,0);
          for ( unsigned int i = 0; i < _partMom.size(); i++ )
        {
          sgn = (int) sign(1,tdi.dot(_partMom[i]));
          tpr += sgn*_partMom[i];
          if (pass == 1) { tpr.setZ(0); } // ###
        } // for i 
          thp = tpr.mag()/tmax;
        } // while 
      // Save good axis. Try new initial axis until enough
      // tries agree.
      if ( thp < dThrust[pass] - dConv ) continue;
      if ( thp > dThrust[pass] + dConv )
        {
          nagree = 0;
          //          if (myrnd.flat() > 0.49999)
          //        {sgn = 1;} else {sgn=-1;}
          sgn = 1;
          TAxes[pass] = sgn*tpr/(tmax*thp);
          dThrust[pass] = thp;
        } // if thp
      nagree++;
    } // for iw (2)
    } // for pass ...
    // Find minor axis and value by orthogonality.
  if (myrnd.flat() > 0.49999)
    {sgn = 1;} else {sgn=-1;}
  TAxes[2].set( -sgn*TAxes[1].y(), sgn*TAxes[1].x(), 0);
  thp = 0.;
  for ( unsigned int i = 0; i < _partMom.size(); i++ )
    {
      thp += fabs(TAxes[2].dot(_partMom[i]) );
    } // for i 
  dThrust[2] = thp/tmax;

  // Rotate back to original coordinate system.
  for ( unsigned int i = 0;i < TAxes.size(); i++)
    {
      TAxes[i].rotateY(theta);
      TAxes[i].rotateZ(phi);
    }
  dOblateness = dThrust[1] - dThrust[2];

  _principleThrustValue = dThrust[0];
  _majorThrustValue     = dThrust[1];
  _minorThrustValue     = dThrust[2];
  _principleThrustAxis  =   TAxes[0];
  _majorThrustAxis      =   TAxes[1];
  _minorThrustAxis      =   TAxes[2];


  return  0;
}

//______________________________________________________________
// helper function to get sign of b
double Thrust::sign(double a, double b)
{
  if ( b < 0 )
    { return -fabs(a); } else { return fabs(a); }
}
//______________________________________________________________
double Thrust::min(double a, double b)
{
  if ( a < b )
    { return a; } else { return b; }
}
//______________________________________________________________


