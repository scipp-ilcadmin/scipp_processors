#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0
/* 
 * Ok, so I like C++11. Unfortunately,
 * Marlin is uilt with ansi C, so the processor
 * constructor freaks out about the string that is
 * passed to it as an argument. The above two lines
 * fix that issue, allowing our code to be compatible
 * with ansi C class declarations.
 * Big thanks to Daniel Bittman for helping me fix this.
 */

/*
 * author Gregory Blockinger
 * April 5, 2019
 */

#include "copyofexample.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>

#include <UTIL/ILDConf.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TGraph2D.h>


// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;


copyofexample copyofexample;

typedef vector<vector<int>> PixelGrid;
typedef vector<PixelGrid> Layers;                      // creating structures to store the pixels and their hits in
typedef vector<string> PixIDs;
typedef vector<PixIDs> layerpixIDs;

static TFile* _rootfile;
static int _nEvt=0; 
static int nhit = 0;
static TH2D* totes;
static vector<TH2D*> ogmods;
static vector<TH2D*> newmods;
static layerpixIDs modpixids(12, vector<string>());
static layerpixIDs moduniqueids(12, vector<string>());
static Layers mods(12, PixelGrid(1000, vector<int>(1000,0)));    //Creates a 2D vector for each module, the size is variable depending on how many pixels we want
static vector<double> zvals;                                     // Since these are rectangles it should be longer in one side than the other.
static TH2D* panels;

static double rsuba [12][3]; //vector to hold average values of the radius from the origin
static double thetavals [12][1]; //theta values for rotations, one for each module

template<typename T>
static T getMax(vector<T> &vec)
{
  return *max_element(vec.begin(), vec.end());
}
template<typename T>
static T getMin(vector<T> &vec)
{
  return *min_element(vec.begin(), vec.end());
}

copyofexample::copyofexample() : Processor("copyofexample") 
{
  // modify processor description
    _description = "Protype Processor" ;
}


void copyofexample::init()
{ 
    streamlog_out(DEBUG) << "   init called  " << std::endl ;
    _rootfile = new TFile("bread.root","RECREATE");
    totes = new TH2D("totes", "magotes", 1000, -20, 20, 1000, -20, 20);
    panels = new TH2D("panels", "fuckJAKE", 1000, -200, 200, 1000, -200, 200);
    for (int i = 1; i < 13; i++)
      {
	ogmods.push_back(new TH2D(Form("ogmod%d ", i), "mods", 1000, -20, 20, 1000, -20, 20));
	newmods.push_back(new TH2D(Form("newmod%d ", i), "mods", 1000, -20, 20, 1000, -20, 20));
      }

    _nEvt = 0 ;
}



void copyofexample::processRunHeader( LCRunHeader* run)
{ 
} 


void copyofexample::processEvent( LCEvent * evt ) { 
    LCCollection* barrelhits = evt->getCollection( "SiVertexBarrelHits" );
    _nEvt++;
    static const double xmin = 100;
    static const double ymin = 100;
    static const double zmin = 100;
    static const int step =1;
    for(int i=0; i < barrelhits->getNumberOfElements(); ++i)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec( hit )[ILDCellID0::layer];
	//int side = idDec ( hit )[ILDCellID0::side];
	int module = idDec (hit)[ILDCellID0::module];
	//int sensor = idDec (hit)[ILDCellID0::sensor];
	//int subdet = idDec (hit)[ILDCellID0::subdet];
	if (layer == 1)
	  {
	    double posx = hit->getPosition()[0];
	    double posy = hit->getPosition()[1];
	    double posz = hit->getPosition()[2];
	    zvals.push_back(posz);
	    double radval = sqrt(posx*posx + posy*posy);
	    if (radval > 12)
	      {
		totes->Fill(posx, posy);
		++nhit;
		//[&] ()
		//{
		    int index = module;
		    rsuba[index][0] += posx;
		    rsuba[index][1] += posy;
		    rsuba[index][2] += posz;
		    ogmods[index]->Fill(posx, posy);
		    //}();
	      }
	  }
      }
    for (int j = 0; j < 12; ++j)
      {
	rsuba[j][0] /= nhit;
	rsuba[j][1] /= nhit;
	thetavals[j][0] = atan2(rsuba[j][0], rsuba[j][1]);
	//cout << "rsuba x: " << rsuba[j][0] << " rsuba y: " << rsuba[j][1];
	//cout << " theta for module " << j << ": " << thetavals[j][0] << endl; 
      }
    for (int i = 0; i < barrelhits->getNumberOfElements(); ++i)
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec ( hit ) [ILDCellID0::layer];
	int module = idDec (hit ) [ILDCellID0::module];
	if (layer == 1)
	  {
	    double posx = hit->getPosition()[0];
	    double posy = hit->getPosition()[1];
	    double posz = hit->getPosition()[2] + zmin;
	    double radval = sqrt(posx*posx + posy*posy);
	    if (radval > 12)
	      {
		double newx = (posx * cos(thetavals[module][0]) - posy * sin(thetavals[module][0])) + 6.0;
		double newy = (posy * cos(thetavals[module][0]) + posx * sin(thetavals[module][0]));
		panels->Fill(newx,posz);
		//[&]()
		//{
		    int index = module;
		    newmods[index]->Fill(newx, newy);
		    string id = to_string(newx/step) + to_string(posz/step);
		    modpixids[index].push_back(id);
		    if(std::find(moduniqueids[index].begin(), moduniqueids[index].end(), id) == moduniqueids[index].end())
		      {
			moduniqueids[index].push_back(id);
		      }
		    if ((newx >= 0 && posz >=0))
		      {
			mods[index][newx/step][posz/step]++;
		      }
		      else
			{
			  cout << "ERROR: " << newx << " , " << posz << endl;
			}
		    // }();
	      }
	  }
      }
}


void copyofexample::check( LCEvent * evt )
{ 
    // nothing to check here - could be used to fill checkplots in reconstruction processor
}



void copyofexample::end()
{

  int matters = 5*63*2;
  double totality = 0;
  for (int i = 0; i < 12; i++)
    {
      cout << "total hits in module " << i+1 << ": " << modpixids[i].size() << endl;
      cout << "area of modules (mm) " << i+1 << ": " << matters << endl;
      cout << "unique pixels hit in module " << i+1 << ": " << moduniqueids[i].size() << endl;
      double hitperc = static_cast<double>(moduniqueids[i].size()) / static_cast<double>(matters) * 100;
      cout << "Percent of pixels hit in layer " << i+1 << ": " << hitperc << endl;
      cout << "Average percent of unique pixels hits per event: " << hitperc/matters << endl;
      totality+=modpixids[i].size();
      cout << endl << endl << endl;
    }
  cout << " added mod vals for hits, should equal number below " << totality << endl;
  cout << "number of hits total: " << nhit << endl;
  //cout << getMax(zvals) << "    " << getMin(zvals) << endl;
  _rootfile->Write("Update");

}
