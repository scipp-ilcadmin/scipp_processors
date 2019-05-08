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
typedef vector<PixelGrid> Layers; // creating structures to store the pixels and their hits in
typedef vector<Layers> Barrel;
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
static Layers mods(12, PixelGrid(9.8*10, vector<int>(126*10,0)));    //Creates a 2D vector for each module, the size is variable depending on how many pixels we want
static Barrel barrel{};
static vector<double> zvals;                                     // Since these are rectangles it should be longer in one side than the other.
static vector<TH2D*> panels;

static double rsuba [12][3]; //vector to hold average values of the radius from the origin
static double thetavals [12][1]; //theta values for rotations, one for each module

template<typename T>
static T getMax(vector<T> &vec) //this gets the greatest value in a vector
{
  return *max_element(vec.begin(), vec.end());   
}
template<typename T>
static T getMin(vector<T> &vec)  //this one gets the smallest value
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
    _rootfile = new TFile("bread.root","RECREATE");  //Any histograms would be created in this block
    barrel.emplace_back(12, PixelGrid(9.8*10, vector<int>(126*10,0)));

    totes = new TH2D("totes", "magotes", 1000, -20, 20, 1000, -20, 20);
    //panels = new TH2D("panels", "fuckJAKE", 1000, -200, 200, 1000, -200, 200);
    for (int i = 1; i < 13; i++)
      {
	ogmods.push_back(new TH2D(Form("ogmod%d ", i), "mods", 1000, -20, 20, 1000, -20, 20));
	newmods.push_back(new TH2D(Form("newmod%d ", i), "mods", 1000, -20, 20, 1000, -20, 20));
	panels.push_back(new TH2D(Form("panels%d " , i), "panels", 200, -20, 20, 200, 0, 200));
      }

    _nEvt = 0 ;
}



void copyofexample::processRunHeader( LCRunHeader* run)
{ 
} 


void copyofexample::processEvent( LCEvent * evt ) { 
  LCCollection* barrelhits = evt->getCollection( "SiVertexBarrelHits" ); //Whatever we put in quotes is the type of collection we will be looking at
  _nEvt++;                                                               //e.g. SiVertexBarrelHits, SiVertexEndcaphits, MCParticle, etc...
    static const double xmin = 100;  //these are all set to ensure we can make all values positive
    static const double ymin = 100;  //this makes it easier to actually pixelize the detector
    static const double zmin = 100;  //can be modified if the geometry of the barrel or whatever is changed
    static const int step =1;
    for(int i=0; i < barrelhits->getNumberOfElements(); ++i)  //loop for each hit in the collection
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec( hit )[ILDCellID0::layer];   //info about hits stored in these types
	//int side = idDec ( hit )[ILDCellID0::side];
	int module = idDec (hit)[ILDCellID0::module];
	//int sensor = idDec (hit)[ILDCellID0::sensor];
	//int subdet = idDec (hit)[ILDCellID0::subdet];
	if (layer == 1)
	  {
	    double posx = hit->getPosition()[0]*1; //position stored in array, one for each coordinate
	    double posy = hit->getPosition()[1]*1;
	    double posz = hit->getPosition()[2]*1;
	    zvals.push_back(posz);
	    double radval = sqrt(posx*posx + posy*posy); //using this to try and eliminate superfluous hits, this can be improved upon.
	    if (radval > 12)
	      {
		totes->Fill(posx, posy);
		++nhit;
		int index = module;
		rsuba[index][0] += posx; //summing up hits to use in calculating our angular rotation later
		rsuba[index][1] += posy; //just need for x and y in this method, but hypothetically could do any combo of coords 
		rsuba[index][2] += posz;
		ogmods[index]->Fill(posx, posy); //plotting our original hits, 1 plot for each module
	      }
	  }
      }
    for (int j = 0; j < 12; ++j) //calculating angle for rotation, 1 for each module
      {
	rsuba[j][0] /= nhit;
	rsuba[j][1] /= nhit;
	thetavals[j][0] = atan2(rsuba[j][0], rsuba[j][1]); //simple geometry to determine the angle of rotatin
	//cout << "rsuba x: " << rsuba[j][0] << " rsuba y: " << rsuba[j][1];
	//cout << " theta for module " << j << ": " << thetavals[j][0] << endl; 
      }
    for (int i = 0; i < barrelhits->getNumberOfElements(); ++i) //needed a second loop to actually run through each hit and rotate to new coords
      {
	SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelhits->getElementAt(i));
	CellIDDecoder<SimTrackerHit> idDec(barrelhits);
	int layer = idDec ( hit ) [ILDCellID0::layer];
	int module = idDec (hit ) [ILDCellID0::module];
	if (layer == 1)
	  {
	    double posx = hit->getPosition()[0]*1;
	    double posy = hit->getPosition()[1]*1;
	    double posz = (hit->getPosition()[2] + zmin)*1;      //rotating the values so they exist in a flat plane in 2D make sit easier to pixelate
	    double radval = sqrt(posx*posx + posy*posy);     //I choose to analyze in xz dimensions, to maintain length in detectors
	    if (radval > 12)
	      {
		double newx = (posx * cos(thetavals[module][0]) - posy * sin(thetavals[module][0]))+6;  //newx and newy are rotated values
		double newy = (posy * cos(thetavals[module][0]) + posx * sin(thetavals[module][0]));        
		int index = module; // each hit has a module tag and I separate analysis based on module in this code
		panels[index]->Fill(newx,posz);
		newmods[index]->Fill(newx, newy);
		string id = to_string(static_cast<int>(newx/step)) + to_string(static_cast<int>(posz/step)); //create an 'id' for the hit to pick off unique and repeat hits
		modpixids[index].push_back(id);
		if(std::find(moduniqueids[index].begin(), moduniqueids[index].end(), id) == moduniqueids[index].end()) //find unique id vals, to just see what pixels are hit
		  {
		    moduniqueids[index].push_back(id);
		  }
		if ((newx >= 0 && posz >=0))
		  {
		    mods[index][static_cast<int>(newx/step)][static_cast<int>(posz/step)]++;
		    //cout << "new /step " << static_cast<int>(newx/step) << endl;
		  }
		else
		  {
		    cout << "ERROR,   newx:  " << newx << " , zval" << posz << endl; // another way to eliminate superfluous hits. can probably be improved
		  }
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

  double area = 9.8*126; //amount of pixels per module
  int numpix = static_cast<double>(area/.1);
  cout << "area: " << area << ", numpix: " << numpix << endl;
  double totality = 0; //counting hits in my roated analysis, should match nhit at the end of processing
  for (int i = 0; i < 12; i++) //final analysis for hits in rotated coords
    {
      cout << "total hits in module " << i+1 << ": " << modpixids[i].size() << endl;
      cout << "Pixels per module" << i+1 << ": " << numpix << endl;
      cout << "unique pixels hit in module " << i+1 << ": " << moduniqueids[i].size() << endl;
      double hitperc = static_cast<double>(moduniqueids[i].size()) / static_cast<double>(numpix) * 100;
      cout << "Percent of pixels hit in layer " << i+1 << ": " << hitperc << endl;
      cout << "Average percent of unique pixels hits per event: " << hitperc/numpix << endl;
      totality+=modpixids[i].size();
      cout << endl << endl << endl;
    }
  cout << " added mod vals for hits, should equal number below " << totality << endl;
  cout << "number of hits total: " << nhit << endl;
  //cout << getMax(zvals) << "    " << getMin(zvals) << endl;
  _rootfile->Write("Update");

}
