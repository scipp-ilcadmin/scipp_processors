#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * TrackerOccupancyAnalysis.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
 */

#include "TrackerOccupancyAnalysis.h"
#include "scipp_ilc_utilities.h"
#include <iostream>
#include <algorithm>

#include <UTIL/ILDConf.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>
#include <TH3D.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;

TrackerOccupancyAnalysis TrackerOccupancyAnalysis;


static TFile* _rootfile;
static TH1D* _tileHits;
static TH1I* _layers;
static TH1D* _xPos;
static TH1D* _yPos;
static TH2D* _xyPos;
static TH3D* _xyzPos;
static TH1D* _eDep;
static TH1D* _cells;
static TH1F* _xmom;
static TH1F* _ymom;
static TH1F* _zmom;
static TH2F* _xymom, *_yzmom, *_xzmom;
static TH3F* _xyzmom;
static int _nEvt = 0;
static double pixelSize = 30.0; // Pixel size is in microns
static double tileSize; // Tile size also in microns
static double brlLayerArea[5] = {14817.6, 21168, 31752, 42336, 52920};
static int tileNumBrl[5];
int shingleError = 0;
int brlLyrN = 5;
int bunchCrossigns = 1;
bool aligned;
static vector< map<string,long> > BrlTiles;

static vector<int> layers;
static vector<int> subDets;
static vector<int> modules;
static vector<int> sensors;
static vector<int> sides;
static vector<float> eDepVals;
static vector<double> posxVals;
static vector<double> posyVals;
static vector<double> poszVals;
static vector<int> cell0Vals;
static vector<int> cell1Vals;
static vector<int> nmccontsVals;
static vector<int> tileIDVals;
static vector<float> xmomVals;
static vector<float> ymomVals;
static vector<float> zmomVals;
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

TrackerOccupancyAnalysis::TrackerOccupancyAnalysis() : Processor("TrackerOccupancyAnalysis") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void TrackerOccupancyAnalysis::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("TOA.root", "RECREATE");
  _tileHits = new TH1D("BrlEventsVCsHit", "BrlEventsVCsHit", 57, 997500, 1002000);
  _xPos = new TH1D("xposhits", "xposhits", 57, -62.0, 55.0);
  _yPos = new TH1D("yposhits", "yposhits", 57, -41.0, 48.0);
  _xyPos = new TH2D("xypos", "xypos", 57, -100.0, 100.0, 57, -100.0, 100.0);
  _xyzPos = new TH3D("xyzpos", "xyzpos", 57, -100.0, 100.0, 57, -100.0, 100.0, 57, -100.0, 100.0);
  _layers = new TH1I("layers", "layers", 57, -1, 7);
  _eDep = new TH1D("energy", "energy", 57, -1000, 1000);
  _xmom = new TH1F("xmom", "xmom", 57, -85.0, 85.0);
  _ymom = new TH1F("ymom", "ymom", 57, -85.0, 85.0);
  _zmom = new TH1F("zmom", "zmom", 57, -85.0, 85.0);
  _xymom = new TH2F("xymom", "xymom", 57, -85, 85, 57, -85, 85);
  _xyzmom = new TH3F("xyzmom", "xyzmom", 57, -85, 85, 57, -85, 85, 57, -85, 85);
  _nEvt = 0;
  tileSize = pixelSize/1000.0; // No clue why this is set like this
                               // It just is
                               // Deal with it
  for (int i=0; i<brlLyrN; ++i)
    {
      BrlTiles.push_back(map<string,long>());
    }
  for (int i=0; i<brlLyrN; ++i)
    {
      tileNumBrl[i] = (int)(brlLayerArea[i]/(tileSize * tileSize));
    }

}

void TrackerOccupancyAnalysis::processRunHeader( LCRunHeader* run)
{

}

void TrackerOccupancyAnalysis::processEvent( LCEvent * evt)
{
  _nEvt++;
  int check_layer = 0;
  int hit_count_limit = 100;
  bool use_limit = false;
  bool reject_negative = false;
  int hit_count = 0;
  double brlXseg = tileSize;
  double brlYseg = tileSize;
  int offsetX = 500000;
  int offsetY = 500000;
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  //LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHuts"); will do endcaps after barrelHits

  _nEvt++;
  for (int i = 0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( barrelHits );
      int layer = idDec( hit )[ILDCellID0::layer];
      layers.push_back(layer);
      int subdet = idDec( hit )[ILDCellID0::subdet];
      subDets.push_back(subdet);
      int module = idDec( hit )[ILDCellID0::module];
      modules.push_back(module);
      int sensor = idDec( hit )[ILDCellID0::sensor];
      sensors.push_back(sensor);
      int side = idDec( hit )[ILDCellID0::side];
      sides.push_back(side);
      double posx = hit->getPosition()[0]; // indecies for x,y,z components;
      posxVals.push_back(posx);
      double posy = hit->getPosition()[1];
      posyVals.push_back(posy);
      double posz = hit->getPosition()[2];
      poszVals.push_back(posz);
      int cell0 = hit->getCellID0();
      cell0Vals.push_back(cell0);
      int cell1 = hit->getCellID1();
      cell1Vals.push_back(cell1);
      float eDep = hit->getEDep();
      eDepVals.push_back(eDep);
      float xmom = hit->getMomentum()[0];
      xmomVals.push_back(xmom);
      float ymom = hit->getMomentum()[1];
      ymomVals.push_back(ymom);
      float zmom = hit->getMomentum()[2];
      zmomVals.push_back(zmom);
      int IDx = (int) ((posx/brlXseg) + 500000);
      int IDy = (int) ((posy/brlYseg) + 500000);
      int tileID = (module + IDx + IDy);
      tileIDVals.push_back(tileID);
      _tileHits->Fill(tileID);
      _xPos->Fill(posx);
      _yPos->Fill(posy);
      _xyPos->Fill(posx, posy);
      _xyzPos->Fill(posx, posy, posz);
      _layers->Fill(layer);
      _eDep->Fill(eDep);
      _xmom->Fill(xmom);
      _ymom->Fill(ymom);
      _zmom->Fill(zmom);
      _xymom->Fill(xmom, ymom);
      _xyzmom->Fill(xmom, ymom, zmom);
    }
  }


void TrackerOccupancyAnalysis::check( LCEvent * evt)
{

}

void TrackerOccupancyAnalysis::end()
{
  cout << " max energy: " << getMax(eDepVals) << " MinEnergy: " << getMin(eDepVals) << endl;
  cout << " max xmom: " << getMax(xmomVals) << " min xmom: " << getMin(xmomVals) << endl;
  cout << " max ymom: " << getMax(ymomVals) << " min ymom: " <<getMin(ymomVals) << endl;
  cout << " max zmom: " << getMax(xmomVals) << " min zmom: " <<getMin(zmomVals) << endl;
  _rootfile->Write();
}
