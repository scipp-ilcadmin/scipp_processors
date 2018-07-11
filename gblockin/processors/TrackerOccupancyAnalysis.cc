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

#include <UTIL/ILDConf.h>
#include <EVENT/LCCollection.h>
#include <EVENT/SimCalorimeterHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/MCParticle.h>

#include <TFile.h>
#include <TH2D.h>

// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;

TrackerOccupancyAnalysis TrackerOccupancyAnalysis;

static TFile* _rootfile;
static int _nEvt = 0;
static double pixelSize = 5.0; // Pixel size is in microns
static double tileSize = 5.0; // Tile size also in microns
static double brlLayerArea[5] = {14817.6, 21168, 31752, 42336, 52920};
static int tileNumBrl[5];
int shingleError = 0;
int brlLyrN = 5;
int bunchCrossigns = 1;
bool aligned;


static vector<int> barrelLayers;
static vector<int> barrelSubDets;
static vector<int> barrelModules;
static vector<int> barrelSensors;
static vector<int> barrelSides;
static vector<float> barrelEnergyVals;
static vector<double> barrelPosVals;
static vector<int> barrelCell0Vals;
static vector<int> barrelCell1Vals;
static vector<int> barrelNmccontsVals;


TrackerOccupancyAnalysis::TrackerOccupancyAnalysis() : Processor("TrackerOccupancyAnalysis") 
{
  _description = "Protype Processor";
    registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void TrackerOccupancyAnalysis::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("TrackerOccupancyAnalysis.root", "RECREATE");
  _nEvt = 0;

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
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  //LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHuts"); will do endcaps after barrelHits

  _nEvt++;
  for (int i = 0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec( barrelHits );
      int layer = idDec( hit )[ILDCellID0::layer];
      barrelLayers.push_back(layer);
      int subdet = idDec( hit )[ILDCellID0::subdet];
      barrelSubDets.push_back(subdet);
      int module = idDec( hit )[ILDCellID0::module];
      barrelModules.push_back(module);
      int sensor = idDec( hit )[ILDCellID0::sensor];
      barrelSensors.push_back(sensor);
      int side = idDec( hit )[ILDCellID0::side];
      barrelSides.push_back(side);
      double pos = *hit->getPosition();
      barrelPosVals.push_back(pos);
      cout << i << "   layer: " << barrelLayers[i] << "   subdet: " << barrelSubDets[i] << "   ";
      cout << "   module: " << barrelModules[i] << " sensor: " << barrelSensors[i] << endl;
    }


}


void TrackerOccupancyAnalysis::check( LCEvent * evt)
{

}

void TrackerOccupancyAnalysis::end()
{
  
}
