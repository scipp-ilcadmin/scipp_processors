#undef _GLIBCXX_USE_CXX11_ABI
#define _GLIBCXX_USE_CXX11_ABI 0

/*
 *
 * example.cc
 * @author Gregory Blockinger
 * July 6th, 2018
 *
 */

#include "example.h"
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
#include <THStack.h>
#include <TGraph2D.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"


using namespace lcio;
using namespace marlin;
using namespace std;

example example;

static const int bl1pix = 37044000;
static const int bl2pix = 52920000;
static const int bl3pix = 79380000;
static const int bl4pix = 105840000;
static const int bl5pix = 132300000;
static const int el1pix = 41230947;
static const int el2pix = 41179477;
static const int el3pix = 40695487;
static const int el4pix = 37756617;
static const double e1gap = 23.427/59;
static const double e2gap = 47.883/59;
static const double e3gap = 83.152/59;
static const double e4gap = 179.793/59;

static TFile* _rootfile;
static TH3D* threedim;
static TH2D* twodimexy;
static TH2D* twodimbxy;
static TH1D* allzvals;
static TH1D* etotzvals;
static TH1D* onedimbzvals;
static TH1D* onedimbphivals;
static TH1D* onedimeposrvals;
static TH1D* onedimenegrvals;
static TH1D* onedimepostvals;
static TH1D* onedimenegtvals;
static TH1D* onedimeposzvals;
static TH1D* onedimenegzvals;
static TH1D* bphigone[5];
static TH1D* bzgone[5];
static TH1D* ethetagonepos[4];
static TH1D* ergonepos[4];
static TH1D* ethetagoneneg[4];
static TH1D* ergoneneg[4];
static double modhits[5][30];
static TH2D* exy[4];
static int bposhits[5];
static TH1D* poseradvals[4];
static TH1D* negeradvals[4];
static TH1D* bocc[5];
static TH1D* btemp[5];
static TH1D* predbocc[5];
static TH1D* poseocc[4];
static TH1D* negeocc[4];
static TH1D* posetemp[4];
static TH1D* negetemp[4];
static TH1D* predeocc[4];
static TH1D* onedimeocc;
static TH1D* negbphigone[5];
static TH1D* predrvals[4];

example::example() : Processor("example") 
{
  _description = "Protype Processor";
  //registerProcessorParameter("RootOutputName", "output file", _root_file_name, std::string("output.root"));

}

void example::init()
{
  streamlog_out(DEBUG) << " init called " << endl;
  cout << "Initialized "  << endl;
  _rootfile = new TFile("example.root", "RECREATE");
  threedim = new TH3D("threedim", "3-D Model of Occupancy in Vertex Detector", 100, -80, 80, 100, -80, 80, 100, -200, 200);
  twodimbxy = new TH2D("twodimbxy", "2-D View of Back Scatter hits in the Barrel Array;x (mm); y (mm)", 1000, -80, 80, 1000, -80, 80);
  twodimexy = new TH2D("twodimexy", "2-D View of Back Scatter hits in the End-Cap Arrays;x (mm); y (mm)", 160, -80, 80, 160, -80, 80);
  allzvals = new TH1D("allzvals", "Z Distribution of hits; z (mm)", 200, -200, 200);
  onedimbzvals = new TH1D("onedimbzvals", "Z Component of Barrel Array Hits; (mm)", 200, -65, 65);
  onedimbphivals = new TH1D("onedimbphivals", "Angular Distribution of Barrel Array Hits", 200, -3.5, 3.5);
  onedimeposrvals = new TH1D("onedimeposrvals", "Positive End-Cap Array Radial Distribution of Hits; (mm) from origin", 80, 0, 80);
  onedimenegrvals = new TH1D("onedimenegrvals", "Negative End-Cap Array Radial Distribution of Hits; (mm) from origin", 80, 0, 80);
  onedimepostvals = new TH1D("onedimepostvals", "Positive End-Cap Array Angular Distribution of Hits; Angle in Radians", 200, -3.5, 3.5);
  onedimenegtvals = new TH1D("onedimenegtvals", "Negative End-Cap Array Angular Distribution of Hits; Angle in Radians", 200, -3.5, 3.5);
  onedimeposzvals = new TH1D("onedimeposzvals", "Positive End-Cap Array Z Component Distribution of Hits", 200, -10, 250);
  onedimenegzvals = new TH1D("onedimenegzvals", "Negative End-Cap Array Z Component Distribution of Hits", 200, -250, 10);
  onedimeocc = new TH1D("onedimeocc","BLERGH", 80, 0, 80);
  etotzvals = new TH1D("etotzvals", "End-Cap Array Z Component Distribution of Hits; (mm)", 200, -250, 250);
  for (int i=0; i < 5; ++i)
    {
      bphigone[i] = new TH1D(Form("bphigone%d", i+1), "Collapsed in Phi; z (mm)", 140, -70, 70);
      bocc[i] = new TH1D(Form("bocc%d", i+1), "Occupancy of Barrel Layers; z (mm)", 140, -70, 70);
      negbphigone[i] = new TH1D(Form("negbphigone%d", i+1), "Predicted Occupancy of Barrel Layers; z (mm)", 140, -70, 70);
      predbocc[i] = new TH1D(Form("predbocc%d", i+1), "Predicted Occupancy of Barrel Layer; z (mm)", 140, -70, 70);
      bzgone[i] = new TH1D(Form("bzgone%d", i+1), "Collapsed in Z; Angular Value of Hits", 470, -3.5, 3.5);
    }
  bphigone[0]->SetFillColor(kBlue);
  bphigone[0]->SetFillStyle(3001);
  bphigone[1]->SetFillColor(kCyan);
  bphigone[1]->SetFillStyle(3004);
  bphigone[2]->SetFillColor(kGreen);
  bphigone[2]->SetFillStyle(3007);
  bphigone[3]->SetFillColor(kOrange);
  bphigone[3]->SetFillStyle(3010);
  bphigone[4]->SetFillColor(kRed);
  bphigone[4]->SetFillStyle(3013);
  for (int i =0; i <4; ++i)
    {
      ethetagonepos[i] =  new TH1D(Form("ethetagonepos%d", i+1), "Collapsed in Theta; Z value of Hits", 100, -200, 200);
      ethetagoneneg[i] =  new TH1D(Form("ethetagoneneg%d", i+1), "Collapsed in Theta; Z value of Hits", 100, -200, 200);
      ergonepos[i] =  new TH1D(Form("ergonepos%d", i+1), "Collapsed in Radial Value; Angular Distribution of Hits", 500, -3.5, 3.5);      
      ergoneneg[i] =  new TH1D(Form("ergoneneg%d", i+1), "Collapsed in Radial Value; Angular Distribution of Hits", 500, -3.5, 3.5);
      exy[i] = new TH2D(Form("exy%d", i+1), "2-D View of End-Cap", 500, -80, 80, 500, -80, 80);
      poseradvals[i] = new TH1D(Form("poseradvals%d", i+1), "rad vals; r (mm)", 80, 0, 80);
      negeradvals[i] = new TH1D(Form("negeradvals%d", i+1), "rad vals; r (mm)", 80, 0, 80);
      poseocc[i] = new TH1D(Form("poseocc%d", i+1), "Occupancy of End-Cap Layers; r (mm)", 80, 0, 80);
      negeocc[i] = new TH1D(Form("negeocc%d", i+1), "Occupancy of End-Cap Layers; r (mm)", 80, 0, 80);
      posetemp[i] = new TH1D(Form("posetemp%d", i+1), "Occupancy of End-Cap Layers; r (mm)", 80, 0, 80);
      negetemp[i] = new TH1D(Form("negetemp%d", i+1), "Occupancy of End-Cap Layers; r (mm)", 80, 0, 80);
      predeocc[i] = new TH1D(Form("predeocc%d", i+1), "Predicted Bin Occupancy of End-Cap Layers; r (mm)", 80, 0, 80);
      predrvals[i] = new TH1D(Form("predrvals%d", i+1), "predicted amount of rvals", 80, 0, 80);
    }
  _nEvt = 0;
}

void example::processRunHeader( LCRunHeader* run)
{

}

void example::processEvent( LCEvent * evt)
{
  ++_nEvt;
  LCCollection* barrelHits = evt->getCollection("SiVertexBarrelHits");
  LCCollection* endcapHits = evt->getCollection("SiVertexEndcapHits");
  for (int i = 0; i < endcapHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(endcapHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec(endcapHits);
      int layer = idDec( hit )[ILDCellID0::layer];
      int module = idDec(hit)[ILDCellID0::module];
      double posx = hit->getPosition()[0];
      double posy = hit->getPosition()[1];
      double posz = hit->getPosition()[2];
      double radval = sqrt((posx)*(posx)+(posy)*(posy));
      double theta = atan2(posy,posx);
      allzvals->Fill(posz);
      etotzvals->Fill(posz);
      //threedim->Fill(posx, posy, posz);
      if (posz > 0)
	{
	  threedim->Fill(posx, posy, posz);
	  onedimeposrvals->Fill(radval);
	  onedimepostvals->Fill(theta);
	  onedimeposzvals->Fill(posz);
	  ethetagonepos[layer-1]->Fill(radval);
	  ergonepos[layer-1]->Fill(theta);
	  poseradvals[layer-1]->Fill(radval);
	  predrvals[layer-1]->Fill(radval);
	}
      if (posz < 0)
	{
	  threedim->Fill(posx, posy, posz);
	  onedimenegrvals->Fill(radval);
	  onedimenegtvals->Fill(theta);
	  onedimenegzvals->Fill(posz);
	  ethetagoneneg[layer-1]->Fill(radval);
	  ergoneneg[layer-1]->Fill(theta);
	  negeradvals[layer-1]->Fill(radval);
	  predrvals[layer-1]->Fill(radval);
	}
      if (radval > 12)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimexy->Fill(posx, posy);
	  exy[layer-1]->Fill(posx, posy);
	}
    }
  for (int i =0; i < barrelHits->getNumberOfElements(); ++i)
    {
      SimTrackerHit* hit = dynamic_cast<SimTrackerHit*>(barrelHits->getElementAt(i));
      CellIDDecoder<SimTrackerHit> idDec(endcapHits);
      int layer = idDec( hit )[ILDCellID0::layer];
      int module = idDec(hit)[ILDCellID0::module];
      double posx = hit->getPosition()[0];
      double posy = hit->getPosition()[1];
      double posz = hit->getPosition()[2];
      double radval = sqrt((posx)*(posx)+(posy)*(posy));
      double phi = atan2(posy,posx);
      threedim->Fill(posx,posy,posz);
      onedimbzvals->Fill(posz);
      onedimbphivals->Fill(phi);
      allzvals->Fill(posz);
      modhits[layer-1][module]++;
      if (radval > 13.0 && radval < 16.0)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	  bphigone[layer-1]->Fill(posz);
	  negbphigone[layer-1]->Fill(-posz);
	  negbphigone[layer-1]->Fill(posz);
	  bzgone[layer-1]->Fill(phi);
	  if (posz < 0)
	    {
	      bposhits[layer-1]++;
	    }
	}
      if (radval > 21.4 && radval < 24.99)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	  bphigone[layer-1]->Fill(posz);
          bzgone[layer-1]->Fill(phi);
	  negbphigone[layer-1]->Fill(-posz);
	  negbphigone[layer-1]->Fill(posz);
	  if (posz < 0)
	    {
	      bposhits[layer-1]++;
	    }
	}
      if (radval > 33.5 && radval < 37.75)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	  bphigone[layer-1]->Fill(posz);
          bzgone[layer-1]->Fill(phi);
	  negbphigone[layer-1]->Fill(-posz);
	  negbphigone[layer-1]->Fill(posz);
	  if (posz < 0)
	    {
	      bposhits[layer-1]++;
	    }
	}
      if (radval > 45.5 && radval < 50.3)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	  bphigone[layer-1]->Fill(posz);
          bzgone[layer-1]->Fill(phi);
 	  negbphigone[layer-1]->Fill(-posz);
	  negbphigone[layer-1]->Fill(posz);
	  if (posz < 0)
	    {
	      bposhits[layer-1]++;
	    }
	}
      if (radval > 58.0 && radval < 63.7)
	{
	  //threedim->Fill(posx, posy, posz);
	  twodimbxy->Fill(posx, posy);
	  bphigone[layer-1]->Fill(posz);
          bzgone[layer-1]->Fill(phi);
	  negbphigone[layer-1]->Fill(-posz);
	  negbphigone[layer-1]->Fill(posz);
	  if (posz < 0)
	    {
	      bposhits[layer-1]++;
	    }
	}

    }
}

void example::check( LCEvent * evt)
{

}

void example::end()
{
  
  for (int i = 1; i < 141; ++i)
    {
      bocc[0]->SetBinContent(i, (bphigone[0]->GetBinContent(i) / 9) / (bl1pix / 140));
      bocc[1]->SetBinContent(i, (bphigone[1]->GetBinContent(i) / 9) / (bl2pix / 140));
      bocc[2]->SetBinContent(i, (bphigone[2]->GetBinContent(i) / 9) / (bl3pix / 140));
      bocc[3]->SetBinContent(i, (bphigone[3]->GetBinContent(i) / 9) / (bl4pix / 140));
      bocc[4]->SetBinContent(i, (bphigone[4]->GetBinContent(i) / 9) / (bl5pix / 140));
      predbocc[0]->SetBinContent(i, (negbphigone[0]->GetBinContent(i) / 9) / (bl1pix / 140));
      predbocc[1]->SetBinContent(i, (negbphigone[1]->GetBinContent(i) / 9) / (bl2pix / 140));
      predbocc[2]->SetBinContent(i, (negbphigone[2]->GetBinContent(i) / 9) / (bl3pix / 140));
      predbocc[3]->SetBinContent(i, (negbphigone[3]->GetBinContent(i) / 9) / (bl4pix / 140));
      predbocc[4]->SetBinContent(i, (negbphigone[4]->GetBinContent(i) / 9) / (bl5pix / 140));

    }
  for (int i = 1; i < 81; ++i)
    {
      int big = (i+1)*(i+1);
      int little = i*i;
      int diff = (big-little);
      negeocc[0]->SetBinContent(i, ((negeradvals[0]->GetBinContent(i) / 9) * 0.0004) / ((M_PI*diff)-(16*e1gap)));
      negeocc[1]->SetBinContent(i, ((negeradvals[1]->GetBinContent(i) / 9) * 0.0004) / ((M_PI*diff)-(16*e2gap)));
      negeocc[2]->SetBinContent(i, ((negeradvals[2]->GetBinContent(i) / 9) * 0.0004) / ((M_PI*diff)-(16*e3gap)));
      negeocc[3]->SetBinContent(i, ((negeradvals[3]->GetBinContent(i) / 9) * 0.0004) / ((M_PI*diff)-(16*e4gap)));
      poseocc[0]->SetBinContent(i, ((poseradvals[0]->GetBinContent(i) / 9) * 0.0004) / ((M_PI*diff)-(16*e1gap)));
      poseocc[1]->SetBinContent(i, ((poseradvals[1]->GetBinContent(i) / 9) * 0.0004) / ((M_PI*diff)-(16*e2gap)));
      poseocc[2]->SetBinContent(i, ((poseradvals[2]->GetBinContent(i) / 9) * 0.0004) / ((M_PI*diff)-(16*e3gap)));
      poseocc[3]->SetBinContent(i, ((poseradvals[3]->GetBinContent(i) / 9) * 0.0004) / ((M_PI*diff)-(16*e4gap)));
      onedimeocc->SetBinContent(i, ((onedimenegrvals->GetBinContent(i) /9) * 0.0004) / ((M_PI*diff)-(16*e2gap)));
      predeocc[0]->SetBinContent(i, ((predrvals[0]->GetBinContent(i)/9) * 0.0004) / ((M_PI*diff)-(16*e1gap)));
      predeocc[1]->SetBinContent(i, ((predrvals[1]->GetBinContent(i)/9) * 0.0004) / ((M_PI*diff)-(16*e1gap)));
      predeocc[2]->SetBinContent(i, ((predrvals[2]->GetBinContent(i)/9) * 0.0004) / ((M_PI*diff)-(16*e1gap)));
      predeocc[3]->SetBinContent(i, ((predrvals[3]->GetBinContent(i)/9) * 0.0004) / ((M_PI*diff)-(16*e1gap)));
    }
  _rootfile->Write();
}
