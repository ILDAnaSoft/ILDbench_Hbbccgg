#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <string>
#include <fstream>
#include "../myUtil_new.hh"
#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TPaveText.h"
#include "TPaveLabel.h"
#include "TLegend.h"
#include "TFrame.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

//Double_t x[220];

using namespace std;
using namespace TMVA;

//forward declaration
double calBjet(TLorentzVector jj, jetdata data);

int makeplot_and_rootfile_m08p03(){
//Int_t main(){
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.0);

  //luminosity and correction factor for polarization
  const double lumi = 1600.0;
  const double eplus = 0.1;
  const double eminus = 0.9;
  const double pplus = 0.65;
  const double pminus = 0.35;
  const double bplus = 0.5;
  const double bminus = 0.5;
  const double wplus = 0.5;
  const double wminus = 0.5;

  jetdata data;

  //define histograms
  const int cats = 9;
  TH1F *mh[cats], *cosjj[cats], *pt[cats];
  TH1F *j1e[cats], *j2e[cats], *j1cos[cats], *j2cos[cats], *j1m[cats], *j2m[cats];
  TH1F *j1Bjet[cats], *j2Bjet[cats], *j1nchg[cats], *j2nchg[cats];
  TH2F *flavortag[cats];
  TH1F *npfo[cats], *evis[cats], *missmass[cats];  //without beam backgrounds
  TH1F *thrust0[cats], *thrust1[cats], *thrust2[cats];  //without beam backgrounds
  TH1F *y23[cats], *y34[cats], *y45[cats], *y56[cats];  
  TH1F *wzmass1[cats], *wzmass2[cats];  
  TH1F *mva0[cats], *mva1[cats], *mva2[cats]; 
  TH3F *temp[cats];

  string histtxt,labeltxt;
  for(int i=0;i<cats;i++){
    //mh
    histtxt = "mh_" + itos(i);
    labeltxt = "m(H) (GeV/c^{2})";
    //mh[i] = new TH1F(histtxt.c_str(),"",60,0.0,300.0);
    mh[i] = new TH1F(histtxt.c_str(),"",60,90.0,150.0);
    mh[i]->SetMarkerStyle(20);
    mh[i]->SetFillColor(0);
    if(i<3) mh[i]->SetLineColor(2 + i);
    else mh[i]->SetLineColor(i*10 + 1);
    //setColor(mh[i], j);
    mh[i]->SetLineWidth(1.0);
    mh[i]->GetYaxis()->SetTitle("Number of Events /10GeV/c^{2}");
    mh[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    mh[i]->SetMinimum(0.0);

    //cosjj
    histtxt = "cosjj_" + itos(i);
    labeltxt = "cos#theta_{jj}";
    cosjj[i] = new TH1F(histtxt.c_str(),"",40,-1.0,1.0);
    cosjj[i]->SetMarkerStyle(20);
    cosjj[i]->SetFillColor(0);
    if(i<3) cosjj[i]->SetLineColor(2 + i);
    else cosjj[i]->SetLineColor(i*10 + 1);
    //setColor(cosjj[i], j);
    cosjj[i]->SetLineWidth(1.0);
    cosjj[i]->GetYaxis()->SetTitle("Number of Events /0.05");
    cosjj[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    cosjj[i]->SetMinimum(0.0);

    //pt
    histtxt = "pt_" + itos(i);
    labeltxt = "Pt";
    pt[i] = new TH1F(histtxt.c_str(),"",60,0.0,300.0);
    pt[i]->SetMarkerStyle(20);
    pt[i]->SetFillColor(0);
    if(i<3) pt[i]->SetLineColor(2 + i);
    else pt[i]->SetLineColor(i*10 + 1);
    //setColor(pt[i], j);
    pt[i]->SetLineWidth(1.0);
    pt[i]->GetYaxis()->SetTitle("Number of Events /5GeV/c");
    pt[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    pt[i]->SetMinimum(0.0);

    //j1e
    histtxt = "j1e_" + itos(i);
    labeltxt = "E_{j1}";
    j1e[i] = new TH1F(histtxt.c_str(),"",60,0.0,300.0);
    j1e[i]->SetMarkerStyle(20);
    j1e[i]->SetFillColor(0);
    if(i<3) j1e[i]->SetLineColor(2 + i);
    else j1e[i]->SetLineColor(i*10 + 1);
    //setColor(j1e[i], j);
    j1e[i]->SetLineWidth(1.0);
    j1e[i]->GetYaxis()->SetTitle("Number of Events /5GeV");
    j1e[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j1e[i]->SetMinimum(0.0);

    //j2e
    histtxt = "j2e_" + itos(i);
    labeltxt = "E_{j2}";
    j2e[i] = new TH1F(histtxt.c_str(),"",60,0.0,300.0);
    j2e[i]->SetMarkerStyle(20);
    j2e[i]->SetFillColor(0);
    if(i<3) j2e[i]->SetLineColor(2 + i);
    else j2e[i]->SetLineColor(i*10 + 1);
    //setColor(j2e[i], j);
    j2e[i]->SetLineWidth(1.0);
    j2e[i]->GetYaxis()->SetTitle("Number of Events /5GeV");
    j2e[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j2e[i]->SetMinimum(0.0);

    //j1cos
    histtxt = "j1cos_" + itos(i);
    labeltxt = "cos#theta_{j1}";
    j1cos[i] = new TH1F(histtxt.c_str(),"",40,-1.0,1.0);
    j1cos[i]->SetMarkerStyle(20);
    j1cos[i]->SetFillColor(0);
    if(i<3) j1cos[i]->SetLineColor(2 + i);
    else j1cos[i]->SetLineColor(i*10 + 1);
    //setColor(j1cos[i], j);
    j1cos[i]->SetLineWidth(1.0);
    j1cos[i]->GetYaxis()->SetTitle("Number of Events /0.05");
    j1cos[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j1cos[i]->SetMinimum(0.0);

    //j2cos
    histtxt = "j2cos_" + itos(i);
    labeltxt = "cos#theta_{j2}";
    j2cos[i] = new TH1F(histtxt.c_str(),"",40,-1.0,1.0);
    j2cos[i]->SetMarkerStyle(20);
    j2cos[i]->SetFillColor(0);
    if(i<3) j2cos[i]->SetLineColor(2 + i);
    else j2cos[i]->SetLineColor(i*10 + 1);
    //setColor(j2cos[i], j);
    j2cos[i]->SetLineWidth(1.0);
    j2cos[i]->GetYaxis()->SetTitle("Number of Events /0.05");
    j2cos[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j2cos[i]->SetMinimum(0.0);

    //j1m
    histtxt = "j1m_" + itos(i);
    labeltxt = "m_{j1} (GeV/c^{2})";
    j1m[i] = new TH1F(histtxt.c_str(),"",50,0.0,150.0);
    j1m[i]->SetMarkerStyle(20);
    j1m[i]->SetFillColor(0);
    if(i<3) j1m[i]->SetLineColor(2 + i);
    else j1m[i]->SetLineColor(i*10 + 1);
    //setColor(j1cos[i], j);
    j1m[i]->SetLineWidth(1.0);
    j1m[i]->GetYaxis()->SetTitle("Number of Events /3.0GeV/c^{2}");
    j1m[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j1m[i]->SetMinimum(0.0);

    //j2m
    histtxt = "j2m_" + itos(i);
    labeltxt = "m_{j2} (GeV/c^{2})";
    j2m[i] = new TH1F(histtxt.c_str(),"",50,0.0,150.0);
    j2m[i]->SetMarkerStyle(20);
    j2m[i]->SetFillColor(0);
    if(i<3) j2m[i]->SetLineColor(2 + i);
    else j2m[i]->SetLineColor(i*10 + 1);
    //setColor(j2cos[i], j);
    j2m[i]->SetLineWidth(1.0);
    j2m[i]->GetYaxis()->SetTitle("Number of Events /3.0GeV/c^{2}");
    j2m[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j2m[i]->SetMinimum(0.0);

    //j1Bjet
    histtxt = "j1Bjet_" + itos(i);
    labeltxt = "Bjet_{j1}";
    j1Bjet[i] = new TH1F(histtxt.c_str(),"",50,0.0,1.0);
    j1Bjet[i]->SetMarkerStyle(20);
    j1Bjet[i]->SetFillColor(0);
    if(i<3) j1Bjet[i]->SetLineColor(2 + i);
    else j1Bjet[i]->SetLineColor(i*10 + 1);
    //setColor(j1cos[i], j);
    j1Bjet[i]->SetLineWidth(1.0);
    j1Bjet[i]->GetYaxis()->SetTitle("Number of Events /0.02");
    j1Bjet[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j1Bjet[i]->SetMinimum(0.0);

    //j2Bjet
    histtxt = "j2Bjet_" + itos(i);
    labeltxt = "Bjet_{j2}";
    j2Bjet[i] = new TH1F(histtxt.c_str(),"",50,0.0,1.0);
    j2Bjet[i]->SetMarkerStyle(20);
    j2Bjet[i]->SetFillColor(0);
    if(i<3) j2Bjet[i]->SetLineColor(2 + i);
    else j2Bjet[i]->SetLineColor(i*10 + 1);
    //setColor(j2cos[i], j);
    j2Bjet[i]->SetLineWidth(1.0);
    j2Bjet[i]->GetYaxis()->SetTitle("Number of Events /0.02");
    j2Bjet[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j2Bjet[i]->SetMinimum(0.0);

    //j1nchg
    histtxt = "j1nchg_" + itos(i);
    labeltxt = "Nchg_{j1}";
    j1nchg[i] = new TH1F(histtxt.c_str(),"",70,0.0,70.0);
    j1nchg[i]->SetMarkerStyle(20);
    j1nchg[i]->SetFillColor(0);
    if(i<3) j1nchg[i]->SetLineColor(2 + i);
    else j1nchg[i]->SetLineColor(i*10 + 1);
    //setColor(j1cos[i], j);
    j1nchg[i]->SetLineWidth(1.0);
    j1nchg[i]->GetYaxis()->SetTitle("Number of Events /1");
    j1nchg[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j1nchg[i]->SetMinimum(0.0);

    //j2nchg
    histtxt = "j2nchg_" + itos(i);
    labeltxt = "Nchg_{j2}";
    j2nchg[i] = new TH1F(histtxt.c_str(),"",70,0.0,70.0);
    j2nchg[i]->SetMarkerStyle(20);
    j2nchg[i]->SetFillColor(0);
    if(i<3) j2nchg[i]->SetLineColor(2 + i);
    else j2nchg[i]->SetLineColor(i*10 + 1);
    //setColor(j2cos[i], j);
    j2nchg[i]->SetLineWidth(1.0);
    j2nchg[i]->GetYaxis()->SetTitle("Number of Events /1");
    j2nchg[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    j2nchg[i]->SetMinimum(0.0);

    //npfo
    histtxt = "npfo_" + itos(i);
    labeltxt = "N_{pfo}";
    npfo[i] = new TH1F(histtxt.c_str(),"",100,0.0,100.0);
    npfo[i]->SetMarkerStyle(20);
    npfo[i]->SetFillColor(0);
    if(i<3) npfo[i]->SetLineColor(2 + i);
    else npfo[i]->SetLineColor(i*10 + 1);
    //setColor(npfo[i], j);
    npfo[i]->SetLineWidth(1.0);
    npfo[i]->GetYaxis()->SetTitle("Number of Events /1");
    npfo[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    npfo[i]->SetMinimum(0.0);

    //evis
    histtxt = "evis_" + itos(i);
    labeltxt = "E_{vis}";
    evis[i] = new TH1F(histtxt.c_str(),"",50,0.0,500.0);
    evis[i]->SetMarkerStyle(20);
    evis[i]->SetFillColor(0);
    if(i<3) evis[i]->SetLineColor(2 + i);
    else evis[i]->SetLineColor(i*10 + 1);
    //setColor(evis[i], j);
    evis[i]->SetLineWidth(1.0);
    evis[i]->GetYaxis()->SetTitle("Number of Events /10GeV");
    evis[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    evis[i]->SetMinimum(0.0);

    //missmass
    histtxt = "missmass_" + itos(i);
    labeltxt = "m(miss) (GeV/c^{2})";
    missmass[i] = new TH1F(histtxt.c_str(),"",50,0.0,500.0);
    missmass[i]->SetMarkerStyle(20);
    missmass[i]->SetFillColor(0);
    if(i<3) missmass[i]->SetLineColor(2 + i);
    else missmass[i]->SetLineColor(i*10 + 1);
    //setColor(missmass[i], j);
    missmass[i]->SetLineWidth(1.0);
    missmass[i]->GetYaxis()->SetTitle("Number of Events /10GeV/c^{2}");
    missmass[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    missmass[i]->SetMinimum(0.0);

    //flavortag
    histtxt = "flavortag_" + itos(i);
    labeltxt = "b-likeliness";
    flavortag[i] = new TH2F(histtxt.c_str(),"",10,0.0,1.0,10,0.0,1.0);
    flavortag[i]->SetMarkerStyle(20);
    flavortag[i]->SetFillColor(0);
    //if(i<3) flavortag[i]->SetLineColor(2 + i);
    //else flavortag[i]->SetLineColor(i*10 + 1);
    //setColor(flavortag[i], j);
    flavortag[i]->SetLineWidth(1.0);
    flavortag[i]->GetYaxis()->SetTitle("c-likeliness");
    flavortag[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    //flavortag[i]->SetMinimum(0.0);

    //temp
    histtxt = "template_" + itos(i);
    labeltxt = "b-likeliness";
    temp[i] = new TH3F(histtxt.c_str(),"",50,0.0,1.0,50,0.0,1.0,50,0.0,1.0);
    //temp[i]->SetMarkerStyle(20);
    //temp[i]->SetFillColor(0);
    //if(i<3) temp[i]->SetLineColor(2 + i);
    //else temp[i]->SetLineColor(i*10 + 1);
    //setColor(temp[i], j);
    //temp[i]->SetLineWidth(1.0);
    temp[i]->GetZaxis()->SetTitle("bc-likeliness");
    temp[i]->GetYaxis()->SetTitle("c-likeliness");
    temp[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    //temp[i]->SetMinimum(0.0);

    //thrust0
    histtxt = "thrust0_" + itos(i);
    labeltxt = "Principal Thrust";
    thrust0[i] = new TH1F(histtxt.c_str(),"",50,0.0,1.0);
    thrust0[i]->SetMarkerStyle(20);
    thrust0[i]->SetFillColor(0);
    if(i<3) thrust0[i]->SetLineColor(2 + i);
    else thrust0[i]->SetLineColor(i*10 + 1);
    //setColor(thrust0[i], j);
    thrust0[i]->SetLineWidth(1.0);
    thrust0[i]->GetYaxis()->SetTitle("Number of Events /0.02");
    thrust0[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    thrust0[i]->SetMinimum(0.0);

    //thrust1
    histtxt = "thrust1_" + itos(i);
    labeltxt = "Major Thrust";
    thrust1[i] = new TH1F(histtxt.c_str(),"",50,0.0,1.0);
    thrust1[i]->SetMarkerStyle(20);
    thrust1[i]->SetFillColor(0);
    if(i<3) thrust1[i]->SetLineColor(2 + i);
    else thrust1[i]->SetLineColor(i*10 + 1);
    //setColor(thrust1[i], j);
    thrust1[i]->SetLineWidth(1.0);
    thrust1[i]->GetYaxis()->SetTitle("Number of Events /0.02");
    thrust1[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    thrust1[i]->SetMinimum(0.0);

    //thrust2
    histtxt = "thrust2_" + itos(i);
    labeltxt = "Minor Thrust";
    thrust2[i] = new TH1F(histtxt.c_str(),"",50,0.0,1.0);
    thrust2[i]->SetMarkerStyle(20);
    thrust2[i]->SetFillColor(0);
    if(i<3) thrust2[i]->SetLineColor(2 + i);
    else thrust2[i]->SetLineColor(i*10 + 1);
    //setColor(thrust2[i], j);
    thrust2[i]->SetLineWidth(1.0);
    thrust2[i]->GetYaxis()->SetTitle("Number of Events /0.02");
    thrust2[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    thrust2[i]->SetMinimum(0.0);

    //y23
    histtxt = "y23_" + itos(i);
    labeltxt = "Log(y23)";
    y23[i] = new TH1F(histtxt.c_str(),"",150,-15.0,0.0);
    y23[i]->SetMarkerStyle(20);
    y23[i]->SetFillColor(0);
    if(i<3) y23[i]->SetLineColor(2 + i);
    else y23[i]->SetLineColor(i*10 + 1);
    //setColor(y23[i], j);
    y23[i]->SetLineWidth(1.0);
    y23[i]->GetYaxis()->SetTitle("Number of Events /0.1");
    y23[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    y23[i]->SetMinimum(0.0);

    //y34
    histtxt = "y34_" + itos(i);
    labeltxt = "Log(y34)";
    y34[i] = new TH1F(histtxt.c_str(),"",150,-15.0,0.0);
    y34[i]->SetMarkerStyle(20);
    y34[i]->SetFillColor(0);
    if(i<3) y34[i]->SetLineColor(2 + i);
    else y34[i]->SetLineColor(i*10 + 1);
    //setColor(y34[i], j);
    y34[i]->SetLineWidth(1.0);
    y34[i]->GetYaxis()->SetTitle("Number of Events /0.1");
    y34[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    y34[i]->SetMinimum(0.0);

    //y45
    histtxt = "y45_" + itos(i);
    labeltxt = "Log(y45)";
    y45[i] = new TH1F(histtxt.c_str(),"",150,-15.0,0.0);
    y45[i]->SetMarkerStyle(20);
    y45[i]->SetFillColor(0);
    if(i<3) y45[i]->SetLineColor(2 + i);
    else y45[i]->SetLineColor(i*10 + 1);
    //setColor(y45[i], j);
    y45[i]->SetLineWidth(1.0);
    y45[i]->GetYaxis()->SetTitle("Number of Events /0.1");
    y45[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    y45[i]->SetMinimum(0.0);

    //y56
    histtxt = "y56_" + itos(i);
    labeltxt = "Log(y56)";
    y56[i] = new TH1F(histtxt.c_str(),"",150,-15.0,0.0);
    y56[i]->SetMarkerStyle(20);
    y56[i]->SetFillColor(0);
    if(i<3) y56[i]->SetLineColor(2 + i);
    else y56[i]->SetLineColor(i*10 + 1);
    //setColor(y56[i], j);
    y56[i]->SetLineWidth(1.0);
    y56[i]->GetYaxis()->SetTitle("Number of Events /0.1");
    y56[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    y56[i]->SetMinimum(0.0);

    //wzmass1
    histtxt = "wzmass1_" + itos(i);
    labeltxt = "m_{2jets1} (GeV/c^{2})";
    wzmass1[i] = new TH1F(histtxt.c_str(),"",50,0.0,200.0);
    wzmass1[i]->SetMarkerStyle(20);
    wzmass1[i]->SetFillColor(0);
    if(i<3) wzmass1[i]->SetLineColor(2 + i);
    else wzmass1[i]->SetLineColor(i*10 + 1);
    //setColor(y45[i], j);
    wzmass1[i]->SetLineWidth(1.0);
    wzmass1[i]->GetYaxis()->SetTitle("Number of Events /4GeV/c^{2}");
    wzmass1[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    wzmass1[i]->SetMinimum(0.0);

    //wzmass2
    histtxt = "wzmass2_" + itos(i);
    labeltxt = "m_{2jets2} (GeV/c^{2})";
    wzmass2[i] = new TH1F(histtxt.c_str(),"",50,0.0,200.0);
    wzmass2[i]->SetMarkerStyle(20);
    wzmass2[i]->SetFillColor(0);
    if(i<3) wzmass2[i]->SetLineColor(2 + i);
    else wzmass2[i]->SetLineColor(i*10 + 1);
    //setColor(y45[i], j);
    wzmass2[i]->SetLineWidth(1.0);
    wzmass2[i]->GetYaxis()->SetTitle("Number of Events /4GeV/c^{2}");
    wzmass2[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    wzmass2[i]->SetMinimum(0.0);

    //mva0
    histtxt = "mva0_" + itos(i);
    labeltxt = "MVA_{bb}";
    mva0[i] = new TH1F(histtxt.c_str(),"",100,-1.0,1.0);
    mva0[i]->SetMarkerStyle(20);
    mva0[i]->SetFillColor(0);
    if(i<3) mva0[i]->SetLineColor(2 + i);
    else mva0[i]->SetLineColor(i*10 + 1);
    //setColor(mva0[i], j);
    mva0[i]->SetLineWidth(1.0);
    mva0[i]->GetYaxis()->SetTitle("Number of Events /0.02");
    mva0[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    mva0[i]->SetMinimum(0.0);

    //mva1
    histtxt = "mva1_" + itos(i);
    labeltxt = "MVA_{cc}";
    mva1[i] = new TH1F(histtxt.c_str(),"",100,-1.0,1.0);
    mva1[i]->SetMarkerStyle(20);
    mva1[i]->SetFillColor(0);
    if(i<3) mva1[i]->SetLineColor(2 + i);
    else mva1[i]->SetLineColor(i*10 + 1);
    //setColor(mva1[i], j);
    mva1[i]->SetLineWidth(1.0);
    mva1[i]->GetYaxis()->SetTitle("Number of Events /0.02");
    mva1[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    mva1[i]->SetMinimum(0.0);

    //mva2
    histtxt = "mva2_" + itos(i);
    labeltxt = "MVA_{gg}";
    mva2[i] = new TH1F(histtxt.c_str(),"",100,-1.0,1.0);
    mva2[i]->SetMarkerStyle(20);
    mva2[i]->SetFillColor(0);
    if(i<3) mva2[i]->SetLineColor(2 + i);
    else mva2[i]->SetLineColor(i*10 + 1);
    //setColor(mva2[i], j);
    mva2[i]->SetLineWidth(1.0);
    mva2[i]->GetYaxis()->SetTitle("Number of Events /0.02");
    mva2[i]->GetXaxis()->SetTitle(labeltxt.c_str());
    mva2[i]->SetMinimum(0.0);
  }

  cout << "histograms done" << endl;

  //define rootfile
  double vars[100];
  TFile *tf = new TFile("./vars_for_tmva_m08p03test3.root","RECREATE");
  TTree *tt[cats+4];
  tt[0] = new TTree("nnh_bb","H->bb");
  tt[1] = new TTree("nnh_cc","H->cc");
  tt[2] = new TTree("nnh_gg","H->gg");
  tt[3] = new TTree("2f","2 fermion");
  tt[4] = new TTree("4f","4 fermion");
  tt[5] = new TTree("5f","5 fermion");
  tt[6] = new TTree("6f","6 fermion");
  tt[7] = new TTree("aa","gg->hadrons");
  tt[8] = new TTree("higgs_other","other higgs processes");
  tt[9] = new TTree("background","all the backgrounds");
  tt[10] = new TTree("signal","signal");
  tt[11] = new TTree("zh_other","other zh processes");
  tt[12] = new TTree("nnh_other","other nnh processes");

  for(int i=0;i<cats+4;i++){
    tt[i]->Branch("mh",&vars[0],"mh /D");
    tt[i]->Branch("j1px",&vars[1],"j1px /D");
    tt[i]->Branch("j1py",&vars[2],"j1py /D");
    tt[i]->Branch("j1pz",&vars[3],"j1pz /D");
    tt[i]->Branch("j1e",&vars[4],"j1e /D");
    tt[i]->Branch("j1cos",&vars[5],"j1cos /D");
    tt[i]->Branch("j1btag",&vars[6],"j1btag /D");
    tt[i]->Branch("j1ctag",&vars[7],"j1ctag /D");
    tt[i]->Branch("j1nchg",&vars[8],"j1nchg /D");
    tt[i]->Branch("j1Bjet",&vars[9],"j1Bjet /D");
    tt[i]->Branch("j1m",&vars[10],"j1m /D");
    tt[i]->Branch("j2px",&vars[11],"j2px /D");
    tt[i]->Branch("j2py",&vars[12],"j2py /D");
    tt[i]->Branch("j2pz",&vars[13],"j2pz /D");
    tt[i]->Branch("j2e",&vars[14],"j2e /D");
    tt[i]->Branch("j2cos",&vars[15],"j2cos /D");
    tt[i]->Branch("j2btag",&vars[16],"j2btag /D");
    tt[i]->Branch("j2ctag",&vars[17],"j2ctag /D");
    tt[i]->Branch("j2nchg",&vars[18],"j2nchg /D");
    tt[i]->Branch("j2Bjet",&vars[19],"j2Bjet /D");
    tt[i]->Branch("j2m",&vars[20],"j2m /D");
    tt[i]->Branch("cosjj",&vars[21],"cosjj /D");
    tt[i]->Branch("evis",&vars[22],"evis /D");
    tt[i]->Branch("eviscorr",&vars[23],"eviscorr /D");
    tt[i]->Branch("npfo",&vars[24],"npfo /D");
    tt[i]->Branch("metpx",&vars[25],"metpx /D");
    tt[i]->Branch("metpy",&vars[26],"metpy /D");
    tt[i]->Branch("metpz",&vars[27],"metpz /D");
    tt[i]->Branch("mete",&vars[28],"mete /D");
    tt[i]->Branch("metecorr",&vars[29],"metcorr /D");
    tt[i]->Branch("missmass",&vars[30],"missmass /D");
    tt[i]->Branch("cbtag",&vars[31],"cbtag /D");
    tt[i]->Branch("cctag",&vars[32],"cctag /D");
    tt[i]->Branch("cbctag",&vars[33],"cbctag /D");
    tt[i]->Branch("thrust0",&vars[34],"thrust0 /D");
    tt[i]->Branch("thrust1",&vars[35],"thrust1 /D");
    tt[i]->Branch("thrust2",&vars[36],"thrust2 /D");
    tt[i]->Branch("wzmass1",&vars[37],"wzmass1 /D");
    tt[i]->Branch("wzmass2",&vars[38],"wzmass2 /D");
    tt[i]->Branch("wzmass3",&vars[39],"wzmass3 /D");
    tt[i]->Branch("wzmass4",&vars[40],"wzmass4 /D");
    tt[i]->Branch("y23",&vars[41],"y23 /D");
    tt[i]->Branch("y34",&vars[42],"y34 /D");
    tt[i]->Branch("y45",&vars[43],"y45 /D");
    tt[i]->Branch("y56",&vars[44],"y56 /D");
    tt[i]->Branch("wz1npfo",&vars[49],"wz1npfo /D");
    tt[i]->Branch("wz2npfo",&vars[50],"wz2npfo /D");
    tt[i]->Branch("pt",&vars[51],"pt /D");
    tt[i]->Branch("hcos",&vars[52],"hcos /D");

    tt[i]->Branch("weight",&vars[99],"weight /D");
  } //all?

  cout << "rootfiles done" << endl;

  //define TMVA
  // create the Reader object
  TMVA::Reader *reader;
  reader= new TMVA::Reader( "!Color:!Silent" );

  float fvars[100];
  //Float_t var1, var2, var3, var4;
  reader->AddVariable( "mh", &fvars[0] );
  reader->AddVariable( "cosjj", &fvars[21] );
  reader->AddVariable( "pt", &fvars[51] );
  reader->AddVariable( "hcos", &fvars[52] );
  reader->AddVariable( "j1e", &fvars[4] );
  reader->AddVariable( "j2e", &fvars[14] );
  reader->AddVariable( "j1cos", &fvars[5] );
  reader->AddVariable( "j2cos", &fvars[15] );
  //reader->AddVariable( "npfo", &fvars[24] );
  //reader->AddVariable( "eviscorr", &fvars[17] );
  reader->AddVariable( "missmass", &fvars[30] );
  // reader->AddVariable( "j1btag", &fvars[6] );
  // reader->AddVariable( "j1ctag", &fvars[7] );
  // reader->AddVariable( "j2btag", &fvars[13] );
  // reader->AddVariable( "j2ctag", &fvars[14] );
  reader->AddVariable( "thrust0", &fvars[34] );
  reader->AddVariable( "thrust1", &fvars[35] );
  //reader->AddVariable( "thrust2", &fvars[29] );
  //reader->AddVariable( "wzmass3", &fvars[32] );
  //reader->AddVariable( "y23", &fvars[41] );
  //reader->AddVariable( "y34", &fvars[42] );
  //reader->AddVariable( "wzmass1", &fvars[37] );
  //reader->AddVariable( "wzmass2", &fvars[38] );
  //reader->AddVariable( "wzmass4", &fvars[33] );
  reader->AddVariable( "y23 * npfo", &fvars[45] );
  reader->AddVariable( "y34 * npfo", &fvars[46] );
  //reader->AddVariable( "j1m / j1nchg", &fvars[47] );
  //reader->AddVariable( "j2m / j2nchg", &fvars[48] );
  
  string dir    = "TMVA/tmva/dataset/weights/";
  string prefix = "TMVAClassification_BDTG_m08p03_";
  string weightfile = "";
  weightfile = dir + prefix + "1.weights.xml";
  reader->BookMVA("BDTG_m08p03_1", weightfile.c_str());
  weightfile = dir + prefix + "2.weights.xml";
  reader->BookMVA("BDTG_m08p03_2", weightfile.c_str());
  weightfile = dir + prefix + "3.weights.xml";
  reader->BookMVA("BDTG_m08p03_3", weightfile.c_str());

  //get parameters here
  ifstream fin("../numofevents/weightlists_l.txt");
  string procdir, fileprefix;
  double cs=0.0, wgt=0.0;
  int filenum=0, evtnum=0;
  vector<string> proclists, filelists;
  vector<double> cslists, weightlists;
  vector<int> filenumlists, numofeventlists;
  double totevtnum=0.0;
  while(!fin.eof()){
    fin >> procdir >> fileprefix >> cs >> wgt >> filenum >> evtnum;

    proclists.push_back(procdir);
    filelists.push_back(fileprefix);
    cslists.push_back(cs);
    weightlists.push_back(wgt);
    filenumlists.push_back(filenum);
    numofeventlists.push_back(evtnum);
  }
  //done

  cout << "read file done" << endl;

  //scan them
  //cut table
  ofstream fout("./cuttable_-0.8+0.3test3.txt");
  fout << "     nnh_bb   nnh_cc   nnh_gg   2f      4f      5f      6f      aa" << endl;
  //for cut table
  const int category = 3;
  double narray[cats][10]={};  //proctype:cuts
  int catnum0[cats]={},catnum[cats]={};
  float mvasn[category][2][51][51][51]={};
  string dirprefix = "../rootfile_l_retrain111/";
  for(unsigned int i=0;i<proclists.size();i++){
    //check whether files are there
    string dir;
    
    // if(!(proclists[i].find("nnh_cc")!=std::string::npos
    // 	 // 	 proclists[i].find("nnh_cc")!=std::string::npos
    // 	 // 	 || proclists[i].find("nnh_gg")!=std::string::npos
    // 	 // 	 //|| proclists[i].find("2f")!=std::string::npos
    // 	 // 	 //|| proclists[i].find("4f")!=std::string::npos
    // 	 )){
    //   dirprefix = "../rootfile_l_retrain111/";
    //   dir = dirprefix + proclists[i];
    // }else{
    //   dirprefix = "../rootfile_l/";
    //   dir = dirprefix + proclists[i];
    // }

    struct stat statBuf;
    if (stat(dir.c_str(), &statBuf) != 0){    
      cout << proclists[i] << " " << "no directory there" << endl;
      continue;
    }    

    if(filenumlists[i]!=0){
      cout << proclists[i] << " starts" << endl;
      TChain *jetchain=new TChain("jetTree");
      TChain *lepchain=new TChain("lepTree");
      TChain *metchain=new TChain("metTree");
      TChain *jettrkchain=new TChain("jettrkTree");
       for(int j=0;j<filenumlists[i]+1;j++){
	string filename = dir + "/" + filelists[i] + "_" + itos(j) + ".root";
	jetchain->Add(filename.c_str());
	lepchain->Add(filename.c_str());
	metchain->Add(filename.c_str());
	jettrkchain->Add(filename.c_str());
      }

      //jetchain
      jetchain->SetBranchAddress("njets", &data.njets);
      jetchain->SetBranchAddress("jetid", &data.jetid[0]);
      jetchain->SetBranchAddress("jete", &data.jete[0]);
      jetchain->SetBranchAddress("jetet", &data.jetet[0]);
      jetchain->SetBranchAddress("jetpx", &data.jetpx[0]);
      jetchain->SetBranchAddress("jetpy", &data.jetpy[0]);
      jetchain->SetBranchAddress("jetpz", &data.jetpz[0]);
      jetchain->SetBranchAddress("btag", &data.btag[0]);
      jetchain->SetBranchAddress("ctag", &data.ctag[0]);
      //jetchain->SetBranchAddress("evis", &data.evis);
      jetchain->SetBranchAddress("ntrack", &data.ntrack[0]);
      jetchain->SetBranchAddress("y23", &data.y23[0]);
      jetchain->SetBranchAddress("y34", &data.y34[0]);
      jetchain->SetBranchAddress("y45", &data.y45[0]);
      jetchain->SetBranchAddress("y56", &data.y56[0]);

      //lepchain
      lepchain->SetBranchAddress("nlep", &data.nlep);

      //metchain
      metchain->SetBranchAddress("metpx", &data.metpx);
      metchain->SetBranchAddress("metpy", &data.metpy);
      metchain->SetBranchAddress("metpz", &data.metpz);
      metchain->SetBranchAddress("thrust0", &data.thrust0);
      metchain->SetBranchAddress("thrust1", &data.thrust1);
      metchain->SetBranchAddress("thrust2", &data.thrust2);
      metchain->SetBranchAddress("wzmass1", &data.wzmass1);
      metchain->SetBranchAddress("wzmass2", &data.wzmass2);
      metchain->SetBranchAddress("wzmass3", &data.wzmass3);
      metchain->SetBranchAddress("wzmass4", &data.wzmass4);
      metchain->SetBranchAddress("wz1npfo", &data.wz1npfo);
      metchain->SetBranchAddress("wz2npfo", &data.wz2npfo);
      
      //jettrkchain
      jettrkchain->SetBranchAddress("tr_npart", &data.tr_npart);
      jettrkchain->SetBranchAddress("nu_npart", &data.nu_npart);
      jettrkchain->SetBranchAddress("tr_e", &data.tr_e[0]);
      jettrkchain->SetBranchAddress("tr_px", &data.tr_px[0]);
      jettrkchain->SetBranchAddress("tr_py", &data.tr_py[0]);
      jettrkchain->SetBranchAddress("tr_pz", &data.tr_pz[0]);
      jettrkchain->SetBranchAddress("nu_e", &data.nu_e[0]);
      jettrkchain->SetBranchAddress("nu_px", &data.nu_px[0]);
      jettrkchain->SetBranchAddress("nu_py", &data.nu_py[0]);
      jettrkchain->SetBranchAddress("nu_pz", &data.nu_pz[0]);

      //start selection
      //set correction factor
      double corrfact = 1.0;
      if(proclists[i].find("eL.pR")!=std::string::npos)
	corrfact *= eminus * pplus;
      if(proclists[i].find("eL.pL")!=std::string::npos)
	corrfact *= eminus * pminus;
      if(proclists[i].find("eR.pL")!=std::string::npos)
	corrfact *= eplus * pminus;
      if(proclists[i].find("eR.pR")!=std::string::npos)
	corrfact *= eminus * pminus;
      if(proclists[i].find("eL.pW")!=std::string::npos)
	corrfact *= eminus * wplus;
      if(proclists[i].find("eL.pB")!=std::string::npos)
	corrfact *= eminus * bplus;
      if(proclists[i].find("eR.pW")!=std::string::npos)
	corrfact *= eplus * wplus;
      if(proclists[i].find("eR.pB")!=std::string::npos)
	corrfact *= eplus * bplus;
      if(proclists[i].find("eW.pL")!=std::string::npos)
	corrfact *= wplus * pminus;
      if(proclists[i].find("eW.pR")!=std::string::npos)
	corrfact *= wplus * pplus;
      if(proclists[i].find("eB.pL")!=std::string::npos)
	corrfact *= bplus * pminus;
      if(proclists[i].find("eB.pR")!=std::string::npos)
	corrfact *= bplus * pplus;

      //set process type and cat
      int cat = -1;
      string procname;
      if(proclists[i].find("nnh_bb")!=std::string::npos){
	procname = "nnh_bb";
	cat = 0;
      }
      if(proclists[i].find("nnh_cc")!=std::string::npos){
	procname = "nnh_cc";
	cat = 1;
      }
      if(proclists[i].find("nnh_gg")!=std::string::npos){
	procname = "nnh_gg";
	cat = 2;
      }
      if(proclists[i].find("2f")!=std::string::npos){
	procname = "2f";
	cat = 3;
      }
      if(proclists[i].find("4f")!=std::string::npos){
	procname = "4f";
	cat = 4;
      }
      if(proclists[i].find("5f")!=std::string::npos){
	procname = "5f";
	cat = 5;
      }
      if(proclists[i].find("6f")!=std::string::npos){
	procname = "6f";
	cat = 6;
      }
      if(proclists[i].find("aa")!=std::string::npos){
	procname = "gghadrons";
	cat = 7;
      }
      if(proclists[i].find("nnh_ww")!=std::string::npos
	 || proclists[i].find("qqh")!=std::string::npos
	 || proclists[i].find("e1e1h")!=std::string::npos
	 || proclists[i].find("e2e2h")!=std::string::npos
	 || proclists[i].find("e3e3h")!=std::string::npos
	 || proclists[i].find("ffh_mumu")!=std::string::npos
	 ){
	procname = "other Higgs";
	cat = 8;
      }

      if(cat<0) continue;

      totevtnum += numofeventlists[i];
      //evtnum = min(0,numofeventlists[i]);  //for debug
      //if(cat==1)
	evtnum = numofeventlists[i];

      int npos=0;
      for(int j=0;j<evtnum;j++){
	catnum0[cat]++;
	if(j>0){
	  npos += data.njets;
	}
	jetchain->GetEntry(j);
	lepchain->GetEntry(j);
	metchain->GetEntry(j);

	//cut0;
	narray[cat][0] += lumi * corrfact * weightlists[i];

	//cut1
	if(data.nlep!=0) continue;
	narray[cat][1] += lumi * corrfact * weightlists[i];

	//cut2
	if(data.njets!=2) continue;
	narray[cat][2] += lumi * corrfact * weightlists[i];

	//cut3
	if(data.y23[0]!=data.y23[0]) continue;
	narray[cat][3] += lumi * corrfact * weightlists[i];

	//cut4
	if(data.jete[0] + data.jete[1]>300.0) continue;
	narray[cat][4] += lumi * corrfact * weightlists[i];
	
	//cut5
	TLorentzVector j1(data.jetpx[0],
			  data.jetpy[0],
			  data.jetpz[0],
			  data.jete[0]);
	
	TLorentzVector j2(data.jetpx[1],
			  data.jetpy[1],
			  data.jetpz[1],
			  data.jete[1]);
	
	if(!((j1+j2).M()>40.0 && (j1+j2).M()<200.0)) continue;
	narray[cat][5] += lumi * corrfact * weightlists[i];
	
	//cut6
	if(data.ntrack[0] + data.ntrack[1] <5) continue;
	narray[cat][6] += lumi * corrfact * weightlists[i];

	//fill histos and rootfile
	//rootfile
	//jet1
	vars[0] = (j1+j2).M();
	vars[1] = j1.Px();
	vars[2] = j1.Py();
	vars[3] = j1.Pz();
	vars[4] = j1.E();
	vars[5] = j1.CosTheta();
	vars[6] = data.btag[0];
	vars[7] = data.ctag[0];
	//quark - gluon separation
	jettrkchain->GetEntry(npos);
	vars[8] = data.tr_npart;
	vars[9] = calBjet(j1, data);
	vars[10] = j1.M();
	//jet2
	vars[11] = j2.Px();
	vars[12] = j2.Py();
	vars[13] = j2.Pz();
	vars[14] = j2.E();
	vars[15] = j2.CosTheta();
	vars[16] = data.btag[1];
	vars[17] = data.ctag[1];
	//quark - gluon separation
	jettrkchain->GetEntry(npos+1);
	vars[18] = data.tr_npart;
	vars[19] = calBjet(j2, data);
	vars[20] = j2.M();
	//other
	vars[21] = TMath::Cos(j1.Vect().Angle(j2.Vect()));
	vars[22] = data.evis;  //with beam backgrounds?
	vars[23] = j1.E() + j2.E();  //without beam backgrounds
	vars[24] = data.ntrack[0] + data.ntrack[1];  //without beam backgrounds
	vars[25] = -(j1+j2).Px();  //without
	vars[26] = -(j1+j2).Py();  //without
	vars[27] = -(j1+j2).Pz();  //without
	vars[28] = 500.0 - vars[22];  //with
	vars[29] = 500.0 - vars[23];  //without

	TLorentzVector miss(vars[25],
			    vars[26],
			    vars[27],
			    vars[29]);
	vars[30] = miss.M();
	vars[31] = vars[6]*vars[16]/(vars[6]* vars[16] + (1-vars[6])*(1-vars[16]));
	vars[32] = vars[7]*vars[17]/(vars[7]* vars[17] + (1-vars[7])*(1-vars[17]));
	double j1bc = vars[7] / (vars[6] + vars[7]);
	double j2bc = vars[17] / (vars[16] + vars[17]);
	vars[33] = j1bc*j2bc / (j1bc*j2bc + (1.0 - j1bc) * (1.0 - j2bc));
	vars[34] = data.thrust0;
	vars[35] = data.thrust1;
	vars[36] = data.thrust2;
	vars[37] = data.wzmass1;
	vars[38] = data.wzmass2;
	vars[39] = data.wzmass3;
	vars[40] = data.wzmass4;
	vars[41] = TMath::Log(data.y23[0]);
	vars[42] = TMath::Log(data.y34[0]);
	vars[43] = TMath::Log(data.y45[0]);
	vars[44] = TMath::Log(data.y56[0]);
	vars[45] = TMath::Log(data.y23[0]) * vars[24];
	vars[46] = TMath::Log(data.y34[0]) * vars[24];
	vars[47] = vars[10] / vars[8];
	vars[48] = vars[20] / vars[18];
	vars[49] = data.wz1npfo;
	vars[50] = data.wz2npfo;
	vars[51] = (j1 + j2).Pt();
	vars[52] = (j1 + j2).CosTheta();

	vars[99] = lumi * corrfact * weightlists[i];

	//cut 7(trivial cut, don't count)
	if(!(miss.M()>0.0 && vars[10]<120.0 && vars[20]<120.0)) continue;
	narray[cat][7] += vars[99];

	float mvaoutput_1 = 0.0;   //reader->EvaluateMulticlass("BDTG")[0];
	float mvaoutput_2 = 0.0;   //reader->EvaluateMulticlass("BDTG")[1];
	float mvaoutput_3 = 0.0;   //reader->EvaluateMulticlass("BDTG")[2];
	if(true){   //start MVA
	  for(int k=0;k<100;k++)
	    fvars[k] = (float) vars[k];
	  const float cut_1 = 0.02 * 17 - 1.0;
	  const float cut_2 = 0.02 * 48;
	  const float cut_3 = 0.02 * 39;

	  mvaoutput_1 = reader->EvaluateMVA("BDTG_m08p03_1");
	  mvaoutput_2 = reader->EvaluateMVA("BDTG_m08p03_2");
	  mvaoutput_3 = reader->EvaluateMVA("BDTG_m08p03_3");

	  for(int k=0;k<50;k++){
	    for(int l=0;l<50;l++){
	      for(int m=0;m<50;m++){
		if(mvaoutput_1>k*0.02-1.0 && mvaoutput_2>l*0.02 && mvaoutput_3>m*0.02){
		  if(cat<3)
		    mvasn[0][0][k][l][m] += vars[99];
		  else //if(cat>=3 && cat<8)
		    mvasn[0][1][k][l][m] += vars[99];
		}
	      }
	    }
	  }
	  
	  mva0[cat]->Fill(mvaoutput_1, vars[99]);
	  mva1[cat]->Fill(mvaoutput_2, vars[99]);
	  mva2[cat]->Fill(mvaoutput_3, vars[99]);

	  if(mvaoutput_1 > cut_1 && mvaoutput_2 > cut_2 && mvaoutput_3 > cut_3)	  
	    narray[cat][8] += vars[99];
	  else
	    continue;
	}

	catnum[cat]++;

	tt[cat]->Fill();
	if(cat>=3 && cat<8)  //fill backgrounds
	  tt[9]->Fill();
	if(cat<3)  //fill signal
	  tt[10]->Fill();
	if(cat==8 && 
	   (proclists[i].find("qqh")!=std::string::npos
	   || proclists[i].find("e1e1h")!=std::string::npos
	   || proclists[i].find("e2e2h")!=std::string::npos
	   || proclists[i].find("e3e3h")!=std::string::npos
	    || proclists[i].find("ffh_mumu")!=std::string::npos)
	   )  //other ZH
	  tt[11]->Fill();
	if(cat==8 &&
	   proclists[i].find("nnh_ww")!=std::string::npos)  //other nnh
	  tt[12]->Fill();

	//histograms
	mh[cat]->Fill(vars[0], vars[99]);
	j1e[cat]->Fill(vars[4], vars[99]);
	j2e[cat]->Fill(vars[14], vars[99]);
	j1cos[cat]->Fill(vars[5], vars[99]);
	j2cos[cat]->Fill(vars[15], vars[99]);
	j1m[cat]->Fill(vars[10], vars[99]);
	j2m[cat]->Fill(vars[20], vars[99]);
	j1Bjet[cat]->Fill(vars[9], vars[99]);
	j2Bjet[cat]->Fill(vars[19], vars[99]);
	j1nchg[cat]->Fill(vars[8], vars[99]);
	j2nchg[cat]->Fill(vars[18], vars[99]);
	cosjj[cat]->Fill(vars[21], vars[99]);
	//flavortag[cat]->Fill(vars[6], vars[7], vars[99]);
	//flavortag[cat]->Fill(vars[16], vars[17], vars[99]);
	flavortag[cat]->Fill(vars[31], vars[32], vars[99]);
	//if(j%2==0)
	  temp[cat]->Fill(vars[31], vars[32], vars[33], vars[99]);
	npfo[cat]->Fill(vars[24], vars[99]);
	evis[cat]->Fill(vars[23], vars[99]);
	missmass[cat]->Fill(vars[30], vars[99]);
	thrust0[cat]->Fill(vars[34], vars[99]);
	thrust1[cat]->Fill(vars[35], vars[99]);
	thrust2[cat]->Fill(vars[36], vars[99]);
	y23[cat]->Fill(vars[41], vars[99]);
	y34[cat]->Fill(vars[42], vars[99]);
	y45[cat]->Fill(vars[43], vars[99]);
	y56[cat]->Fill(vars[44], vars[99]);
	wzmass1[cat]->Fill(vars[37], vars[99]);
	wzmass2[cat]->Fill(vars[38], vars[99]);
	pt[cat]->Fill(vars[51], vars[99]);

	//mva0[cat]->Fill(mvaoutput_1, vars[99]);
	//mva1[cat]->Fill(mvaoutput_2, vars[99]);
	//mva2[cat]->Fill(mvaoutput_3, vars[99]);
      }

      delete jetchain;
      delete lepchain;
      delete metchain;
    }
  }

  for(int i=0;i<cats+4;i++)
    tt[i]->Write();
  tf->Close();

  TPaveText *titletxt=new TPaveText(0.10,0.905,0.40,1.003,"NDC");
  titletxt->SetFillColor(0);
  titletxt->SetBorderSize(0);
  //titletxt->AddText("Higgs Self Coupling Analysis");
  titletxt->AddText("Higgs Branching Ratio Analysis");
  titletxt->AddText("e^{+}e^{-}#rightarrow#nu#nuH");

  TPaveLabel *lumtxt=new TPaveLabel(0.40,0.91,1.0,0.975,"ILD Preliminary","NDC");  //   #intL=9.4fb^{-1}","NDC");
  lumtxt->SetFillColor(0);
  lumtxt->SetBorderSize(0);

  //make legend - need to change here!
  TLegend *ilcleg=new TLegend(0.75,0.58,0.97,0.91,"","brNDC");
  ilcleg->SetFillColor(0);
  ilcleg->SetLineWidth(0.5);

  string label, opt;
  opt="F";
  label="nnh_bb";
  ilcleg->AddEntry(mh[0],label.c_str(),opt.c_str());
  label="nnh_cc";
  ilcleg->AddEntry(mh[1],label.c_str(),opt.c_str());
  label="nnh_gg";
  ilcleg->AddEntry(mh[2],label.c_str(),opt.c_str());
  label="2f";
  ilcleg->AddEntry(mh[3],label.c_str(),opt.c_str());
  label="4f";
  ilcleg->AddEntry(mh[4],label.c_str(),opt.c_str());
  label="5f";
  ilcleg->AddEntry(mh[5],label.c_str(),opt.c_str());
  label="6f";
  ilcleg->AddEntry(mh[6],label.c_str(),opt.c_str());
  label="#gamma#gamma";
  ilcleg->AddEntry(mh[7],label.c_str(),opt.c_str());
  label="other Higgs";
  ilcleg->AddEntry(mh[8],label.c_str(),opt.c_str());

  //make plots
  TCanvas *c1 = new TCanvas("c1","Dynamic Filling Example",300,30,1600,1200);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);
  c1->Divide(2,2);
 
  for(Int_t i=0;i<cats;i++){
    c1->cd(1);
    if(i==0)
      mh[i]->DrawNormalized("hist");
    else
      mh[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c1->cd(2);
    if(i==0)
      cosjj[i]->DrawNormalized("hist");
    else
      cosjj[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c1->cd(3);
    if(i==0)
      pt[i]->DrawNormalized("hist");
    else
      pt[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
  }

  c1->cd(0);
  c1->SaveAs("plots/vars_2jets_preselection_m08p03test3.eps");
  c1->SaveAs("plots/vars_2jets_preselection_m08p03test3.pdf");

  TCanvas *c2 = new TCanvas("c2","Dynamic Filling Example",300,30,2400,1200);
  c2->GetFrame()->SetBorderSize(6);
  c2->GetFrame()->SetBorderMode(-1);
  c2->Divide(3,2);

  for(Int_t i=0;i<cats;i++){
   c2->cd(1);
    if(i==0)
      j1e[i]->DrawNormalized("hist");
    else
      j1e[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c2->cd(2);
    if(i==0)
      j1cos[i]->DrawNormalized("hist");
    else
      j1cos[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c2->cd(3);
    if(i==0)
      j1m[i]->DrawNormalized("hist");
    else
      j1m[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c2->cd(4);
    if(i==0)
      j1Bjet[i]->DrawNormalized("hist");
    else
      j1Bjet[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c2->cd(5);
    if(i==0)
      j1nchg[i]->DrawNormalized("hist");
    else
      j1nchg[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
  }

  c2->cd(0);
  c2->SaveAs("plots/vars_jet1_preselection_m08p03test3.eps");
  c2->SaveAs("plots/vars_jet1_preselection_m08p03test3.pdf");

  TCanvas *c3 = new TCanvas("c3","Dynamic Filling Example",300,30,2400,1200);
  c3->GetFrame()->SetBorderSize(6);
  c3->GetFrame()->SetBorderMode(-1);
  c3->Divide(3,2);

  for(Int_t i=0;i<cats;i++){
   c3->cd(1);
    if(i==0)
      j2e[i]->DrawNormalized("hist");
    else
      j2e[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c3->cd(2);
    if(i==0)
      j2cos[i]->DrawNormalized("hist");
    else
      j2cos[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c3->cd(3);
    if(i==0)
      j2m[i]->DrawNormalized("hist");
    else
      j2m[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c3->cd(4);
    if(i==0)
      j2Bjet[i]->DrawNormalized("hist");
    else
      j2Bjet[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c3->cd(5);
    if(i==0)
      j2nchg[i]->DrawNormalized("hist");
    else
      j2nchg[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
  }

  c3->cd(0);
  c3->SaveAs("plots/vars_jet2_preselection_m08p03test3.eps");
  c3->SaveAs("plots/vars_jet2_preselection_m08p03test3.pdf");

  TCanvas *c4 = new TCanvas("c4","Dynamic Filling Example",300,30,1600,1200);
  c4->GetFrame()->SetBorderSize(6);
  c4->GetFrame()->SetBorderMode(-1);
  c4->Divide(2,2);

  for(Int_t i=0;i<cats;i++){
   c4->cd(1);
    if(i==0)
      npfo[i]->DrawNormalized("hist");
    else
      npfo[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c4->cd(2);
    if(i==0)
      evis[i]->DrawNormalized("hist");
    else
      evis[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c4->cd(3);
    if(i==0)
      missmass[i]->DrawNormalized("hist");
    else
      missmass[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
  }

  c4->cd(0);
  c4->SaveAs("plots/vars_evtvars_preselection_m08p03test3.eps");
  c4->SaveAs("plots/vars_evtvars_preselection_m08p03test3.pdf");


  TCanvas *c5 = new TCanvas("c5","Dynamic Filling Example",300,30,1600,1200);
  c5->GetFrame()->SetBorderSize(6);
  c5->GetFrame()->SetBorderMode(-1);
  c5->Divide(2,2);

  for(Int_t i=0;i<cats;i++){
   c5->cd(1);
    if(i==0)
      thrust0[i]->DrawNormalized("hist");
    else
      thrust0[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c5->cd(2);
    if(i==0)
      thrust1[i]->DrawNormalized("hist");
    else
      thrust1[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c5->cd(3);
    if(i==0)
      thrust2[i]->DrawNormalized("hist");
    else
      thrust2[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
  }

  c5->cd(0);
  c5->SaveAs("plots/vars_thrust_preselection_m08p03test3.eps");
  c5->SaveAs("plots/vars_thrust_preselection_m08p03test3.pdf");

  TCanvas *c6 = new TCanvas("c6","Dynamic Filling Example",300,30,1600,1200);
  c6->GetFrame()->SetBorderSize(6);
  c6->GetFrame()->SetBorderMode(-1);
  c6->Divide(2,2);

  for(Int_t i=0;i<cats;i++){
   c6->cd(1);
    if(i==0)
      y23[i]->DrawNormalized("hist");
    else
      y23[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c6->cd(2);
    if(i==0)
      y34[i]->DrawNormalized("hist");
    else
      y34[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c6->cd(3);
    if(i==0)
      y45[i]->DrawNormalized("hist");
    else
      y45[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c6->cd(4);
    if(i==0)
      y56[i]->DrawNormalized("hist");
    else
      y56[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
  }

  c6->cd(0);
  c6->SaveAs("plots/vars_yvalue_preselection_m08p03test3.eps");
  c6->SaveAs("plots/vars_yvalue_preselection_m08p03test3.pdf");
  
  TCanvas *c7 = new TCanvas("c7","Dynamic Filling Example",300,30,1600,1200);
  c7->GetFrame()->SetBorderSize(6);
  c7->GetFrame()->SetBorderMode(-1);
  c7->Divide(2,2);

  c7->cd(1);
  flavortag[0]->DrawNormalized("lego0");
  lumtxt->Draw();
  titletxt->Draw();
  ilcleg->Draw();
  
  c7->cd(2);
  flavortag[1]->DrawNormalized("lego0");
  lumtxt->Draw();
  titletxt->Draw();
  ilcleg->Draw();
  
  c7->cd(3);
  flavortag[2]->DrawNormalized("lego0");
  lumtxt->Draw();
  titletxt->Draw();
  ilcleg->Draw();

  c7->cd(0);
  c7->SaveAs("plots/vars_flavortag_preselection_m08p03test3.eps");
  c7->SaveAs("plots/vars_flavortag_preselection_m08p03test3.pdf");

  TCanvas *c8 = new TCanvas("c8","Dynamic Filling Example",300,30,1600,600);
  c8->GetFrame()->SetBorderSize(6);
  c8->GetFrame()->SetBorderMode(-1);
  c8->Divide(2,1);

  for(Int_t i=0;i<cats;i++){
   c8->cd(1);
    if(i==0)
      wzmass1[i]->DrawNormalized("hist");
    else
      wzmass1[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();

   c8->cd(2);
    if(i==0)
      wzmass2[i]->DrawNormalized("hist");
    else
      wzmass2[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
  }

  c8->cd(0);
  c8->SaveAs("plots/vars_wzmass_preselection_m08p03test3.eps");
  c8->SaveAs("plots/vars_wzmass_preselection_m08p03test3.pdf");

  TCanvas *c9 = new TCanvas("c9","Dynamic Filling Example",300,30,1600,1200);
  c9->GetFrame()->SetBorderSize(6);
  c9->GetFrame()->SetBorderMode(-1);
  c9->Divide(2,2);
 
  for(Int_t i=0;i<cats;i++){
    c9->cd(1);
    if(i==0)
      mva0[i]->DrawNormalized("hist");
    else
      mva0[i]->DrawNormalized("histsame");

    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
    
    c9->cd(2);
    if(i==0)
      mva1[i]->DrawNormalized("hist");
    else
      mva1[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
    
    c9->cd(3);
    if(i==0)
      mva2[i]->DrawNormalized("hist");
    else
      mva2[i]->DrawNormalized("histsame");
    lumtxt->Draw();
    titletxt->Draw();
    ilcleg->Draw();
  }
  c9->cd(0);
  c9->SaveAs("plots/mvaoutput_m08p03test3.eps");
  c9->SaveAs("plots/mvaoutput_m08p03test3.pdf");

  //make plots for sn
  TCanvas *c10 = new TCanvas("c10","Dynamic Filling Example",300,30,800,600);
  c10->GetFrame()->SetBorderSize(6);
  c10->GetFrame()->SetBorderMode(-1);
  c10->Divide(1,1);

  TH3F *sn = new TH3F("mvasn","",100,-1.0,1.0,100,-1.0,1.0,100,-1.0,1.0);

  double oksn = 0.0;
  int binx =0, biny=0, binz=0;
  for(int i=0;i<50;i++){
    for(int j=0;j<50;j++){
      for(int k=0;k<50;k++){
	double y=0.0;
	if(mvasn[0][0][i][j][k] + mvasn[0][1][i][j][k]!=0.0)
	  y = mvasn[0][0][i][j][k] / sqrt(mvasn[0][0][i][j][k] + mvasn[0][1][i][j][k]);
	
	sn->SetBinContent(i+1, j+1, k+1, y);
	
	if(y>oksn){
	  oksn = y;
	  binx = i;
	  biny = j;
	  binz = k;
	}
      }
    }
  }

  c10->cd(1);
  sn->Draw("colz");
  lumtxt->Draw();
  titletxt->Draw();

  c10->SaveAs("plots/mva_sn_m08p03test3.eps");
  c10->SaveAs("plots/mva_sn_m08p03test3.pdf");

  //save histos
  TFile *tout = new TFile("vars_histograms_m08p03test3.root","RECREATE");
  for(int i=0;i<cats;i++){
    mh[i]->Write();
    cosjj[i]->Write();
    j1e[i]->Write();
    j2e[i]->Write();
    j1cos[i]->Write();
    j2cos[i]->Write();
    j1m[i]->Write();
    j2m[i]->Write();
    j1Bjet[i]->Write();
    j2Bjet[i]->Write();
    j1nchg[i]->Write();
    j2nchg[i]->Write();
    npfo[i]->Write();
    evis[i]->Write();
    missmass[i]->Write();
    thrust0[i]->Write();
    thrust1[i]->Write();
    thrust2[i]->Write();
    y23[i]->Write();
    y34[i]->Write();
    y45[i]->Write();
    y56[i]->Write();
    wzmass1[i]->Write();
    wzmass2[i]->Write();
    flavortag[i]->Write();
    mva0[i]->Write();
    mva1[i]->Write();
    mva2[i]->Write();
  }
  tout->Close();

  TFile *tout2 = new TFile("template_flavortag_m08p03test3.root","RECREATE");
  for(int i=0;i<cats;i++){
    temp[i]->Write();
  }
  tout2->Close();

  //cut table
  fout << "cut0 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][0] << " ";
  fout << endl;
  fout << "cut 1 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][1] << " ";
  fout << endl;
  fout << "cut 2 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][2] << " ";
  fout << endl;
  fout << "cut 3 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][3] << " ";
  fout << endl;
  fout << "cut 4 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][4] << " ";
  fout << endl;
  fout << "cut 5 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][5] << " ";
  fout << endl;
  fout << "cut 6 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][6] << " ";
  fout << endl;
  fout << "cut 7 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][7] << " ";
  fout << endl;
  fout << "cut 8 ";
  for(int i=0;i<cats;i++)
    fout << narray[i][8] << " ";
  fout << endl;

  fout << "cut position: " << binx << " " << biny << " " << binz << " " << oksn << endl;

  fout.close();

  cout << "totevtnum: " << totevtnum << endl;
  cout << "catnum0: " << catnum0[0] << " " 
       << catnum0[1] << " " 
       << catnum0[2] << " " 
       << catnum0[3] << " " 
       << catnum0[4] << " " 
       << catnum0[5] << " " 
       << catnum0[6] << " " 
       << catnum0[7] << " " 
       << catnum0[8] << " " 
       << endl;
  cout << "catnum: " << catnum[0] << " " 
       << catnum[1] << " " 
       << catnum[2] << " " 
       << catnum[3] << " " 
       << catnum[4] << " " 
       << catnum[5] << " " 
       << catnum[6] << " " 
       << catnum[7] << " " 
       << catnum[8] << " " 
       << endl;

  return 0;
}

double calBjet(TLorentzVector jj, jetdata data){
  //charged
  double denom = 0.0, nume =0.0;
  for(int i=0;i<data.tr_npart;i++){
    TLorentzVector tr(data.tr_px[i],
		      data.tr_py[i],
		      data.tr_pz[i],
		      data.tr_e[i]);

    TVector3 kj = jj.Vect().Unit().Cross(tr.Vect());

    denom += tr.P();
    nume += kj.Mag();
  }

  //neutral
  for(int i=0;i<data.nu_npart;i++){
    TLorentzVector nu(data.nu_px[i],
		      data.nu_py[i],
		      data.nu_pz[i],
		      data.nu_e[i]);

    TVector3 kj = jj.Vect().Unit().Cross(nu.Vect());

    denom += nu.P();
    nume += kj.Mag();
  }

  return nume / denom;
}
