#ifndef _myUtil_new_hh_
#define _myUtil_new_hh_

#include <vector>
#include <string>
#include <cstring>
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraph.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#include "TMVA/Tools.h"

struct jetdata
{
  // basic variables for jets
  int processid;
  int evtnum;
  int njets;  // including neutral clusters
  int evtvtx;  // number of vertex in the event
  int npjets; // == 1
  int nnjets; // == 1
  int nmjets; // == 1
  // selection parameters
  int nprong;
  int ntrack[20];
  int nneu[20];
  int pidtrack[20];
  double lpidtrack[20];
  double hitbytrack[20];
  double fracecal[20];
  double evis;
  double jetevis;
  int jetid[20];
  double jete[20];
  double jetet[20];
  double jetpx[20];
  double jetpy[20];
  double jetpz[20];
  double jetphi[20];
  double jeteta[20];
  int jetcharge[20];
  double vtxe[20];
  double vtxpx[20];
  double vtxpy[20];
  double vtxpz[20];
  double vtxx[20];
  double vtxy[20];
  double vtxz[20];
  
  //for b jet correction
  double jetmt[20];
  double trkptsum[20];
  double trkpsum[20];
  double trkmaxpt[20];
  double trkmaxp[20];
  double trkmaxptrel[20];
  double trkptrelsum[20];
  double jetneuEnergy[20];
  double neumaxpt[20];
  double neumaxp[20];
  double neumaxptrel[20];
  double neuptrelsum[20];
  double nearcos[20];
  int nearjetid[20];
  
  //btagging
  double btag[20];
  double ctag[20];
  int category[20];
  double vtxprob1[20];
  double d0bprob[20];
  double d0cprob[20];
  double d0qprob[20];
  double jprobr[20];
  double jprobr5sigma[20];
  double jprobz[20];
  double jprobz5sigma[20];
  int nelectron[20];
  int nmuon[20];
  int ntrk[20];
  int ntrkwithoutv0[20];
  int nvtx[20];
  int nvtxall[20];
  double sphericity[20];
  double trk1d0sig[20];
  double trk1pt[20];
  double trk1pt_jete[20];
  double trk1z0sig[20];
  double trk2d0sig[20];
  double trk2pt[20];
  double trk2pt_jete[20];
  double trk2z0sig[20];
  double trkmass[20];
  double vtxdirang1[20];
  double vtxdirang12[20];
  double vtxdirang12_jete[20];
  double vtxdirang1_jete[20];
  double vtxdirang2[20];
  double vtxdirang2_jete[20];
  double vtxlen1[20];
  double vtxlen1_jete[20];
  double vtxlen2[20];
  double vtxlen2_jete[20];
  double vtxlen12[20];
  double vtxlen12_jete[20];
  double vtxlen12all[20];
  double vtxlen12all_jete[20];
  double vtxmass[20];
  double vtxmass1[20];
  double vtxmass2[20];
  double vtxmassall[20];
  double vtxmasspc[20];
  double vtxmom[20];
  double vtxmom1[20];
  double vtxmom1_jete[20];
  double vtxmom2[20];
  double vtxmom2_jete[20];
  double vtxmom_jete[20];
  double vtxprob[20];
  double vtxsig1[20];
  double vtxsig12[20];
  double vtxsig12_jete[20];
  double vtxsig1_jete[20];
  double vtxsig2[20];
  double vtxsig2_jete[20];
  double z0bprob[20];
  double z0cprob[20];
  double z0qprob[20];
    double CorrVtxMass1[20];
    double CorrVtxMass2[20];
    double CorrVtxMassAll[20];
    double CorrVtxMomentum1[20];
    double CorrVtxMomentum2[20];
    double CorrVtxMomentumAll[20];
    double Pi0Momentum1[20];
    double Pi0Momentum2[20];
    double Pi0MomentumAll[20];
    int NPi0s1[20];
    int NPi0s2[20];

  double bness0[20];
  double bness1[20];
  double bness2[20];
  double bnessmass[20];
  int nbness[20];

    //for my flavor tag study
    double Vtx1Px[20];
    double Vtx1Py[20];
    double Vtx1Pz[20];
    double Vtx1E[20];
    double Vtx2Px[20];
    double Vtx2Py[20];
    double Vtx2Pz[20];
    double Vtx2E[20];
    double Vtx1x[20];
    double Vtx1y[20];
    double Vtx1z[20];
    double Vtx2x[20];
    double Vtx2y[20];
    double Vtx2z[20];
    int Vtx1e_PID[20];
    int Vtx1mu_PID[20];
    int Vtx1pi_PID[20];
    int Vtx1K_PID[20];
    int Vtx1p_PID[20];
    int Vtx1echarge_PID[20];
    int Vtx1mucharge_PID[20];
    int Vtx1picharge_PID[20];
    int Vtx1Kcharge_PID[20];
    int Vtx1pcharge_PID[20];
    int Vtx2e_PID[20];
    int Vtx2mu_PID[20];
    int Vtx2pi_PID[20];
    int Vtx2K_PID[20];
    int Vtx2p_PID[20];
    int Vtx2echarge_PID[20];
    int Vtx2mucharge_PID[20];
    int Vtx2picharge_PID[20];
    int Vtx2Kcharge_PID[20];
    int Vtx2pcharge_PID[20];
    double PriVtxx;
    double PriVtxy;
    double PriVtxz;
  
  //truth information
  double y01[20];
  double y12[20];
  double y23[20];
  double y34[20];
  double y45[20];
  double y56[20];
  double y67[20];
  double y78[20];
  double y89[20];
  double y910[20];

  
  
  //basic variables for lepton
  int nlep;  // including neutral clusters
  int nplep; // == 1
  int nmlep; // == 1
  double lepevis;
  int leptype[5];
  double lepe[5];
  double lepet[5];
  double leppx[5];
  double leppy[5];
  double leppz[5];
  double lepphi[5];
  double lepeta[5];
  int lepcharge[5];
  double ecal[5];
  double hcal[5];
  double mucal[5];
  double coneEnergy[5];
  double ep[5];
  double ehad[5];
  double deltax[5];
  double deltay[5];
  double deltaz[5];
  double d0[5];
  double z0[5];
  double dedx[5];
  double dedxhit[5];
  double hitdedx[300];
  double chisq[5];
  int Ndf[5];
  double jetmindr[5];
  int minjetid[5];
  double minjete[5];
  double chi2[5];
  double max_Ed[5];
  double showerMax[5];
  double absorptionLength[5];
  double showerMax_photon[5];
  double showerMax_ratio[5];
  double Rm[5];
  double shift[5];
  double smax[5];
  double xl20[5];
  double colzmass;

  // positive / negative track
  int charge;
  
  //basic variables for missing Et
  double metpx;
  double metpy;
  double metpz;
  double met;
  double metphi;
  double thrust0;
  double thrust1;
  double thrust2;
  double wzmass1;
  double wzmass2;
  double wzmass3;
  double wzmass4;
  double wz1npfo;
  double wz2npfo;
  
  //cal. W mass
  double testw;
  
  // MC truth(for jet)
  int mc_eventtype;
  int mc_nchild;
  int mc_nb;
  int mc_nq;
  int mc_nz;
  int mc_jetid[20]; // first 5 decay daughters
  int mc_pdg1[20]; // first 5 decay daughters
  int mc_pdgid1[20]; // first 5 decay daughters
  int mc_parentpdg1[20]; // first 5 decay daughters
  double mc_dr1[20]; // first 5 decay daughters
  double mc_drb1[20]; // first 5 decay daughters
  double mc_e1[20];
  double mc_et1[20];
  double mc_px1[20];
  double mc_py1[20];
  double mc_pz1[20];
  double mc_phi1[20];
  double mc_eta1[20];
  int mc_B1[20];
  int mc_B1charge[20];
  int mc_B1ntot[20];
  int mc_B1nok[20];
  int mc_B1ngood[20];
  int mc_B1nbad[20];
  int mc_B1nplus[20];
  int mc_B1nminus[20];
  int mc_B1nzero[20];
  int mc_B1ne[20];
  int mc_B1nmu[20];
  int mc_B1npi[20];
  int mc_B1nk[20];
  int mc_B1np[20];
  double mc_B1badpx[20][20];
  double mc_B1badpy[20][20];
  double mc_B1badpz[20][20];
  double mc_B1bade[20][20];
  double mc_B1badvtxx[20][20];
  double mc_B1badvtxy[20][20];
  double mc_B1badvtxz[20][20];
  int mc_B1badlink[20][20];
  double mc_B1badchi2[20][20];
  double mc_B1baddist[20][20];
  double mc_B1badldist[20][20];
  double mc_B1badpid[20][20];

  double mc_B1goodpx[20][30];
  double mc_B1goodpy[20][30];
  double mc_B1goodpz[20][30];
  double mc_B1goode[20][30];
  double mc_B1goodvtxx[20][30];
  double mc_B1goodvtxy[20][30];
  double mc_B1goodvtxz[20][30];
  int mc_B1goodstat[20][30];
  double mc_B1goodchi2[20][30];
  double mc_B1gooddist[20][30];
  double mc_B1goodldist[20][30];
  double mc_B1goodpid[20][30];

  int mc_pdg2[20]; // first 5 decay daughters
  int mc_pdgid2[20]; // first 5 decay daughters
  int mc_parentpdg2[20]; // first 5 decay daughters
  double mc_dr2[20]; // first 5 decay daughters
  double mc_drb2[20]; // first 5 decay daughters
  double mc_e2[20];
  double mc_et2[20];
  double mc_px2[20];
  double mc_py2[20];
  double mc_pz2[20];
  double mc_phi2[20];
  double mc_eta2[20];
  int mc_B2[20];
  int mc_B2charge[20];
  int mc_B2ntot[20];
  int mc_B2nok[20];
  int mc_B2ngood[20];
  int mc_B2nbad[20];
  int mc_B2nplus[20];
  int mc_B2nminus[20];
  int mc_B2nzero[20];
  int mc_B2ne[20];
  int mc_B2nmu[20];
  int mc_B2npi[20];
  int mc_B2nk[20];
  int mc_B2np[20];
  double mc_B2badpx[20][20];
  double mc_B2badpy[20][20];
  double mc_B2badpz[20][20];
  double mc_B2bade[20][20];
  double mc_B2badvtxx[20][20];
  double mc_B2badvtxy[20][20];
  double mc_B2badvtxz[20][20];
  int mc_B2badlink[20][20];
  double mc_B2badchi2[20][20];
  double mc_B2baddist[20][20];
  double mc_B2badldist[20][20];
  double mc_B2badpid[20][20];

  double mc_B2goodpx[20][30];
  double mc_B2goodpy[20][30];
  double mc_B2goodpz[20][30];
  double mc_B2goode[20][30];
  double mc_B2goodvtxx[20][30];
  double mc_B2goodvtxy[20][30];
  double mc_B2goodvtxz[20][30];
  int mc_B2goodstat[20][30];
  double mc_B2goodchi2[20][30];
  double mc_B2gooddist[20][30];
  double mc_B2goodldist[20][30];
  double mc_B2goodpid[20][30];

  int mcz_jetid[20]; // first 5 decay daughters
  int mcz_pdg1[20]; // first 5 decay daughters
  int mcz_pdgid1[20]; // first 5 decay daughters
  int mcz_parentpdg1[20]; // first 5 decay daughters
  double mcz_dr1[20]; // first 5 decay daughters
  double mcz_e1[20];
  double mcz_et1[20];
  double mcz_px1[20];
  double mcz_py1[20];
  double mcz_pz1[20];
  double mcz_phi1[20];
  double mcz_eta1[20];
  int mcz_pdg2[20]; // first 5 decay daughters
  int mcz_pdgid2[20]; // first 5 decay daughters
  int mcz_parentpdg2[20]; // first 5 decay daughters
  double mcz_dr2[20]; // first 5 decay daughters
  double mcz_e2[20];
  double mcz_et2[20];
  double mcz_px2[20];
  double mcz_py2[20];
  double mcz_pz2[20];
  double mcz_phi2[20];
  double mcz_eta2[20];
  int mcz_B1[20];
  int mcz_B1charge[20];
  
  int mcq_pdg1[20]; // first 5 decay daughters
  int mcq_pdgid1[20]; // first 5 decay daughters
  int mcq_parentpdg1[20]; // first 5 decay daughters
  double mcq_dr1[20]; // first 5 decay daughters
  double mcq_e1[20];
  double mcq_et1[20];
  double mcq_px1[20];
  double mcq_py1[20];
  double mcq_pz1[20];
  double mcq_phi1[20];
  double mcq_eta1[20];
  int mcq_pdg2[20]; // first 5 decay daughters
  int mcq_pdgid2[20]; // first 5 decay daughters
  int mcq_parentpdg2[20]; // first 5 decay daughters
  double mcq_dr2[20]; // first 5 decay daughters
  double mcq_e2[20];
  double mcq_et2[20];
  double mcq_px2[20];
  double mcq_py2[20];
  double mcq_pz2[20];
  double mcq_phi2[20];
  double mcq_eta2[20];

  // MC truth(for lepton)
  int mc2_nchild;
  int eventlep;
  int mc2_pdg1[5]; // first 5 decay daughters
  double mc2_dr1[5]; // first 5 decay daughters
  double mc2_e1[5];
  double mc2_et1[5];
  double mc2_px1[5];
  double mc2_py1[5];
  double mc2_pz1[5];
  double mc2_phi1[5];
  double mc2_eta1[5];
  double mc2_x1[5];
  double mc2_y1[5];
  double mc2_z1[5];
  int mc2_genstat1[5];
  int mc2_pdg2[5]; // first 5 decay daughters
  double mc2_dr2[5]; // first 5 decay daughters
  double mc2_e2[5];
  double mc2_et2[5];
  double mc2_px2[5];
  double mc2_py2[5];
  double mc2_pz2[5];
  double mc2_phi2[5];
  double mc2_eta2[5];
  double mc2_x2[5];
  double mc2_y2[5];
  double mc2_z2[5];
  int mc2_genstat2[5];
  
  //MC truth for(missing Et)
  int mc3_pdg; // first 5 decay daughters
  double mc3_dr; // first 5 decay daughters
  double mc3_e;
  double mc3_et;
  double mc3_px;
  double mc3_py;
  double mc3_pz;
  double mc3_phi;
  double mc3_eta;
  
  // track variables
  int tr_jetid;
  int tr_npart;
  double tr_esum;
  double emcal;
  double hadcal;
  int tr_charge[50];
  int tr_pid[50];
  double tr_fracecal[50];
  double tr_e[50];
  double tr_px[50];
  double tr_py[50];
  double tr_pz[50];
  double tr_x[50];
  double tr_y[50];
  double tr_z[50];
  double tr_d0[50];
  double tr_z0[50];
  double tr_dedx[50];
  double tr_dedxhit[50];
  double tr_dedxerror[50];
  double tr_omega[50];
  double tr_phi[50];
  double tr_tanlambda[50];
  double tr_d0sig[50];
  double tr_z0sig[50];
  int tr_nhit[50];
  int tr_pdg1[50];
  double tr_pdge[50];
  int tr_mother[50];
  double tr_modist[50];
  double tr_mopos[50];
  double tr_dr1[50];
  int tr_qid[50];
  int tr_isoverlay[50];
  
  double cluweight[800];
  double trkweight[800];
  double clupos[800];
  double trkpos[800];
  double cludr[800];
  double trkdr[800];
  double cluphi[800];
  double clucos[800];
    
  //shower profile
  double tr_chi2[50];
  double tr_max_Ed[50];
  double tr_showerMax[50];
  double tr_absorptionLength[50];
  double tr_showerMax_photon[50];
  double tr_showerMax_ratio[50];
  double tr_Rm[50];
  double tr_shift[50];
  double tr_xl20[50];
  double tr_photonEnergy[50];
  double tr_ep[50];
  double tr_ehad[50];
  double tr_mucal[50];
  double tr_deltax[50];
  double tr_deltaz[50];
  
  // neutral parameters
  int nu_npart;
  int nu_npart1;// larger than 1GeV
  int nu_pid[50];
  double nu_fracecal[50];
  double nu_e[50];
  double nu_px[50];
  double nu_py[50];
  double nu_pz[50];
  double nu_x[50];
  double nu_y[50];
  double nu_z[50];
  double nu_esum;
  int nu_pdg1[50];
  double nu_pdge[50];
  double nu_dr1[50];
  int nu_qid[50];
  int nu_isoverlay[50];
   
  //shower profile
  double nu_chi2[50];
  double nu_max_Ed[50];
  double nu_showerMax[50];
  double nu_absorptionLength[50];
  double nu_showerMax_photon[50];
  double nu_showerMax_ratio[50];
  double nu_Rm[50];
  double nu_shift[50];
  double nu_xl20[50];
  double nu_ehad[50];

  //vertex study for b-tagging
  int vtx_nvtx;
  int vtx_id;
  int vtx_charge;
  int vtx_npart;
  double vtx_px;
  double vtx_py;
  double vtx_pz;
  double vtx_e;
  double vtx_length;
  double vtx_x;
  double vtx_y;
  double vtx_z;
  int vtx_nearestPDG;
  int vtx_nearestID;
  double vtx_nearestdist;
  double vtx_nearestdistxy;
  double vtx_nearestdistz;
  double vtx_nearestdistdelphi;
  double vtx_nearestdistdeleta;
  int vtx_nearestPDG2;
  int vtx_nearestID2;
  double vtx_nearestdist2;
  double vtx_nearestdist2xy;
  double vtx_nearestdist2z;
  double vtx_nearestdist2delphi;
  double vtx_nearestdist2deleta;
  double vtx_tr_charge[15];
  int vtx_tr_pid[15];
  double vtx_tr_px[15];
  double vtx_tr_py[15];
  double vtx_tr_pz[15];
  double vtx_tr_x[15];
  double vtx_tr_y[15];
  double vtx_tr_z[15];
  double vtx_tr_e[15];
  double vtx_tr_ep[15];
  double vtx_tr_ehad[15];
  double vtx_tr_mucal[15];
  double vtx_tr_deltax[15];
  double vtx_tr_deltaz[15];
  double vtx_tr_dEdx[15];
  double vtx_tr_dEdxhit[15];
  double vtx_tr_showerchi2[15];
  double vtx_tr_showerMax[15];
  double vtx_tr_showerMax_ratio[15];
  double vtx_tr_absLength[15];
  double vtx_tr_xl20[15];
  double vtx_mc_pdg[15];
  double vtx_mc_pdg_mother[15];
  double vtx_mc_pdg_mother2[15];
  double vtx_mc_dr[15];
  double vtx_mc_neupx;
  double vtx_mc_neupy;
  double vtx_mc_neupz;
  double vtx_mc_neue;
  
  //pi0 study for b-tagging
  int npi0;
  int recoflg[70];
  double g1_px[70];
  double g1_py[70];
  double g1_pz[70];
  double g1_e[70];
  double g1_x[70];
  double g1_y[70];
  double g1_z[70];
  double g2_px[70];
  double g2_py[70];
  double g2_pz[70];
  double g2_e[70];
  double g2_x[70];
  double g2_y[70];
  double g2_z[70];
  
  double g1_ehad[70];
  double g1_showerchi2[70];
  double g1_showerMax[70];
  double g1_showerMax_ratio[70];
  double g1_absLength[70];
  double g1_xl20[70];
  double g2_ehad[70];
  double g2_showerchi2[70];
  double g2_showerMax[70];
  double g2_showerMax_ratio[70];
  double g2_absLength[70];
  double g2_xl20[70];
  double g1_pdg_px[70];
  double g1_pdg_py[70];
  double g1_pdg_pz[70];
  double g1_pdg_e[70];
  double g2_pdg_px[70];
  double g2_pdg_py[70];
  double g2_pdg_pz[70];
  double g2_pdg_e[70];
  double pi0_x[70];
  double pi0_y[70];
  double pi0_z[70];

  //neu study for b-tagging
  int neu_nneu;
  double neu_px[50];
  double neu_py[50];
  double neu_pz[50];
  double neu_e[50];
  double neu_x[50];
  double neu_y[50];
  double neu_z[50];
  double neu_ehad[50];
  double neu_showerchi2[50];
  double neu_showerMax[50];
  double neu_showerMax_ratio[50];
  double neu_absLength[50];
  double neu_xl20[50];
 
  //B study for b-tagging
  int Btype;
  int nD;
  int nJpsi;
  int nLambda;
  int resopdg;
  int nBpip;
  int nBpim;
  int nBpi0;
  int nBKp;
  int nBKm;
  int nBK0;
  int nBe;
  int nBmu;
  int nBtau;
  int nReso;
  int Ddaughter1[10];
  int Ddaughter2[10];
  int Jpsidaughter[10];
  int Lambdadaughter[10];
  int Dprong1;
  int Dprong2;
  double BdecayLength;
  double D1decayLength;
  double D2decayLength;



};

class CombinationSolver : public TObject{
public:
  CombinationSolver(){};
  CombinationSolver(string fstr);
  ~CombinationSolver();
  
  double get_combination(jetdata data,int nbtagnum, int ljetid[], int *wid, int *zid);
  double cal_chisq(vector<TLorentzVector> wvect, vector<TLorentzVector> zvect, TLorentzVector lep, TLorentzVector neu);

private:
  double get_likelihood(double costheta, int costype, int frametype);
  double get_likelihood_WZ(double costheta);
  double get_likelihood_WW(double costheta){return 0;}
  double get_likelihood_jj(double costheta, int jettype);

  TFile *fpdf[2];
  TH1F *pdf[4][2];
  TH1F *pdf_WZ;
  TH1F *pdf_WW;
  TH1F *pdf_jj[2];

  const double mmw=80.385;
  const double wgamma=2.085;
  const double zmeas=91.01917;
  const double hmeas=119.051;
  const double zgammameas=9.13036;
  const double hgammameas=12.1111;

  ClassDef(CombinationSolver, 1);
};

class CombinationSolver_bbbb : public TObject{
public:
  CombinationSolver_bbbb(){};
  CombinationSolver_bbbb(string fstr);
  ~CombinationSolver_bbbb();
  
  double get_combination(jetdata data,int btagnum, int bjetid[], int *hid, int *zid);
  double cal_chisq(vector<TLorentzVector> hvect, vector<TLorentzVector> zvect);

private:
  double get_likelihood(double costheta, int costype);
  double get_likelihood_ZH(double costheta);

  TFile *fpdf;
  TH1F *pdf[4];

  const double mmw=80.385;
  const double wgamma=2.085;
  const double zmeas=87.9638;
  const double hmeas=120.0303;
  double zgammameas;  //=8.72260;
  double hgammameas;  //=10.8044;

  ClassDef(CombinationSolver_bbbb, 1);
};

class CombinationSolver_bbb : public TObject{
public:
  CombinationSolver_bbb(){};
  CombinationSolver_bbb(string fstr);
  ~CombinationSolver_bbb();
  
  double get_combination(jetdata data,int btagnum, int nbtagnum, int bjetid[], int ljetid[], int *hid, int *zid);
  double cal_chisq(vector<TLorentzVector> hvect,
		   vector<TLorentzVector> zvect,
		   vector<TLorentzVector> wvect,
		   TLorentzVector lep, TLorentzVector neu
		   );

private:
  double get_likelihood1(double costheta, int costype);
  double get_likelihood2(double costheta, int costype, int frametype);
  double get_likelihood_ZH(double costheta);
  double get_likelihood_WZ(double costheta);

  TFile *fpdf[3];
  TH1F *pdf[4][3];
  
  const double mmw=80.385;
  const double wgamma=2.085;
  const double zmeas=87.9638;
  const double hmeas=120.0303;
  const double zgammameas=8.72260;
  const double hgammameas=10.8044;
  const double hmeas2=119.051;
  const double hgammameas2=12.1111;

  ClassDef(CombinationSolver_bbb, 1);
};

class CombinationSolver_bbbb_allhad : public TObject{
public:
  CombinationSolver_bbbb_allhad(){};
  CombinationSolver_bbbb_allhad(string fstr);
  ~CombinationSolver_bbbb_allhad();
  
  double get_combination(jetdata data,int btagnum, int nbtagnum, int bjetid[], int ljetid[], int *hid, int *zid, int *w1id, int *w2id);
  double cal_chisq(vector<TLorentzVector> hvect, vector<TLorentzVector> zvect,
		   vector<TLorentzVector> w1vect,vector<TLorentzVector> w2vect);
  
private:
  double get_likelihood1(double costheta, int costype);
  double get_likelihood2(double costheta, int costype);
  double get_likelihood_ZH(double costheta);
  
  TFile *fpdf[2];
  TH1F *pdf[6][3];
  
  const double mmw=80.385;
  const double wgamma=2.085;
  const double zmeas=91.2;   //87.9638;
  const double hmeas=125.00;   //120.0303;
  double zgammameas;  //=8.72260;
  double hgammameas;  //=10.8044;
  const double hmeas2=125.00;   //125.562;
  const double hgammameas2=14.5718;

  ClassDef(CombinationSolver_bbbb_allhad, 1);
};

class CombinationSolver_dilepton : public TObject{
public:
  CombinationSolver_dilepton(){};
  CombinationSolver_dilepton(string fstr);
  ~CombinationSolver_dilepton();
  
  double get_combination(jetdata data,int nbtagnum, int ljetid[], int *w1id, int *w2id);
  double cal_chisq(vector<TLorentzVector> w1vect,vector<TLorentzVector> w2vect);
  
private:
  double get_likelihood2(double costheta, int costype);
  
  TFile *fpdf;
  TH1F *pdf[6];
  
  const double mmw=80.385;
  const double wgamma=2.085;
  const double hmeas=125.562;
  const double hgammameas=14.5718;

  ClassDef(CombinationSolver_dilepton, 1);
};

class mvafile : public TObject{
public:
  mvafile(){};
  mvafile(string fstr, int fieldtype);
  ~mvafile();
  void fill(double *par, int field);
  void write(int fieldtype);
  void close();

private:
  TFile* ftmva;
  TTree* tts[50];
  double var[200];

  ClassDef(mvafile, 1);
};


class jetEnergyScale : public TObject{
public:
  jetEnergyScale();
  ~jetEnergyScale();
  Double_t getScale(Int_t type, Float_t *val, Int_t valtype);

private:
  TMVA::Reader *reader[8]; 
  Float_t var[50];

  ClassDef(jetEnergyScale, 1);
};

class expert_ttbar : public TObject{
public:
  expert_ttbar(){};
  expert_ttbar(int cattype);
  ~expert_ttbar();
  double getMVAvalue(float *par);

private:
  TMVA::Reader *reader[50]; 
  float var[200];
  int type;

  ClassDef(expert_ttbar, 1);
};

class expert_ttbb : public TObject{
public:
  expert_ttbb(){};
  expert_ttbb(int cattype);
  ~expert_ttbb();
  double getMVAvalue(float *par);

private:
  TMVA::Reader *reader[50]; 
  float var[200];
  int type;

  ClassDef(expert_ttbb, 1);
};

class softjetFinder : public TObject{
public:
  softjetFinder();
  ~softjetFinder();
  void get_nearjetPull(int jetid, jetdata data);
  void get_jetPull(int jetid, jetdata data);
  void get_Bjet(int jetid, jetdata data);
  void get_Aa(int jetid, jetdata data);
  void get_twomom(int jetid, jetdata data);
  void get_geomom(int jetid, jetdata data);
  void get_trkNN(int jetid, jetdata data);
  double get_MVAvalue(int jetid, jetdata data); //please do getEntry() before calling this function!
  TMVA::Reader *disc_track;
  TMVA::Reader *disc_track2;
  TMVA::Reader *disc_cluster;
  TMVA::Reader *disc_cluster2;

private:
  //TH1F *tmphist[2];
  float fvar[20];
  float sjvar[15];
  double pullangle[2];
  double pullangle2[2];
  double Bjet;
  double Aa;
  double Tbeta;
  double quartic;
  double eccentricity;
  double pflow;


  //TMVA::Reader *disc_cluster;
  //TMVA::Reader *disc_track;
  TMVA::Reader *disc_softjet;
  TMVA::Reader *disc_softjet2;

  void get_Pullvalue(int jetid1, int jetid2, jetdata data, int type);

  ClassDef(softjetFinder, 1);
};

void setColor(TH1F *h, int type);
double calDelPhi(double phi1, double phi2);
string itos(int i);
string ftos(float i);
double calWeight(int leptype, int proctype, double evtnum, double sigcs);
double calWeight_sig(int leptype, int proctype, double evtnum, double sigcs);
double corrEvtWeight(double evtweight, int proctype);
void jetclustering_2jets(int snjets,int cjets, float *var, double *jets);
TH1F* make_SNplot(int proctype, TH1F *s, TH1F *b);
TH1F* make_Splot(int proctype, TH1F *s);
TH1F* make_SNplotRev(int proctype, TH1F *s, TH1F *b);
TH1F* make_SplotRev(int proctype, TH1F *s);
TH1F* make_SNplotAbs(int proctype, TH1F *s, TH1F *b);
TH1F* make_SplotAbs(int proctype, TH1F *s);

double cal_moment(jetdata data, double Evis, int order);
double cal_sphericity(jetdata data, int type);

#endif
