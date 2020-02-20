#include "TMatrixDEigen.h"
#include "myUtil_new.hh"

ClassImp(CombinationSolver);
ClassImp(CombinationSolver_bbbb);
ClassImp(CombinationSolver_bbb);
ClassImp(CombinationSolver_bbbb_allhad);
ClassImp(CombinationSolver_dilepton);
ClassImp(mvafile);
ClassImp(jetEnergyScale);
ClassImp(expert_ttbar);
ClassImp(expert_ttbb);
ClassImp(softjetFinder);

using namespace std;

//combination solver for bb case
CombinationSolver::CombinationSolver(string fstr){
  string ffstr1= fstr +"1.root";
  string ffstr2= fstr +"2.root";

  fpdf[0]=new TFile(ffstr1.c_str());
  fpdf[1]=new TFile(ffstr2.c_str());

  pdf[0][0]=(TH1F*)fpdf[0]->Get("deltacosthetalepj3;1")->Clone("pdf11");
  pdf[1][0]=(TH1F*)fpdf[0]->Get("deltacosthetalepj4;1")->Clone("pdf21");
  pdf[2][0]=(TH1F*)fpdf[0]->Get("deltacosthetalepj5;1")->Clone("pdf31");
  pdf[3][0]=(TH1F*)fpdf[0]->Get("deltacosthetalepj6;1")->Clone("pdf41");
  pdf[0][1]=(TH1F*)fpdf[1]->Get("deltacosthetalepj3;1")->Clone("pdf12");
  pdf[1][1]=(TH1F*)fpdf[1]->Get("deltacosthetalepj4;1")->Clone("pdf22");
  pdf[2][1]=(TH1F*)fpdf[1]->Get("deltacosthetalepj5;1")->Clone("pdf32");
  pdf[3][1]=(TH1F*)fpdf[1]->Get("deltacosthetalepj6;1")->Clone("pdf42");

  pdf_WZ=(TH1F*)fpdf[1]->Get("deltaRljlj1")->Clone("pdfWZ");

  pdf_jj[0]=(TH1F*)fpdf[0]->Get("deltacosthetajj1")->Clone("pdf_jj_W");
  pdf_jj[1]=(TH1F*)fpdf[0]->Get("deltacosthetajj2")->Clone("pdf_jj_Z");

  //normalize histograms
  double weight=1.0;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<2;j++){
      weight=pdf[i][j]->Integral();
      pdf[i][j]->Scale(1/weight);
    }
  }
    
  weight=pdf_WZ->Integral();
  pdf_WZ->Scale(1/weight);
  
  for(Int_t i=0;i<2;i++){
    weight=pdf_jj[i]->Integral();
    pdf_jj[i]->Scale(1/weight);
  }

  return;
}

CombinationSolver::~CombinationSolver(){
  fpdf[0]->Close();
  fpdf[1]->Close();
}

double CombinationSolver::cal_chisq(vector<TLorentzVector> wvect, vector<TLorentzVector> zvect, TLorentzVector lep, TLorentzVector neu){
  
  //cal boosting vector for W&Z
  //boost in expected Higgs Frame
  TLorentzVector b1(-(wvect[0].Px()+wvect[1].Px()+lep.Px()+neu.Px()),
		    -(wvect[0].Py()+wvect[1].Py()+lep.Py()+neu.Py()),
		    -(wvect[0].Pz()+wvect[1].Pz()+lep.Pz()+neu.Pz()),
		    wvect[0].E()+wvect[1].E()+lep.E()+neu.E());
  
  //boost in expected wrong combination
  TLorentzVector b2(-(zvect[0].Px()+zvect[1].Px()+lep.Px()+neu.Px()),
		    -(zvect[0].Py()+zvect[1].Py()+lep.Py()+neu.Py()),
		    -(zvect[0].Pz()+zvect[1].Pz()+lep.Pz()+neu.Pz()),
		    zvect[0].E()+zvect[1].E()+lep.E()+neu.E());

  //make LorentzVector for lepton and jets
  TLorentzVector wjet[2][2];
  wjet[0][0]=wvect[0];
  wjet[0][1]=wvect[0];
  wjet[1][0]=wvect[1];
  wjet[1][1]=wvect[1];

  TLorentzVector zjet[2][2];
  zjet[0][0]=zvect[0];
  zjet[0][1]=zvect[0];
  zjet[1][0]=zvect[1];
  zjet[1][1]=zvect[1];

  TLorentzVector wlep[2];
  wlep[0]=lep;
  wlep[1]=lep;

  //make LorentzVector for W and Z
  TLorentzVector vw[2];
  
  vw[0].SetPxPyPzE(wvect[0].Px()+wvect[1].Px(),
		   wvect[0].Py()+wvect[1].Py(),
		   wvect[0].Pz()+wvect[1].Pz(),
		   wvect[0].E()+wvect[1].E());
  
  vw[1].SetPxPyPzE(wvect[0].Px()+wvect[1].Px(),
		   wvect[0].Py()+wvect[1].Py(),
		   wvect[0].Pz()+wvect[1].Pz(),
		   wvect[0].E()+wvect[1].E());
  
  TLorentzVector vz[2];
  
  vz[0].SetPxPyPzE(zvect[0].Px()+zvect[1].Px(),
		   zvect[0].Py()+zvect[1].Py(),
		   zvect[0].Pz()+zvect[1].Pz(),
		   zvect[0].E()+zvect[1].E());
  
  vz[1].SetPxPyPzE(zvect[0].Px()+zvect[1].Px(),
		   zvect[0].Py()+zvect[1].Py(),
		   zvect[0].Pz()+zvect[1].Pz(),
		   zvect[0].E()+zvect[1].E());
  
  //make LorentzVector for W decaying into lep & neutrino
  TLorentzVector vwlnu[2];

  vwlnu[0].SetPxPyPzE(lep.Px()+neu.Px(),
		      lep.Py()+neu.Py(),
		      lep.Pz()+neu.Pz(),
		      lep.E()+neu.E());
  
  vwlnu[1].SetPxPyPzE(lep.Px()+neu.Px(),
		      lep.Py()+neu.Py(),
		      lep.Pz()+neu.Pz(),
		      lep.E()+neu.E());
  
  //start to calculate chi square
  //first, calculate Higgs and Z mass
  double mh=b1.M();
  double mz=vz[0].M();

 
  //second, cal W(*) mass for Breit-Wigner
  double mw=vw[0].M();

  //third, cal angle
  //(1) boosting various vectors //0:higgs rest frame 1:W(lnu)+Z rest frame
  wjet[0][0].Boost(b1.BoostVector());
  wjet[0][1].Boost(b2.BoostVector());
  wjet[1][0].Boost(b1.BoostVector());
  wjet[1][1].Boost(b2.BoostVector());
  zjet[0][0].Boost(b1.BoostVector());
  zjet[0][1].Boost(b2.BoostVector());
  zjet[1][0].Boost(b1.BoostVector());
  zjet[1][1].Boost(b2.BoostVector());
  wlep[0].Boost(b1.BoostVector());
  wlep[1].Boost(b2.BoostVector());
  vw[0].Boost(b1.BoostVector());
  vw[1].Boost(b2.BoostVector());
  vz[0].Boost(b1.BoostVector());
  vz[1].Boost(b2.BoostVector());
  vwlnu[0].Boost(b1.BoostVector());
  vwlnu[1].Boost(b2.BoostVector());

  double goodcos[4];
  double badcos[4];

  goodcos[0]=TMath::Cos(wjet[0][0].Angle(wlep[0].Vect()));
  goodcos[1]=TMath::Cos(wjet[1][0].Angle(wlep[0].Vect()));
  goodcos[2]=TMath::Cos(zjet[0][0].Angle(wlep[0].Vect()));
  goodcos[3]=TMath::Cos(zjet[1][0].Angle(wlep[0].Vect()));
  badcos[0]=TMath::Cos(wjet[0][1].Angle(wlep[1].Vect()));
  badcos[1]=TMath::Cos(wjet[1][1].Angle(wlep[1].Vect()));
  badcos[2]=TMath::Cos(zjet[0][1].Angle(wlep[1].Vect()));
  badcos[3]=TMath::Cos(zjet[1][1].Angle(wlep[1].Vect()));
  
  double coswz,cosww;
  coswz=TMath::Cos(vz[0].Angle(vwlnu[0].Vect()));  //cos in Higgs rest frame
  cosww=TMath::Cos(vw[1].Angle(vwlnu[1].Vect()));  //cos in W+Z rest frame

  double lh1,lh2,lh3,lh4;

  //get likelihood value far side: using
  if(goodcos[0]>goodcos[1]){
    lh2=get_likelihood(goodcos[1],1,0);
  }else{
    lh2=get_likelihood(goodcos[0],1,0);
  }

  if(goodcos[2]>goodcos[3]){
    lh4=get_likelihood(goodcos[3],3,0);
  }else{
    lh4=get_likelihood(goodcos[2],3,0);
  }
  
  if(badcos[0]>badcos[1]){
    lh1=get_likelihood(badcos[0],0,1);
  }else{
    lh1=get_likelihood(badcos[1],0,1);
  }
  
  if(badcos[2]>badcos[3]){
    lh3=get_likelihood(badcos[2],2,1);
  }else{
    lh3=get_likelihood(badcos[3],2,1);
  }
  
  //get likelihood for the angle of WW or WZ
  double lh5=get_likelihood_WZ(cosww);
  
  //get likelihood between jets in Higgs rest frame
  double lh6=get_likelihood_jj(TMath::Cos(wjet[0][0].Angle(wjet[1][0].Vect())),0);
  double lh7=get_likelihood_jj(TMath::Cos(zjet[0][0].Angle(zjet[1][0].Vect())),1);

  //now calculate chi square!!
  double chi1=pow(mh-hmeas,2.0)/pow(hgammameas,2.0)+TMath::Log(2*TMath::Pi()*hgammameas*hgammameas);
  double chi2=pow(mz-zmeas,2.0)/pow(zgammameas,2.0)+TMath::Log(2*TMath::Pi()*zgammameas*zgammameas);
  double chi3=2*TMath::Log((1/TMath::Pi())*wgamma*mmw
			   /(pow(mw*mw-mmw*mmw,2.0)+mmw*mmw*wgamma*wgamma));
 
  if(lh1*lh4==0.0) return 1.0e+20;

  //cout << "result1: " << -0.5*chi1 << " " << -0.5*chi2  << " " << -0.5*chi3 << endl;
  //cout << "result2: " << lh1 << " " << lh2  << " " << lh3  << " " << lh4  << " " << lh5 << endl;
  //cout << "result3: " << lh6 << " " << lh7 << endl;
  //-chi3
  //double chiall=chi1+chi2-chi3-2*TMath::Log(lh1)-2*TMath::Log(lh4);
  double chiall=chi1+chi2-chi3-2*TMath::Log(lh1)-2*TMath::Log(lh4)-2*TMath::Log(lh5);  //)-2*TMath::Log(lh7);

  return chiall;
}

//costype: whether the jet is near to the lepton   1:nearW 2:farW 3:near Z 4:far Z 
//frametype: Higgs rest frame or W+Z rest frame    1:Higgs 2: W+Z
double CombinationSolver::get_likelihood(double costheta, int costype, int frametype){
  
  int nbins=pdf[costype][frametype]->GetNbinsX();
  double interval=pdf[costype][frametype]->GetBinWidth(1);



  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=-1.0+i*interval && costheta<-1.0+(i+1)*interval){
      val=pdf[costype][frametype]->GetBinContent(i+1);
      //cout << "val_Wjet: " << costheta << " " << val << endl;
      break;
    }
  }

  return val;
}

double CombinationSolver::get_likelihood_WZ(double costheta){

  int nbins=pdf_WZ->GetNbinsX();
  double interval=pdf_WZ->GetBinWidth(1);
  
  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=-1.0+i*interval && costheta<-1.0+(i+1)*interval){
      val=pdf_WZ->GetBinContent(i+1);
      //cout << "val_WZ: " << costheta << " " << val << endl;
      break;
    }
  }
  
  return val;
}

double CombinationSolver::get_likelihood_jj(double costheta, int jettype){

  int nbins=pdf_jj[jettype]->GetNbinsX();
  double interval=pdf_jj[jettype]->GetBinWidth(1);
  
  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=-1.0+i*interval && costheta<-1.0+(i+1)*interval){
      val=pdf_jj[jettype]->GetBinContent(i+1);
      //cout << "val_jj: " << costheta << " " << val << endl;
      break;
    }
  }
  
  return val;
}

double CombinationSolver::get_combination(jetdata data,int nbtagnum, int ljetid[], int *wid, int *zid){
  double chisq=1.0e20;

  if(nbtagnum<4) return chisq;
  
  TLorentzVector jets[4];
  jets[0].SetPxPyPzE(data.jetpx[ljetid[0]],
		     data.jetpy[ljetid[0]],
		     data.jetpz[ljetid[0]],
		     data.jete[ljetid[0]]);
  jets[1].SetPxPyPzE(data.jetpx[ljetid[1]],
		     data.jetpy[ljetid[1]],
		     data.jetpz[ljetid[1]],
		     data.jete[ljetid[1]]);
  jets[2].SetPxPyPzE(data.jetpx[ljetid[2]],
		     data.jetpy[ljetid[2]],
		     data.jetpz[ljetid[2]],
		     data.jete[ljetid[2]]);
  jets[3].SetPxPyPzE(data.jetpx[ljetid[3]],
		     data.jetpy[ljetid[3]],
		     data.jetpz[ljetid[3]],
		     data.jete[ljetid[3]]);

  TLorentzVector lep(data.leppx[0],
		     data.leppy[0],
		     data.leppz[0],
		     data.lepe[0]);

  TLorentzVector neu(data.metpx,
		     data.metpy,
		     data.metpz,
		     sqrt(pow(data.metpx,2.0)+pow(data.metpy,2.0)+pow(data.metpz,2.0)));

  vector<TLorentzVector> wvect;
  vector<TLorentzVector> zvect;
  for(Int_t i=0;i<nbtagnum;i++){
    for(Int_t j=i+1;j<nbtagnum;j++){
      wvect.push_back(jets[i]);
      wvect.push_back(jets[j]);

      //calculate others
      int k,l;
      switch(i){
      case 0:
	switch(j){
	case 1:
	  k=2;
	  l=3;
	  break;
	case 2:
	  k=1;
	  l=3;
	  break;
	case 3:
	  k=1;
	  l=2;
	  break;
	}
	break;
      case 1:
	switch(j){
	case 2:
	  k=0;
	  l=3;
	  break;
	case 3:
	  k=0;
	  l=2;
	  break;
	}
	break;
      case 2:
	switch(j){
	case 3:
	  k=0;
	  l=1;
	  break;
	}
	break;
      }
 
      //cout << "loop: " << ljetid[i] << " " << ljetid[j]
      //	    << " " << ljetid[k] << " " << ljetid[l] << endl;
     
      zvect.push_back(jets[k]);
      zvect.push_back(jets[l]);
      //get chi square
      double tmpchi=cal_chisq(wvect,zvect,lep,neu);
      if(chisq>tmpchi){
	chisq=tmpchi;
	wid[0]=ljetid[i];
	wid[1]=ljetid[j];
	zid[0]=ljetid[k];
	zid[1]=ljetid[l];
	//cout << "tmpchi: " << tmpchi << endl;
      }
      
      wvect.clear();
      zvect.clear();
    }
  }
  
  return chisq;
}

//combination solver for bbbb case
CombinationSolver_bbbb::CombinationSolver_bbbb(string fstr){
  string ffstr1= fstr +".root";

  fpdf=new TFile(ffstr1.c_str());

  pdf[0]=(TH1F*)fpdf->Get("yZH1;1")->Clone("pdf11");
  pdf[1]=(TH1F*)fpdf->Get("yZH2;1")->Clone("pdf21");
  pdf[2]=(TH1F*)fpdf->Get("yZH3;1")->Clone("pdf31");

  //normalize histograms
  double weight=1.0;
  for(Int_t i=0;i<3;i++){
    weight=pdf[i]->Integral();
    pdf[i]->Scale(1/weight);
  }

  return;
}

CombinationSolver_bbbb::~CombinationSolver_bbbb(){
  fpdf->Close();
}

double CombinationSolver_bbbb::cal_chisq(vector<TLorentzVector> hvect, vector<TLorentzVector> zvect){
  
  //make LorentzVector for lepton and jets
  TLorentzVector hjet[2];
  hjet[0]=hvect[0];
  hjet[1]=hvect[1];

  TLorentzVector zjet[2];
  zjet[0]=zvect[0];
  zjet[1]=zvect[1];

  //make LorentzVector for W and Z
  TLorentzVector vh(hvect[0].Px()+hvect[1].Px(),
		    hvect[0].Py()+hvect[1].Py(),
		    hvect[0].Pz()+hvect[1].Pz(),
		    hvect[0].E()+hvect[1].E());
  
  TLorentzVector vz(zvect[0].Px()+zvect[1].Px(),
		    zvect[0].Py()+zvect[1].Py(),
		    zvect[0].Pz()+zvect[1].Pz(),
		    zvect[0].E()+zvect[1].E());
  
  //start to calculate chi square
  //first, calculate Higgs and Z mass
  double mh=vh.M();
  double mz=vz.M();

  //second, cal angle  
  double coszh;
  coszh=TMath::Cos(vz.Angle(vh.Vect()));  //cos in Higgs rest frame

  double lh1,lh2,lh3;

  //get likelihood value of angle
  double yy=0.0;
  yy=TMath::Min(pow(hjet[0].E(),2.0),pow(hjet[1].E(),2.0))
    *(1-TMath::Cos(hjet[0].Angle(hjet[1].Vect())))/pow(500.0, 2.0);
  lh1=get_likelihood(yy,0);
  yy=TMath::Min(pow(zjet[0].E(),2.0),pow(zjet[1].E(),2.0))
    *(1-TMath::Cos(zjet[0].Angle(zjet[1].Vect())))/pow(500.0, 2.0);
  lh2=get_likelihood(yy,1);
  
  //get likelihood for the angle between Z and H
  yy=TMath::Min(pow(vh.E(),2.0),pow(vz.E(),2.0))
    *(1-coszh)/pow(500.0, 2.0);
  lh3=get_likelihood_ZH(yy);


  //cal.sigma of higgs and Z
  hgammameas=hmeas*0.550190*sqrt(pow(1/sqrt(hvect[0].E()),2.0)
				 +pow(1/sqrt(hvect[1].E()),2.0));
  zgammameas=zmeas*0.562543*sqrt(pow(1/sqrt(zvect[0].E()),2.0)
				 +pow(1/sqrt(zvect[1].E()),2.0));

  //now calculate chi square!!
  double chi1=pow(mh-hmeas,2.0)/pow(hgammameas,2.0)+TMath::Log(2*TMath::Pi()*hgammameas*hgammameas);
  double chi2=pow(mz-zmeas,2.0)/pow(zgammameas,2.0)+TMath::Log(2*TMath::Pi()*zgammameas*zgammameas);
 
  if(lh2==0.0) return 1.0e+20;

  // cout << "mass: " << chi1 << " " << chi2 << endl;
  // cout << "likelihood: " << -2*TMath::Log(lh1) << " "
  //      << -2*TMath::Log(lh2) << " " << -2*TMath::Log(lh3) << " " << endl;
  //double chiall=chi1+chi2-2*TMath::Log(lh2);   //-2*TMath::Log(lh3);
  double chiall=chi1+chi2-2*TMath::Log(lh2);   //-2*TMath::Log(lh3);
  //-2*TMath::Log(lh1)-2*TMath::Log(lh2)-2*TMath::Log(lh3);

  return chiall;
}

//costype: whether the jet is near to the lepton   1:nearW 2:farW 3:near Z 4:far Z 
//frametype: Higgs rest frame or W+Z rest frame    1:Higgs 2: W+Z
double CombinationSolver_bbbb::get_likelihood(double costheta, int costype){
  
  int nbins=pdf[costype]->GetNbinsX();
  double interval=pdf[costype]->GetBinWidth(1);

  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=0.0+i*interval && costheta<0.0+(i+1)*interval){
      val=pdf[costype]->GetBinContent(i+1);
      break;
    }
  }
  
  return val;
}

double CombinationSolver_bbbb::get_likelihood_ZH(double costheta){
  
  int nbins=pdf[2]->GetNbinsX();
  double interval=pdf[2]->GetBinWidth(1);
  
  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=0.0+i*interval && costheta<0.0+(i+1)*interval){
      val=pdf[2]->GetBinContent(i+1);
      //cout << "val_ZH: " << costheta << " " << val << endl;
      break;
    }
  }
  
  return val;
}

double CombinationSolver_bbbb::get_combination(jetdata data,int btagnum, int bjetid[], int *hid, int *zid){
  double chisq=1.0e20;

  if(btagnum<4) return chisq;
  
  TLorentzVector jets[4];
  jets[0].SetPxPyPzE(data.jetpx[bjetid[0]],
		     data.jetpy[bjetid[0]],
		     data.jetpz[bjetid[0]],
		     data.jete[bjetid[0]]);
  jets[1].SetPxPyPzE(data.jetpx[bjetid[1]],
		     data.jetpy[bjetid[1]],
		     data.jetpz[bjetid[1]],
		     data.jete[bjetid[1]]);
  jets[2].SetPxPyPzE(data.jetpx[bjetid[2]],
		     data.jetpy[bjetid[2]],
		     data.jetpz[bjetid[2]],
		     data.jete[bjetid[2]]);
  jets[3].SetPxPyPzE(data.jetpx[bjetid[3]],
		     data.jetpy[bjetid[3]],
		     data.jetpz[bjetid[3]],
		     data.jete[bjetid[3]]);

  vector<TLorentzVector> hvect;
  vector<TLorentzVector> zvect;
  for(Int_t i=0;i<btagnum;i++){
    for(Int_t j=i+1;j<btagnum;j++){
      hvect.push_back(jets[i]);
      hvect.push_back(jets[j]);
      
      //calculate others
      int k,l;
      switch(i){
      case 0:
	switch(j){
	case 1:
	  k=2;
	  l=3;
	  break;
	case 2:
	  k=1;
	  l=3;
	  break;
	case 3:
	  k=1;
	  l=2;
	  break;
	}
	break;
      case 1:
	switch(j){
	case 2:
	  k=0;
	  l=3;
	  break;
	case 3:
	  k=0;
	  l=2;
	  break;
	}
	break;
      case 2:
	switch(j){
	case 3:
	  k=0;
	  l=1;
	  break;
	}
	break;
      }
      
      //cout << "loop: " << ljetid[i] << " " << ljetid[j]
      //	    << " " << ljetid[k] << " " << ljetid[l] << endl;
     
      zvect.push_back(jets[k]);
      zvect.push_back(jets[l]);
      //get chi square
      double tmpchi=cal_chisq(hvect,zvect);
      if(chisq>tmpchi){
	chisq=tmpchi;
	hid[0]=bjetid[i];
	hid[1]=bjetid[j];
	zid[0]=bjetid[k];
	zid[1]=bjetid[l];
	//cout << "tmpchi: " << tmpchi << endl;
      }
      
      hvect.clear();
      zvect.clear();
    }
  }
  
  return chisq;
}

//combination solver for bbb case
CombinationSolver_bbb::CombinationSolver_bbb(string fstr){
  string ffstr1= fstr +"1.root";
  string ffstr2= fstr +"2.root";
  string ffstr3= fstr +"_bbbb.root";
  
  fpdf[0]=new TFile(ffstr1.c_str());
  fpdf[1]=new TFile(ffstr2.c_str());
  fpdf[2]=new TFile(ffstr3.c_str());

  pdf[0][0]=(TH1F*)fpdf[0]->Get("deltacosthetalepj3;1")->Clone("pdf11");
  pdf[1][0]=(TH1F*)fpdf[0]->Get("deltacosthetalepj4;1")->Clone("pdf21");
  pdf[2][0]=(TH1F*)fpdf[0]->Get("deltacosthetalepj5;1")->Clone("pdf31");
  pdf[3][0]=(TH1F*)fpdf[0]->Get("deltacosthetalepj6;1")->Clone("pdf41");
  pdf[0][1]=(TH1F*)fpdf[1]->Get("deltacosthetalepj3;1")->Clone("pdf12");
  pdf[1][1]=(TH1F*)fpdf[1]->Get("deltacosthetalepj4;1")->Clone("pdf22");
  pdf[2][1]=(TH1F*)fpdf[1]->Get("deltacosthetalepj5;1")->Clone("pdf32");
  pdf[3][1]=(TH1F*)fpdf[1]->Get("deltacosthetalepj6;1")->Clone("pdf42");

  pdf[0][2]=(TH1F*)fpdf[2]->Get("yZH1;1")->Clone("pdf13");
  pdf[1][2]=(TH1F*)fpdf[2]->Get("yZH2;1")->Clone("pdf23");
  pdf[2][2]=(TH1F*)fpdf[2]->Get("yZH3;1")->Clone("pdf33");

  pdf[3][2]=(TH1F*)fpdf[1]->Get("deltaRljlj1")->Clone("pdfWZ");

  //normalize histograms
  double weight=1.0;
  for(Int_t i=0;i<4;i++){
    for(Int_t j=0;j<3;j++){
      //if(i==3 && j==2) continue;
      weight=pdf[i][j]->Integral();
      pdf[i][j]->Scale(1/weight);
    }
  }
  
  cout << "constructor ok" << endl;
  return;
}

CombinationSolver_bbb::~CombinationSolver_bbb(){
  fpdf[0]->Close();
  fpdf[1]->Close();
  fpdf[2]->Close();
}

double CombinationSolver_bbb::cal_chisq(vector<TLorentzVector> hvect,
					vector<TLorentzVector> zvect,
					vector<TLorentzVector> wvect,
					TLorentzVector lep, TLorentzVector neu
					){
  
  //make LorentzVector for lepton and jets
  TLorentzVector hjet[2];
  hjet[0]=hvect[0];
  hjet[1]=hvect[1];

  TLorentzVector zjet[2];
  zjet[0]=zvect[0];
  zjet[1]=zvect[1];

  TLorentzVector wjet[2];
  wjet[0]=wvect[0];
  wjet[1]=wvect[1];

  TLorentzVector wlep[2];
  wlep[0]=lep;
  wlep[1]=lep;
  
  //make LorentzVector for W and Z
  TLorentzVector vh(hvect[0].Px()+hvect[1].Px(),
		    hvect[0].Py()+hvect[1].Py(),
		    hvect[0].Pz()+hvect[1].Pz(),
		    hvect[0].E()+hvect[1].E());
  
  TLorentzVector vz(zvect[0].Px()+zvect[1].Px(),
		    zvect[0].Py()+zvect[1].Py(),
		    zvect[0].Pz()+zvect[1].Pz(),
		    zvect[0].E()+zvect[1].E());

  TLorentzVector vw=wjet[0]+wjet[1];

  //boost in expected Higgs Frame
  TLorentzVector b1(-(wvect[0].Px()+wvect[1].Px()+lep.Px()+neu.Px()),
		    -(wvect[0].Py()+wvect[1].Py()+lep.Py()+neu.Py()),
		    -(wvect[0].Pz()+wvect[1].Pz()+lep.Pz()+neu.Pz()),
		    wvect[0].E()+wvect[1].E()+lep.E()+neu.E());
  
  //boost in expected wrong combination
  TLorentzVector b2(-(zvect[0].Px()+zvect[1].Px()+lep.Px()+neu.Px()),
		    -(zvect[0].Py()+zvect[1].Py()+lep.Py()+neu.Py()),
		    -(zvect[0].Pz()+zvect[1].Pz()+lep.Pz()+neu.Pz()),
		    zvect[0].E()+zvect[1].E()+lep.E()+neu.E());
  
  //start to calculate chi square
  //first, calculate Higgs and Z mass
  double mh=vh.M();
  double mh2=b1.M();
  double mz=vz.M();
  double mw=vw.M();

  //second, cal angle
  double coszh,coswz;
  coszh=TMath::Cos(vz.Angle(vh.Vect()));  //cos in Higgs rest frame

  double goodcos[2],badcos[2];
  //boost wvect with wrong rest frame
  wjet[0].Boost(b2.BoostVector());
  wjet[1].Boost(b2.BoostVector());
  wlep[0].Boost(b2.BoostVector());
  badcos[0]=TMath::Cos(wjet[0].Angle(wlep[0].Vect()));
  badcos[1]=TMath::Cos(wjet[1].Angle(wlep[0].Vect()));


  //boost zvect with wrong rest frame
  zjet[0].Boost(b2.BoostVector());
  zjet[1].Boost(b2.BoostVector());
  wlep[1].Boost(b2.BoostVector());
  goodcos[0]=TMath::Cos(wjet[0].Angle(wlep[1].Vect()));
  goodcos[1]=TMath::Cos(wjet[1].Angle(wlep[1].Vect()));

  coswz=TMath::Cos(wjet[0].Angle(wjet[0].Vect()));

  double lh1,lh2,lh3,lh4,lh5,lh6;

  //get likelihood value of angle
  double yy=0.0;
  yy=TMath::Min(pow(hjet[0].E(),2.0),pow(hjet[1].E(),2.0))
    *(1-TMath::Cos(hjet[0].Angle(hjet[1].Vect())))/pow(500.0, 2.0);
  lh1=get_likelihood1(yy,0);
  yy=TMath::Min(pow(zjet[0].E(),2.0),pow(zjet[1].E(),2.0))
    *(1-TMath::Cos(zjet[0].Angle(zjet[1].Vect())))/pow(500.0, 2.0);
  lh2=get_likelihood1(yy,1);
  
  //get likelihood for the angle between Z and H
  yy=TMath::Min(pow(vh.E(),2.0),pow(vz.E(),2.0))
    *(1-coszh)/pow(500.0, 2.0);
  lh3=get_likelihood_ZH(yy);

  //get likelihood for angle
  if(badcos[0]>badcos[1]){
    lh4=get_likelihood2(badcos[0],0,1);
  }else{
    lh4=get_likelihood2(badcos[1],0,1);
  }
  
  if(goodcos[0]>goodcos[1]){
    lh5=get_likelihood2(goodcos[1],3,1);
  }else{
    lh5=get_likelihood2(goodcos[0],3,1);
  }
  
  lh6=get_likelihood_WZ(coswz);

  //now calculate chi square!!
  double chi1=pow(mh-hmeas,2.0)/pow(hgammameas,2.0)+TMath::Log(2*TMath::Pi()*hgammameas*hgammameas);
  double chi2=pow(mz-zmeas,2.0)/pow(zgammameas,2.0)+TMath::Log(2*TMath::Pi()*zgammameas*zgammameas);
  double chi3=2*TMath::Log((1/TMath::Pi())*wgamma*mmw
			   /(pow(mw*mw-mmw*mmw,2.0)+mmw*mmw*wgamma*wgamma));
  double chi4=pow(mh2-hmeas2,2.0)/pow(hgammameas2,2.0)+TMath::Log(2*TMath::Pi()*hgammameas2*hgammameas2);
 
  if(lh1*lh3==0.0) return 1.0e+20;

  // cout << "mass: " << chi1 << " " << chi2 << endl;
  // cout << "likelihood: " << -2*TMath::Log(lh1) << " "
  //      << -2*TMath::Log(lh2) << " " << -2*TMath::Log(lh3) << " " << endl;
  //double chiall=chi1+chi2-2*TMath::Log(lh2);   //-2*TMath::Log(lh3);
  double chiall=chi1+chi2-chi3+chi4-2*TMath::Log(lh1)-2*TMath::Log(lh3);

  return chiall;
}

//costype: whether the jet is near to the lepton   1:nearW 2:farW 3:near Z 4:far Z 
//frametype: Higgs rest frame or W+Z rest frame    1:Higgs 2: W+Z
double CombinationSolver_bbb::get_likelihood2(double costheta, int costype, int frametype){
  
  int nbins=pdf[costype][frametype]->GetNbinsX();
  double interval=pdf[costype][frametype]->GetBinWidth(1);

  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=-1.0+i*interval && costheta<-1.0+(i+1)*interval){
      val=pdf[costype][frametype]->GetBinContent(i+1);
      //cout << "val_Wjet: " << costheta << " " << val << endl;
      break;
    }
  }

  return val;
}

//costype: whether the jet is near to the lepton   1:nearW 2:farW 3:near Z 4:far Z 
//frametype: Higgs rest frame or W+Z rest frame    1:Higgs 2: W+Z
double CombinationSolver_bbb::get_likelihood1(double costheta, int costype){
  
  int nbins=pdf[costype][2]->GetNbinsX();
  double interval=pdf[costype][2]->GetBinWidth(1);

  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=0.0+i*interval && costheta<0.0+(i+1)*interval){
      val=pdf[costype][2]->GetBinContent(i+1);
      //cout << "val_1: " << costheta << " " << val << endl;
       break;
    }
  }
  
  return val;
}

double CombinationSolver_bbb::get_likelihood_WZ(double costheta){
  
  int nbins=pdf[3][2]->GetNbinsX();
  double interval=pdf[3][2]->GetBinWidth(1);
  
  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=0.0+i*interval && costheta<0.0+(i+1)*interval){
      val=pdf[3][2]->GetBinContent(i+1);
      //cout << "val_ZH: " << costheta << " " << val << endl;
      break;
    }
  }
  
  return val;
}

double CombinationSolver_bbb::get_likelihood_ZH(double costheta){
  
  int nbins=pdf[2][2]->GetNbinsX();
  double interval=pdf[2][2]->GetBinWidth(1);
  
  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=0.0+i*interval && costheta<0.0+(i+1)*interval){
      val=pdf[2][2]->GetBinContent(i+1);
      //cout << "val_ZH: " << costheta << " " << val << endl;
      break;
    }
  }
  
  return val;
}

double CombinationSolver_bbb::get_combination(jetdata data, int btagnum, int nbtagnum,
					     int bjetid[], int ljetid[],
					     int *hid, int *zid){
  double chisq=1.0e20;
  
  if(btagnum<4) return chisq;
  
  TLorentzVector jets[4];
  jets[0].SetPxPyPzE(data.jetpx[bjetid[0]],
		     data.jetpy[bjetid[0]],
		     data.jetpz[bjetid[0]],
		     data.jete[bjetid[0]]);
  jets[1].SetPxPyPzE(data.jetpx[bjetid[1]],
		     data.jetpy[bjetid[1]],
		     data.jetpz[bjetid[1]],
		     data.jete[bjetid[1]]);
  jets[2].SetPxPyPzE(data.jetpx[bjetid[2]],
		     data.jetpy[bjetid[2]],
		     data.jetpz[bjetid[2]],
		     data.jete[bjetid[2]]);
  jets[3].SetPxPyPzE(data.jetpx[bjetid[3]],
		     data.jetpy[bjetid[3]],
		     data.jetpz[bjetid[3]],
		     data.jete[bjetid[3]]);

  int ljets[2];
  //make w jet
  int tmpnj=0;
  for(Int_t i=0;i<nbtagnum;i++){
    bool ljflg=true;
    for(Int_t j=0;j<btagnum;j++){
      if(bjetid[j]==ljetid[i]){
	ljflg=false;
	break;
      }
    }
    if(ljflg==true){
      ljets[tmpnj]=ljetid[i];
      tmpnj++;
    }
  }

  TLorentzVector wj[2];
  wj[0].SetPxPyPzE(data.jetpx[ljets[0]],
		   data.jetpy[ljets[0]],
		   data.jetpz[ljets[0]],
		   data.jete[ljets[0]]);

  wj[1].SetPxPyPzE(data.jetpx[ljets[1]],
		   data.jetpy[ljets[1]],
		   data.jetpz[ljets[1]],
		   data.jete[ljets[1]]);

  TLorentzVector lep(data.leppx[0],
		     data.leppy[0],
		     data.leppz[0],
		     data.lepe[0]);

  TLorentzVector neu(data.metpx,
		     data.metpy,
		     data.metpz,
		     sqrt(pow(data.metpx,2.0)+pow(data.metpy,2.0)+pow(data.metpz,2.0)));

  vector<TLorentzVector> hvect;
  vector<TLorentzVector> wvect;
  vector<TLorentzVector> zvect;

  wvect.push_back(wj[0]);
  wvect.push_back(wj[1]);

  for(Int_t i=0;i<btagnum;i++){
    for(Int_t j=i+1;j<btagnum;j++){
      hvect.push_back(jets[i]);
      hvect.push_back(jets[j]);
      
      //calculate others
      int k,l;
      switch(i){
      case 0:
	switch(j){
	case 1:
	  k=2;
	  l=3;
	  break;
	case 2:
	  k=1;
	  l=3;
	  break;
	case 3:
	  k=1;
	  l=2;
	  break;
	}
	break;
      case 1:
	switch(j){
	case 2:
	  k=0;
	  l=3;
	  break;
	case 3:
	  k=0;
	  l=2;
	  break;
	}
	break;
      case 2:
	switch(j){
	case 3:
	  k=0;
	  l=1;
	  break;
	}
	break;
      }
 
      //cout << "loop: " << bjetid[i] << " " << bjetid[j]
      //	    << " " << bjetid[k] << " " << bjetid[l] << endl;
     
      zvect.push_back(jets[k]);
      zvect.push_back(jets[l]);
      //get chi square
      double tmpchi=cal_chisq(hvect,zvect,wvect,lep,neu);
      if(chisq>tmpchi){
	chisq=tmpchi;
	hid[0]=bjetid[i];
	hid[1]=bjetid[j];
	zid[0]=bjetid[k];
	zid[1]=bjetid[l];
	//cout << "tmpchi: " << tmpchi << endl;
      }
      
      hvect.clear();
      zvect.clear();
    }
  }
  
  return chisq;
}

//combination solver for bbbb of all hadronic case
CombinationSolver_bbbb_allhad::CombinationSolver_bbbb_allhad(string fstr){
  string ffstr1= fstr +"_bbbb.root";
  string ffstr2= fstr +"_allhad_bbbb.root";
  
  fpdf[0]=new TFile(ffstr1.c_str());
  fpdf[1]=new TFile(ffstr2.c_str());

  //use for b jets
  pdf[0][0]=(TH1F*)fpdf[0]->Get("yZH1;1")->Clone("pdf11");
  pdf[1][0]=(TH1F*)fpdf[0]->Get("yZH2;1")->Clone("pdf21");
  pdf[2][0]=(TH1F*)fpdf[0]->Get("yZH3;1")->Clone("pdf31");

  pdf[0][1]=(TH1F*)fpdf[1]->Get("deltacosthetabjbj21;1")->Clone("pdf12");
  pdf[1][1]=(TH1F*)fpdf[1]->Get("deltacosthetabjbj22;1")->Clone("pdf22");
  pdf[2][1]=(TH1F*)fpdf[1]->Get("deltacosthetabjbj23;1")->Clone("pdf32");
  pdf[3][1]=(TH1F*)fpdf[1]->Get("deltacosthetabjbj24;1")->Clone("pdf42");
  pdf[4][1]=(TH1F*)fpdf[1]->Get("deltacosthetabjbj25;1")->Clone("pdf52");
  pdf[5][1]=(TH1F*)fpdf[1]->Get("deltacosthetabjbj26;1")->Clone("pdf62");

  pdf[0][2]=(TH1F*)fpdf[1]->Get("deltacosthetabjbj1;1")->Clone("pdf13");
  pdf[1][2]=(TH1F*)fpdf[1]->Get("deltacosthetabjbj2;1")->Clone("pdf23");

  //normalize histograms
  double weight=1.0;
  for(Int_t i=0;i<6;i++){
    for(Int_t j=0;j<2;j++){
      if(i>=3 && j==0) continue;
      if(i>=2 && j==2) continue;
      weight=pdf[i][j]->Integral();
      pdf[i][j]->Scale(1/weight);
    }
  }
  
  cout << "constructor ok" << endl;
  return;
}

CombinationSolver_bbbb_allhad::~CombinationSolver_bbbb_allhad(){
  fpdf[0]->Close();
  fpdf[1]->Close();
  return;
}

double CombinationSolver_bbbb_allhad::cal_chisq(vector<TLorentzVector> hvect,
						vector<TLorentzVector> zvect,
						vector<TLorentzVector> w1vect,
						vector<TLorentzVector> w2vect
						){
  
  //make LorentzVector for lepton and jets
  TLorentzVector hjet[2];
  hjet[0]=hvect[0];
  hjet[1]=hvect[1];

  TLorentzVector zjet[2];
  zjet[0]=zvect[0];
  zjet[1]=zvect[1];

  //wjet1[0] should be the highest energy jet always
  TLorentzVector wjet1[2];
  wjet1[0]=w1vect[0];
  wjet1[1]=w1vect[1];

  TLorentzVector wjet2[2];
  wjet2[0]=w2vect[0];
  wjet2[1]=w2vect[1];
  
  //make LorentzVector for W and Z
  TLorentzVector vh(hvect[0].Px()+hvect[1].Px(),
		    hvect[0].Py()+hvect[1].Py(),
		    hvect[0].Pz()+hvect[1].Pz(),
		    hvect[0].E()+hvect[1].E());
  
  TLorentzVector vz(zvect[0].Px()+zvect[1].Px(),
		    zvect[0].Py()+zvect[1].Py(),
		    zvect[0].Pz()+zvect[1].Pz(),
		    zvect[0].E()+zvect[1].E());

  TLorentzVector vw1=wjet1[0]+wjet1[1];
  TLorentzVector vw2=wjet2[0]+wjet2[1];

  //boost in 4 vector rest frame
  TLorentzVector b1(-(w1vect[0].Px()+w1vect[1].Px()+w2vect[0].Px()+w2vect[1].Px()),
		    -(w1vect[0].Py()+w1vect[1].Py()+w2vect[0].Py()+w2vect[1].Py()),
		    -(w1vect[0].Pz()+w1vect[1].Pz()+w2vect[0].Pz()+w2vect[1].Pz()),
		    w1vect[0].E()+w1vect[1].E()+w2vect[0].E()+w2vect[1].E());
  
  //start to calculate chi square
  //first, calculate Higgs and Z mass
  double mh=vh.M();
  double mh2=b1.M();
  double mz=vz.M();
  double mw1=vw1.M();
  double mw2=vw2.M();

  //second, cal angle
  double coszh,coswz;
  coszh=TMath::Cos(vz.Angle(vh.Vect()));  //cos in Higgs rest frame

  double coscos[4];
  //boost wvect with wrong rest frame
  wjet1[0].Boost(b1.BoostVector());
  wjet1[1].Boost(b1.BoostVector());
  wjet2[0].Boost(b1.BoostVector());
  wjet2[1].Boost(b1.BoostVector());

  //cal costheta
  coscos[0]=TMath::Min(TMath::Cos(wjet1[0].Angle(wjet2[0].Vect())), TMath::Cos(wjet1[0].Angle(wjet2[1].Vect())));
  coscos[1]=TMath::Max(TMath::Cos(wjet1[0].Angle(wjet2[0].Vect())), TMath::Cos(wjet1[0].Angle(wjet2[0].Vect())));
  coscos[2]=TMath::Cos(wjet1[0].Angle(wjet1[1].Vect()));
  
  double lh1,lh2,lh3,lh4,lh5,lh6;

  //get likelihood value of angle
  double yy=0.0;
  yy=TMath::Min(pow(hjet[0].E(),2.0),pow(hjet[1].E(),2.0))
    *(1-TMath::Cos(hjet[0].Angle(hjet[1].Vect())))/pow(500.0, 2.0);
  lh1=get_likelihood1(yy,0);
  yy=TMath::Min(pow(zjet[0].E(),2.0),pow(zjet[1].E(),2.0))
    *(1-TMath::Cos(zjet[0].Angle(zjet[1].Vect())))/pow(500.0, 2.0);
  lh2=get_likelihood1(yy,1);
  
  //get likelihood for the angle between Z and H
  yy=TMath::Min(pow(vh.E(),2.0),pow(vz.E(),2.0))
    *(1-coszh)/pow(500.0, 2.0);
  lh3=get_likelihood_ZH(yy);

  //get likelihood for angle
  lh4=get_likelihood2(coscos[0],0);
  lh5=get_likelihood2(coscos[2],4);

  //cal.sigma of higgs and Z
  hgammameas=hmeas*(0.577506*((1/sqrt(hjet[0].E()))+1/sqrt(hjet[1].E())));
  zgammameas=zmeas*(0.508333*((1/sqrt(zjet[0].E()))+1/sqrt(zjet[1].E())));

  //now calculate chi square!!
  double chi1=pow(mh-hmeas,2.0)/pow(hgammameas,2.0)+TMath::Log(2*TMath::Pi()*hgammameas*hgammameas);
  double chi2=pow(mz-zmeas,2.0)/pow(zgammameas,2.0)+TMath::Log(2*TMath::Pi()*zgammameas*zgammameas);
  double chi3=2*TMath::Log((1/TMath::Pi())*wgamma*mmw
			   /(pow(mw1*mw1-mmw*mmw,2.0)+mmw*mmw*wgamma*wgamma));
  double chi4=2*TMath::Log((1/TMath::Pi())*wgamma*mmw
			   /(pow(mw2*mw2-mmw*mmw,2.0)+mmw*mmw*wgamma*wgamma));
  double chi5=pow(mh2-hmeas2,2.0)/pow(hgammameas2,2.0)+TMath::Log(2*TMath::Pi()*hgammameas2*hgammameas2);

  double masym=(mh-mz)/(mh+mz);
  double chi6=pow(masym-0.155365,2.0)/pow(0.157758,2.0)+TMath::Log(2*TMath::Pi()*0.157758*0.157758);
 
  //if(lh1*lh3==0.0) return 1.0e+20;

  // cout << "mass: " << chi1 << " " << chi2 << endl;
  // cout << "likelihood: " << -2*TMath::Log(lh1) << " "
  //      << -2*TMath::Log(lh2) << " " << -2*TMath::Log(lh3) << " " << endl;
  //double chiall=chi1+chi2-2*TMath::Log(lh2);   //-2*TMath::Log(lh3);


  double chiall=chi1+chi2-chi3;   //+chi5-2*TMath::Log(lh2)-2*TMath::Log(lh5);

  //double chiall=chi1+chi2-chi3+chi5-2*TMath::Log(lh2)-2*TMath::Log(lh5);  //chi4 should not be used for the chi square!!!!!  good one


  return chiall;
}

//costype: whether the jet is near to the lepton   1:nearW 2:farW 3:near Z 4:far Z 
//frametype: Higgs rest frame or W+Z rest frame    1:Higgs 2: W+Z
double CombinationSolver_bbbb_allhad::get_likelihood2(double costheta, int costype){
  
  int nbins=pdf[costype][1]->GetNbinsX();
  double interval=pdf[costype][1]->GetBinWidth(1);

  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=-1.0+i*interval && costheta<-1.0+(i+1)*interval){
      val=pdf[costype][1]->GetBinContent(i+1);
      //cout << "val_Wjet: " << costheta << " " << val << endl;
      break;
    }
  }

  return val;
}

//costype: whether the jet is near to the lepton   1:nearW 2:farW 3:near Z 4:far Z 
//frametype: Higgs rest frame or W+Z rest frame    1:Higgs 2: W+Z
double CombinationSolver_bbbb_allhad::get_likelihood1(double costheta, int costype){
  
  int nbins=pdf[costype][0]->GetNbinsX();
  double interval=pdf[costype][0]->GetBinWidth(1);
  
  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=0.0+i*interval && costheta<0.0+(i+1)*interval){
      val=pdf[costype][0]->GetBinContent(i+1);
      //cout << "val_1: " << costheta << " " << val << endl;
      break;
    }
  }
  
  return val;
}

double CombinationSolver_bbbb_allhad::get_likelihood_ZH(double costheta){
  
  int nbins=pdf[2][0]->GetNbinsX();
  double interval=pdf[2][0]->GetBinWidth(1);
  
  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=0.0+i*interval && costheta<0.0+(i+1)*interval){
      val=pdf[2][0]->GetBinContent(i+1);
      //cout << "val_ZH: " << costheta << " " << val << endl;
      break;
    }
  }
  
  return val;
}

double CombinationSolver_bbbb_allhad::get_combination(jetdata data, int btagnum, int nbtagnum,
						      int bjetid[], int ljetid[],
						      int *hid, int *zid, int *w1id, int *w2id){
  double chisq=1.0e20;
  
  if(btagnum<4) return chisq;
  
  TLorentzVector jets[4];
  jets[0].SetPxPyPzE(data.jetpx[bjetid[0]],
		     data.jetpy[bjetid[0]],
		     data.jetpz[bjetid[0]],
		     data.jete[bjetid[0]]);
  jets[1].SetPxPyPzE(data.jetpx[bjetid[1]],
		     data.jetpy[bjetid[1]],
		     data.jetpz[bjetid[1]],
		     data.jete[bjetid[1]]);
  jets[2].SetPxPyPzE(data.jetpx[bjetid[2]],
		     data.jetpy[bjetid[2]],
		     data.jetpz[bjetid[2]],
		     data.jete[bjetid[2]]);
  jets[3].SetPxPyPzE(data.jetpx[bjetid[3]],
		     data.jetpy[bjetid[3]],
		     data.jetpz[bjetid[3]],
		     data.jete[bjetid[3]]);
  
  TLorentzVector ljets[4];
  ljets[0].SetPxPyPzE(data.jetpx[ljetid[0]],
		      data.jetpy[ljetid[0]],
		      data.jetpz[ljetid[0]],
		      data.jete[ljetid[0]]);
  ljets[1].SetPxPyPzE(data.jetpx[ljetid[1]],
		      data.jetpy[ljetid[1]],
		      data.jetpz[ljetid[1]],
		      data.jete[ljetid[1]]);
  ljets[2].SetPxPyPzE(data.jetpx[ljetid[2]],
		      data.jetpy[ljetid[2]],
		      data.jetpz[ljetid[2]],
		      data.jete[ljetid[2]]);
  ljets[3].SetPxPyPzE(data.jetpx[ljetid[3]],
		      data.jetpy[ljetid[3]],
		      data.jetpz[ljetid[3]],
		      data.jete[ljetid[3]]);
  
  
  vector<TLorentzVector> hvect;
  vector<TLorentzVector> w1vect;
  vector<TLorentzVector> w2vect;
  vector<TLorentzVector> zvect;
  
  //first decide b combination
  //set w vectors temporally
  w1vect.push_back(ljets[0]);
  w1vect.push_back(ljets[1]);
  w2vect.push_back(ljets[2]);
  w2vect.push_back(ljets[3]);
  w1id[0]=ljetid[0];
  w1id[1]=ljetid[1];
  w2id[0]=ljetid[2];
  w2id[1]=ljetid[3];

  for(Int_t i1=0;i1<btagnum;i1++){
    for(Int_t j1=i1+1;j1<btagnum;j1++){
      hvect.push_back(jets[i1]);
      hvect.push_back(jets[j1]);
      
      //calculate others
      int k1,l1;
      switch(i1){
      case 0:
	switch(j1){
	case 1:
	  k1=2;
	  l1=3;
	  break;
	case 2:
	  k1=1;
	  l1=3;
	  break;
	case 3:
	  k1=1;
	  l1=2;
	  break;
	}
	break;
      case 1:
	switch(j1){
	case 2:
	  k1=0;
	  l1=3;
	  break;
	case 3:
	  k1=0;
	  l1=2;
	  break;
	}
	break;
      case 2:
	switch(j1){
	case 3:
	  k1=0;
	  l1=1;
	  break;
	}
	break;
      }
      
      //cout << "loop: " << bjetid[i] << " " << bjetid[j]
      //	    << " " << bjetid[k] << " " << bjetid[l] << endl;
      
      zvect.push_back(jets[k1]);
      zvect.push_back(jets[l1]);
      
      //get chi square
      double tmpchi=cal_chisq(hvect,zvect,w1vect,w2vect);
      if(chisq>tmpchi){
	chisq=tmpchi;
	hid[0]=bjetid[i1];
	hid[1]=bjetid[j1];
	zid[0]=bjetid[k1];
	zid[1]=bjetid[l1];
	//cout << "tmpchi: " << tmpchi << endl;
      }
      
      hvect.clear();
      zvect.clear();
    }
  }
  
  //second decide w jets
  //set h and z correctly
  //set bjets again
  jets[0].SetPxPyPzE(data.jetpx[hid[0]],
		     data.jetpy[hid[0]],
		     data.jetpz[hid[0]],
		     data.jete[hid[0]]);
  jets[1].SetPxPyPzE(data.jetpx[hid[1]],
		     data.jetpy[hid[1]],
		     data.jetpz[hid[1]],
		     data.jete[hid[1]]);
  jets[2].SetPxPyPzE(data.jetpx[zid[0]],
		     data.jetpy[zid[0]],
		     data.jetpz[zid[0]],
		     data.jete[zid[0]]);
  jets[3].SetPxPyPzE(data.jetpx[zid[1]],
		     data.jetpy[zid[1]],
		     data.jetpz[zid[1]],
		     data.jete[zid[1]]);
  hvect.clear();
  zvect.clear();
  hvect.push_back(jets[0]);
  hvect.push_back(jets[1]);
  zvect.push_back(jets[2]);
  zvect.push_back(jets[3]);
  w1vect.clear();
  w2vect.clear();
  int i2=0;
  for(int j2=1;j2<nbtagnum;j2++){
    w1vect.push_back(ljets[i2]);
    w1vect.push_back(ljets[j2]);
    
    //calculate others
    int k2,l2;
    switch(j2){
    case 1:
      k2=2;
      l2=3;	
      break;
    case 2:
      k2=1;
      l2=3;	
      break;
    case 3:
      k2=1;
      l2=2;	
      break;
    }
    
    //cout << "loop: " << bjetid[i] << " " << bjetid[j]
    //	    << " " << bjetid[k] << " " << bjetid[l] << endl;
    
    w2vect.push_back(ljets[k2]);
    w2vect.push_back(ljets[l2]);
    //get chi square
    double tmpchi=cal_chisq(hvect,zvect,w1vect,w2vect);
    if(chisq>tmpchi){
      chisq=tmpchi;
      w1id[0]=ljetid[i2];
      w1id[1]=ljetid[j2];
      w2id[0]=ljetid[k2];
      w2id[1]=ljetid[l2];
      //cout << "tmpchi: " << tmpchi << endl;
    }
    w1vect.clear();
    w2vect.clear();
  }
	    
  return chisq;
}

//combination solver for bbbb of all hadronic case
CombinationSolver_dilepton::CombinationSolver_dilepton(string fstr){
  string ffstr= fstr +"_allhad_bbbb.root";
  
  fpdf=new TFile(ffstr.c_str());

  pdf[0]=(TH1F*)fpdf->Get("deltacosthetabjbj21;1")->Clone("pdf12");
  pdf[1]=(TH1F*)fpdf->Get("deltacosthetabjbj22;1")->Clone("pdf22");
  pdf[2]=(TH1F*)fpdf->Get("deltacosthetabjbj23;1")->Clone("pdf32");
  pdf[3]=(TH1F*)fpdf->Get("deltacosthetabjbj24;1")->Clone("pdf42");
  pdf[4]=(TH1F*)fpdf->Get("deltacosthetabjbj25;1")->Clone("pdf52");
  pdf[5]=(TH1F*)fpdf->Get("deltacosthetabjbj26;1")->Clone("pdf62");

  //normalize histograms
  double weight=1.0;
  for(Int_t i=0;i<6;i++){
    weight=pdf[i]->Integral();
    pdf[i]->Scale(1/weight);
  }
  
  cout << "constructor ok" << endl;
  return; 
}

CombinationSolver_dilepton::~CombinationSolver_dilepton(){
  fpdf->Close();
  return;
}

double CombinationSolver_dilepton::cal_chisq(vector<TLorentzVector> w1vect,
					     vector<TLorentzVector> w2vect
					     ){
  
  //wjet1[0] should be the highest energy jet always
  TLorentzVector wjet1[2];
  wjet1[0]=w1vect[0];
  wjet1[1]=w1vect[1];
  
  TLorentzVector wjet2[2];
  wjet2[0]=w2vect[0];
  wjet2[1]=w2vect[1];
  
  TLorentzVector vw1=wjet1[0]+wjet1[1];
  TLorentzVector vw2=wjet2[0]+wjet2[1];
  
  //boost in 4 vector rest frame
  TLorentzVector b1(-(w1vect[0].Px()+w1vect[1].Px()+w2vect[0].Px()+w2vect[1].Px()),
		    -(w1vect[0].Py()+w1vect[1].Py()+w2vect[0].Py()+w2vect[1].Py()),
		    -(w1vect[0].Pz()+w1vect[1].Pz()+w2vect[0].Pz()+w2vect[1].Pz()),
		    w1vect[0].E()+w1vect[1].E()+w2vect[0].E()+w2vect[1].E());
  
  //start to calculate chi square
  //first, calculate Higgs and Z mass
  double mh=b1.M();
  double mw1=vw1.M();
  double mw2=vw2.M();

  //second, cal angle
  double coscos[4];
  //boost wvect with wrong rest frame
  wjet1[0].Boost(b1.BoostVector());
  wjet1[1].Boost(b1.BoostVector());
  wjet2[0].Boost(b1.BoostVector());
  wjet2[1].Boost(b1.BoostVector());
  
  //cal costheta
  coscos[0]=TMath::Min(TMath::Cos(wjet1[0].Angle(wjet2[0].Vect())), TMath::Cos(wjet1[0].Angle(wjet2[1].Vect())));
  coscos[1]=TMath::Max(TMath::Cos(wjet1[0].Angle(wjet2[0].Vect())), TMath::Cos(wjet1[0].Angle(wjet2[1].Vect())));
  coscos[2]=TMath::Cos(wjet1[0].Angle(wjet1[1].Vect()));
  
  double lh1,lh2;
  //get likelihood for angle
  lh1=get_likelihood2(coscos[0],0);
  lh2=get_likelihood2(coscos[2],4);
  
  //now calculate chi square!!
  double chi1=2*TMath::Log((1/TMath::Pi())*wgamma*mmw
			   /(pow(mw1*mw1-mmw*mmw,2.0)+mmw*mmw*wgamma*wgamma));
  double chi2=2*TMath::Log((1/TMath::Pi())*wgamma*mmw
			   /(pow(mw2*mw2-mmw*mmw,2.0)+mmw*mmw*wgamma*wgamma));
  double chi3=pow(mh-hmeas,2.0)/pow(hgammameas,2.0)+TMath::Log(2*TMath::Pi()*hgammameas*hgammameas);
 
  if(lh1==0.0) return 1.0e+20;

  double chiall=-chi1+chi3-2*TMath::Log(lh1); 

  return chiall;
}

//costype: whether the jet is near to the lepton   1:nearW 2:farW 3:near Z 4:far Z 
//frametype: Higgs rest frame or W+Z rest frame    1:Higgs 2: W+Z
double CombinationSolver_dilepton::get_likelihood2(double costheta, int costype){
  
  int nbins=pdf[costype]->GetNbinsX();
  double interval=pdf[costype]->GetBinWidth(1);
  
  double val=1.0e-50;
  for(Int_t i=0;i<nbins;i++){
    if(costheta>=-1.0+i*interval && costheta<-1.0+(i+1)*interval){
      val=pdf[costype]->GetBinContent(i+1);
      //cout << "val_Wjet: " << costheta << " " << val << endl;
      break;
    }
  }
  
  return val;
}

double CombinationSolver_dilepton::get_combination(jetdata data, int nbtagnum, int ljetid[],
						   int *w1id, int *w2id){
  double chisq=1.0e20;
  
  if(nbtagnum<4) return chisq;
  
  TLorentzVector ljets[4];
  ljets[0].SetPxPyPzE(data.jetpx[ljetid[0]],
		      data.jetpy[ljetid[0]],
		      data.jetpz[ljetid[0]],
		      data.jete[ljetid[0]]);
  ljets[1].SetPxPyPzE(data.jetpx[ljetid[1]],
		      data.jetpy[ljetid[1]],
		      data.jetpz[ljetid[1]],
		      data.jete[ljetid[1]]);
  ljets[2].SetPxPyPzE(data.jetpx[ljetid[2]],
		      data.jetpy[ljetid[2]],
		      data.jetpz[ljetid[2]],
		      data.jete[ljetid[2]]);
  ljets[3].SetPxPyPzE(data.jetpx[ljetid[3]],
		      data.jetpy[ljetid[3]],
		      data.jetpz[ljetid[3]],
		      data.jete[ljetid[3]]);
  
  
  vector<TLorentzVector> w1vect;
  vector<TLorentzVector> w2vect;
  
  int i2=0;
  for(int j2=1;j2<nbtagnum;j2++){
    w1vect.push_back(ljets[i2]);
    w1vect.push_back(ljets[j2]);
    
    //calculate others
    int k2,l2;
    switch(j2){
    case 1:
      k2=2;
      l2=3;	
      break;
    case 2:
      k2=1;
      l2=3;	
      break;
    case 3:
      k2=1;
      l2=2;	
      break;
    }
    
    //cout << "loop: " << bjetid[i] << " " << bjetid[j]
    //	    << " " << bjetid[k] << " " << bjetid[l] << endl;
    
    w2vect.push_back(ljets[k2]);
    w2vect.push_back(ljets[l2]);
    //get chi square
    double tmpchi=cal_chisq(w1vect,w2vect);
    if(chisq>tmpchi){
      chisq=tmpchi;
      w1id[0]=ljetid[i2];
      w1id[1]=ljetid[j2];
      w2id[0]=ljetid[k2];
      w2id[1]=ljetid[l2];
      //cout << "tmpchi: " << tmpchi << endl;
    }
    w1vect.clear();
    w2vect.clear();
  }
  
  return chisq;
}

mvafile::mvafile(string fstr, int fieldtype){
  //define the rootfile for MVA training
  ftmva=new TFile(fstr.c_str(), "RECREATE");

  for(Int_t i=0;i<fieldtype;i++){
    switch(i){ //set tree
    case 0:
      tts[i]=new TTree("signal", "signal variables");
      break;
    case 1:
      tts[i]=new TTree("ttbar_lepjets", "ttbar lep+jets variables");
      break;
    case 2:
      tts[i]=new TTree("ttbar_allhadronic", "ttbar all hadronic variables");
      break;
    case 3:
      tts[i]=new TTree("ttbar_dilepton", "ttbar dilepton variables");
      break;
    case 4:
      tts[i]=new TTree("ttbarbb", "ttbar + bb variables");
      break;
    case 5:
      tts[i]=new TTree("ttz", "ttbar + Z variables");
      break;
    case 6:
      tts[i]=new TTree("tth", "ttbar + H variables");
      break;
    case 7:
      tts[i]=new TTree("zzh", "ZZ + H variables");
      break;
    case 8:
      tts[i]=new TTree("zzz", "ZZZ variables");
      break;
    }

    //set address
    tts[i]->Branch("mHbb", &var[0], "mHbb /D");
    tts[i]->Branch("mZbb", &var[1], "mZbb /D");
    tts[i]->Branch("mW1", &var[2], "mW1 /D");
    tts[i]->Branch("mW2", &var[3], "mW2 /D");
    tts[i]->Branch("deltaRHbb", &var[4], "deltaRHbb /D");
    tts[i]->Branch("deltaRZbb", &var[5], "deltaRZbb /D");
    tts[i]->Branch("deltaRW1", &var[6], "deltaRW1 /D");
    tts[i]->Branch("deltaRW2", &var[7], "deltaRW2 /D");
    tts[i]->Branch("costhetaHbb", &var[8], "costhetaHbb /D");
    tts[i]->Branch("costhetaZbb", &var[9], "costhetaZbb /D");
    tts[i]->Branch("costhetaW1", &var[10], "costhetaW1 /D");
    tts[i]->Branch("costhetaW2", &var[11], "costhetaW2 /D");
    tts[i]->Branch("chisq", &var[12], "chisq /D");
    tts[i]->Branch("y56", &var[13], "y56 /D");
    tts[i]->Branch("y67", &var[14], "y67 /D");
    tts[i]->Branch("y78", &var[15], "y78 /D");
    tts[i]->Branch("yj1j2", &var[16], "yj1j2 /D");
    tts[i]->Branch("yj1nearj", &var[17], "yj1nearj /D");
    tts[i]->Branch("yj2nearj", &var[18], "yj2nearj /D");
    tts[i]->Branch("ysum", &var[19], "ysum /D");
    tts[i]->Branch("moment0", &var[20], "moment0 /D");
    tts[i]->Branch("moment1", &var[21], "moment1 /D");
    tts[i]->Branch("moment2", &var[22], "moment2 /D");
    tts[i]->Branch("moment3", &var[23], "moment3 /D");
    tts[i]->Branch("moment4", &var[24], "moment4 /D");
    tts[i]->Branch("moment5", &var[25], "moment5 /D");
    tts[i]->Branch("sphericity", &var[26], "sphericity /D");
    tts[i]->Branch("aplanarity", &var[27], "aplanarity /D");
    tts[i]->Branch("Cvalue", &var[28], "Cvalue /D");
    tts[i]->Branch("Dvalue", &var[29], "Dvalue /D");
    tts[i]->Branch("masym1", &var[30], "masym1 /D");
    tts[i]->Branch("masym2", &var[31], "masym2 /D");
    tts[i]->Branch("m4vect", &var[32], "m4vect /D");
    tts[i]->Branch("Evisible", &var[33], "Evisible /D");
    tts[i]->Branch("bjet1px", &var[34], "bjet1px /D");
    tts[i]->Branch("bjet1py", &var[35], "bjet1py /D");
    tts[i]->Branch("bjet1pz", &var[36], "bjet1pz /D");
    tts[i]->Branch("bjet1e", &var[37], "bjet1e /D");
    tts[i]->Branch("bjet2px", &var[38], "bjet2px /D");
    tts[i]->Branch("bjet2py", &var[39], "bjet2py /D");
    tts[i]->Branch("bjet2pz", &var[40], "bjet2pz /D");
    tts[i]->Branch("bjet2e", &var[41], "bjet2e /D");
    tts[i]->Branch("bjet3px", &var[42], "bjet3px /D");
    tts[i]->Branch("bjet3py", &var[43], "bjet3py /D");
    tts[i]->Branch("bjet3pz", &var[44], "bjet3pz /D");
    tts[i]->Branch("bjet3e", &var[45], "bjet3e /D");
    tts[i]->Branch("bjet4px", &var[46], "bjet4px /D");
    tts[i]->Branch("bjet4py", &var[47], "bjet4py /D");
    tts[i]->Branch("bjet4pz", &var[48], "bjet4pz /D");
    tts[i]->Branch("bjet4e", &var[49], "bjet4e /D");
    tts[i]->Branch("ljet1px", &var[50], "ljet1px /D");
    tts[i]->Branch("ljet1py", &var[51], "ljet1py /D");
    tts[i]->Branch("ljet1pz", &var[52], "ljet1pz /D");
    tts[i]->Branch("ljet1e", &var[53], "ljet1e /D");
    tts[i]->Branch("ljet2px", &var[54], "ljet2px /D");
    tts[i]->Branch("ljet2py", &var[55], "ljet2py /D");
    tts[i]->Branch("ljet2pz", &var[56], "ljet2pz /D");
    tts[i]->Branch("ljet2e", &var[57], "ljet2e /D");
    tts[i]->Branch("ljet3px", &var[58], "ljet3px /D");
    tts[i]->Branch("ljet3py", &var[59], "ljet3py /D");
    tts[i]->Branch("ljet3pz", &var[60], "ljet3pz /D");
    tts[i]->Branch("ljet3e", &var[61], "ljet3e /D");
    tts[i]->Branch("ljet4px", &var[62], "ljet4px /D");
    tts[i]->Branch("ljet4py", &var[63], "ljet4py /D");
    tts[i]->Branch("ljet4pz", &var[64], "ljet4pz /D");
    tts[i]->Branch("ljet4e", &var[65], "ljet4e /D");
    tts[i]->Branch("lep1px", &var[66], "lep1px /D");
    tts[i]->Branch("lep1py", &var[67], "lep1py /D");
    tts[i]->Branch("lep1pz", &var[68], "lep1pz /D");
    tts[i]->Branch("lep1e", &var[69], "lep1e /D");
    tts[i]->Branch("lep2px", &var[70], "lep2px /D");
    tts[i]->Branch("lep2py", &var[71], "lep2py /D");
    tts[i]->Branch("lep2pz", &var[72], "lep2pz /D");
    tts[i]->Branch("lep2e", &var[73], "lep2e /D");
    tts[i]->Branch("metpx", &var[74], "metpx /D");
    tts[i]->Branch("metpy", &var[75], "metpy /D");
    tts[i]->Branch("metpz", &var[76], "metpz /D");
    tts[i]->Branch("mete", &var[77], "mete /D");
    tts[i]->Branch("m6vect", &var[78], "m6vect /D");
    tts[i]->Branch("mtwlnu", &var[79], "mtwlnu /D");
    tts[i]->Branch("costhetaWW", &var[80], "costhetaWW /D");
    tts[i]->Branch("costhetaZH", &var[81], "costhetaZH /D");
    tts[i]->Branch("costhetaZWW", &var[82], "costhetaZWW /D");
    tts[i]->Branch("costhetaHWW", &var[83], "costhetaHWW /D");
    tts[i]->Branch("yW1j1j2", &var[84], "yW1j1j2 /D");
    tts[i]->Branch("yW2j1j2", &var[85], "yW2j1j2 /D");
    tts[i]->Branch("sjvaluezj1", &var[86], "sjvaluezj1 /D");
    tts[i]->Branch("sjvaluezj2", &var[87], "sjvaluezj2 /D");
    tts[i]->Branch("sjvaluehj1", &var[88], "sjvaluehj1 /D");
    tts[i]->Branch("sjvaluehj2", &var[89], "sjvaluehj2 /D");
    tts[i]->Branch("sjvalueW1j1", &var[90], "sjvalueW1j1 /D");
    tts[i]->Branch("sjvalueW1j2", &var[91], "sjvalueW1j2 /D");
    tts[i]->Branch("sjvalueW2j1", &var[92], "sjvalueW2j1 /D");
    tts[i]->Branch("sjvalueW2j2", &var[93], "sjvalueW2j2 /D");
    tts[i]->Branch("btagzj1", &var[94], "btagzj1 /D");
    tts[i]->Branch("btagzj2", &var[95], "btagzj2 /D");
    tts[i]->Branch("btaghj1", &var[96], "btaghj1 /D");
    tts[i]->Branch("btaghj2", &var[97], "btaghj2 /D");
    tts[i]->Branch("btagW1j1", &var[98], "btagW1j1 /D");
    tts[i]->Branch("btagW1j2", &var[99], "btagW1j2 /D");
    tts[i]->Branch("btagW2j1", &var[100], "btagW2j1 /D");
    tts[i]->Branch("btagW2j2", &var[101], "btagW2j2 /D");
    tts[i]->Branch("ctagzj1", &var[102], "ctagzj1 /D");
    tts[i]->Branch("ctagzj2", &var[103], "ctagzj2 /D");
    tts[i]->Branch("ctaghj1", &var[104], "ctaghj1 /D");
    tts[i]->Branch("ctaghj2", &var[105], "ctaghj2 /D");
    tts[i]->Branch("ctagW1j1", &var[106], "ctagW1j1 /D");
    tts[i]->Branch("ctagW1j2", &var[107], "ctagW1j2 /D");
    tts[i]->Branch("ctagW2j1", &var[108], "ctagW2j1 /D");
    tts[i]->Branch("ctagW2j2", &var[109], "ctagW2j2 /D");

    tts[i]->Branch("disc_ttbar", &var[110], "disc_ttbar /D");
    tts[i]->Branch("disc_ttbb", &var[111], "disc_ttbb /D");
    tts[i]->Branch("coneEnergy1", &var[112], "coneEnergy1 /D");
    tts[i]->Branch("coneEnergy2", &var[113], "coneEnergy2 /D");
    tts[i]->Branch("coneEnergy3", &var[114], "coneEnergy3 /D");
    tts[i]->Branch("showerchi21", &var[115], "showerchi21 /D");
    tts[i]->Branch("showerchi22", &var[116], "showerchi22 /D");
    tts[i]->Branch("showerchi23", &var[117], "showerchi23 /D");
    tts[i]->Branch("absLength1", &var[118], "absLength1 /D");
    tts[i]->Branch("absLength2", &var[119], "absLength2 /D");
    tts[i]->Branch("absLength3", &var[120], "absLength3 /D");
    tts[i]->Branch("xl201", &var[121], "xl201 /D");
    tts[i]->Branch("xl202", &var[122], "xl202 /D");
    tts[i]->Branch("xl203", &var[123], "xl203 /D");
    tts[i]->Branch("lepcharge1", &var[124], "lepcharge1 /D");
    tts[i]->Branch("lepcharge2", &var[125], "lepcharge2 /D");
    tts[i]->Branch("lepcharge3", &var[126], "lepcharge3 /D");
    tts[i]->Branch("topveto1", &var[127], "topveto1 /D");
    tts[i]->Branch("topveto2", &var[128], "topveto2 /D");
    tts[i]->Branch("topveto3", &var[129], "topveto3 /D");
    tts[i]->Branch("likeli1", &var[130], "likeli1 /D");
    tts[i]->Branch("likeli2", &var[131], "likeli2 /D");
    tts[i]->Branch("likeli3", &var[132], "likeli3 /D");
    tts[i]->Branch("sigchisq1", &var[133], "sigchisq1 /D");
    tts[i]->Branch("sigchisq2", &var[134], "sigchisq2 /D");
    tts[i]->Branch("scale1", &var[135], "scale1 /D");
    tts[i]->Branch("scale2", &var[136], "scale2 /D");
    tts[i]->Branch("scale3", &var[137], "scale3 /D");
    tts[i]->Branch("scale4", &var[138], "scale4 /D");
    tts[i]->Branch("scale5", &var[139], "scale5 /D");
    tts[i]->Branch("scale6", &var[140], "scale6 /D");
    tts[i]->Branch("scale7", &var[141], "scale7 /D");
    tts[i]->Branch("scale8", &var[142], "scale8 /D");
    tts[i]->Branch("zzhchisq1", &var[143], "zzhchisq1 /D");
    tts[i]->Branch("zzhchisq2", &var[144], "zzhchisq2 /D");
    tts[i]->Branch("z1mass", &var[145], "z1mass /D");
    tts[i]->Branch("z2mass", &var[146], "z2mass /D");
    tts[i]->Branch("hmass", &var[147], "hmass /D");
    tts[i]->Branch("w1mass", &var[148], "w1mass /D");
    tts[i]->Branch("w2mass", &var[149], "w2mass /D");

    tts[i]->Branch("colzmass", &var[197], "colzmass /D");
    tts[i]->Branch("eventtype", &var[198], "eventtype /D");
    tts[i]->Branch("weight", &var[199], "weight /D");
  }

  return;
}

mvafile::~mvafile(){
  return;
}

void mvafile::fill(double *par, int field){
  for(int i=0;i<200;i++){
    var[i]=par[i];
  }
  
  tts[field]->Fill();
  
  return;
}

void mvafile::write(int fieldtype){
  
  for(int i=0;i<fieldtype;i++){
    tts[i]->Write();
  }

  return;
}

void mvafile::close(){
  ftmva->Close();
  
  return;
}

jetEnergyScale::jetEnergyScale(){
  //define TMVA reader
  //(1)2nd vettex
  reader[0]=new TMVA::Reader( "!Color:!Silent" );
  reader[0]->AddVariable( "jete", &var[0] );
  reader[0]->AddVariable( "jetmass", &var[10] );
  reader[0]->AddVariable( "vtxmass", &var[29] );
  //reader[0]->AddVariable( "ntrk+nneu", &var[28] );
  reader[0]->AddVariable( "ntrk", &var[1] );
  reader[0]->AddVariable( "nneu", &var[15] );
  reader[0]->AddVariable( "jetep", &var[11] );
  //reader[0]->AddVariable( "trkpsum", &var[18] );
  reader[0]->AddVariable( "trkptrelsum", &var[6] );
  //reader[0]->AddVariable( "neuEnergy", &var[7] );
  reader[0]->AddVariable( "neuptrelsum", &var[9] );
  // reader[0]->AddVariable( "vtxmom1_jete", &var[2] );
  // reader[0]->AddVariable( "trkmaxp/jete", &var[3] );
  // reader[0]->AddVariable( "trkpsum/ntrack", &var[4] );
  // reader[0]->AddVariable( "trkmaxptrel", &var[5] );
  //reader[0]->AddVariable( "hadem", &var[21] );
  // reader[0]->AddVariable( "jetneuEnergy", &var[7] );
  // reader[0]->AddVariable( "neumaxptrel", &var[8] );
  //reader[0]->AddVariable( "trkpsum", &var[18] );
  reader[0]->AddVariable( "Aa", &var[24] );
  reader[0]->AddVariable( "twomom", &var[30] );
  //reader[0]->AddVariable( "nearcos", &var[19] );
  //reader[0]->AddVariable( "nearjetmass", &var[36] );
  //reader[0]->AddVariable( "nearjete",  &var[31] );
  //reader[0]->AddVariable( "jete+nearjete",  &var[34] );
  //reader[0]->AddVariable( "jete-nearjete",  &var[35] );

  reader[0]->AddSpectator( "mc_e1",  &var[49] );
  //reader[0]->AddSpectator( "nearjete",  &var[31] );
  //reader[0]->AddSpectator( "nearjetmc",  &var[32] );
  //reader[0]->AddSpectator( "nearjetdr",  &var[33] );

  //reader[0]->AddSpectator( "ntrack",  &var[21] );
  //reader[0]->AddSpectator( "nneu",  &var[22] );
  //reader[0]->BookMVA( "BNNMLP2_1", "weights/TMVARegression_BNNMLP2_1.weights.xml" ); 
  //reader[0]->BookMVA( "BNNMLP1_1", "test/weights/TMVARegression_BNNMLP1_1.weights.xml" ); 
  //reader[0]->BookMVA( "MLPjec_1", "test/weights/TMVARegression_MLPjec_1.weights.xml" ); 
  reader[0]->BookMVA( "BDTGjec_1", "test/weights/TMVARegression_BDTGjec_1.weights.xml" ); 
  //reader[0]->BookMVA( "BDTG2_1", "weights/TMVARegression_BDTG2_1.weights.xml" ); 
 
  //(2)nearcos>=0.8 && b-tagging case
  //reader[0]->BookMVA( "BNNMLP2_2", "weights/TMVARegression_BNNMLP2_2.weights.xml" ); 
  //reader[0]->BookMVA( "BDTG", "weights/TMVARegression_BDTG.weights.xml" ); 
 
  //(3)nearcos<0.8 && jet<90 && b-tagging case
  //reader[0]->BookMVA( "BNNMLP1_1", "weights/TMVARegression_BNNMLP1_1.weights.xml" ); 

  //no vertex
  reader[1]=new TMVA::Reader( "!Color:!Silent" );
  //reader[1]->AddVariable( "jete", &var[0] );
  reader[1]->AddVariable( "jetmass", &var[10] );
  reader[1]->AddVariable( "ntrk", &var[1] );
  reader[1]->AddVariable( "nneu", &var[15] );
  //reader[1]->AddVariable( "jetep", &var[11] );
  //reader[1]->AddVariable( "trkptsum", &var[12] );
  //reader[1]->AddVariable( "trkmaxptrel", &var[5] );
  reader[1]->AddVariable( "trkptrelsum", &var[6] );
  //reader[1]->AddVariable( "neumaxptrel", &var[8] );
  reader[1]->AddVariable( "neuptrelsum", &var[9] );
  reader[1]->AddVariable( "Bjet", &var[20] );
  reader[1]->AddVariable( "Aa", &var[24] );
  //reader[1]->AddVariable( "nearcos", &var[19] );
  reader[1]->AddVariable( "nearjetmass", &var[36] );
  //reader[1]->AddVariable( "nearjete",  &var[31] );
  reader[1]->AddVariable( "jete+nearjete",  &var[34] );
  reader[1]->AddVariable( "jete-nearjete",  &var[35] );
  

  reader[1]->AddSpectator( "mc_e1",  &var[49] );
  reader[1]->AddSpectator( "nearjete",  &var[31] );
  reader[1]->AddSpectator( "nearjetmc",  &var[32] );
  reader[1]->AddSpectator( "nearjetdr",  &var[33] );

  //reader[1]->BookMVA( "BDTGjec_2", "test/weights/TMVARegression_BDTGjec_2.weights.xml" ); 
  
  return;
}

jetEnergyScale::~jetEnergyScale(){
  delete reader[0];
  delete reader[1];

  return;
}

Double_t jetEnergyScale::getScale(int type, Float_t *val, Int_t valtype){
  Double_t scale[2]={1.0,1.0};
  
  for(Int_t i=0;i<50;i++){
    var[i]=val[i];
  }
    
  switch(type){
  case 0:
    //scale=reader[0]->EvaluateRegression( "BNNMLP2_1" )[0];
    //scale=reader[0]->EvaluateRegression( "BDTG2_1" )[0];
    //scale=reader[0]->EvaluateRegression( "BNNMLP1_1" )[0];
    //scale=reader[0]->EvaluateRegression( "MLPjec_1" )[0];
    scale[0]=reader[0]->EvaluateRegression( "BDTGjec_1" )[0];
    //scale[0]=reader[0]->EvaluateRegression( "MLPjec_1" )[0];
    //scale[1]=reader[0]->EvaluateRegression( "MLPjec_1" )[1];
    break;
  case 1:
    //scale=reader[0]->EvaluateRegression( "BNNMLP2_2" )[0];
    //scale[0]=reader[1]->EvaluateRegression( "BDTGjec_2" )[0];
    //scale=reader[0]->EvaluateRegression( "BDTG" )[0];
    //scale[0]=reader[1]->EvaluateRegression( "MLPjec_2" )[0];
    //scale[1]=reader[1]->EvaluateRegression( "MLPjec_2" )[1];
    break;
  case 2:
    //scale=reader[0]->EvaluateRegression( "BNNMLP1_2" )[0];
    scale[0]=reader[0]->EvaluateRegression( "BDTG" )[0];
    break;
  case 3:
    scale[0]=reader[1]->EvaluateRegression( "BNNMLP4_1" )[0];
    break;
  case 4:
    scale[0]=reader[1]->EvaluateRegression( "BNNMLP4_2" )[0];
    break;
  case 5:
    scale[0]=reader[1]->EvaluateRegression( "BNNMLP4_1" )[0];
    break;

  }
  return scale[valtype];
}

expert_ttbar::expert_ttbar(int cattype){
  
  //cattype indicates the separated category types
  //0: all hadronic 4btag case
  //1: all hadronic 3btag(?) case
  //2: all hadronic 2btag+ctag(?) case
  //3-10: reserved
  //11: lep+jets, 4btag case
  //12: lep+jets, 3btag case
  //13: lep+jets, 2btag+ctag case
  //14: lep+jets, 2btag case
  //15-20: reserved
  //21: dilepton, 2btag case(reserved)
  type=cattype;
  
  switch(cattype){
  case 0:    //all hadronic
    reader[cattype]=new TMVA::Reader( "Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    reader[cattype]->AddVariable( "mW1", &var[2]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    //reader[cattype]->AddVariable( "moment2", &var[22]);
    //reader[cattype]->AddVariable( "Cvalue", &var[28]);
    //reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    //reader[cattype]->AddVariable( "y78", &var[15]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "sj1", &var[138]);
    reader[cattype]->AddVariable( "sj2", &var[153]);
    reader[cattype]->AddVariable( "sj3", &var[154]);
    //reader[cattype]->AddVariable( "sj4", &var[155]);
    //reader[cattype]->AddVariable( "sj5", &var[142]);
    //reader[cattype]->AddVariable( "sj6", &var[143]);
    //reader[cattype]->AddVariable( "sj7", &var[144]);
    //reader[cattype]->AddVariable( "sj8", &var[145]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "zb", &var[166]);
    reader[cattype]->AddVariable( "hb", &var[167]);
    reader[cattype]->AddVariable( "h1mass", &var[161]);
    reader[cattype]->AddVariable( "h2mass", &var[163]);
    //reader[cattype]->AddVariable( "zmass", &var[162]);
    reader[cattype]->AddVariable( "hmass", &var[147]);
    //reader[cattype]->AddVariable( "zzhchisq1", &var[143]);

    reader[cattype]->BookMVA( "MLP_ttbar_allhad_bbbb_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_allhad_bbbb_new.weights.xml" ); 
    break;
  case 1:    //leptonic(ttbar + zzh)
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    //reader[cattype]->AddVariable( "mW2", &var[3]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    //reader[cattype]->AddVariable( "m4vect", &var[32]);
    //reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    reader[cattype]->AddVariable( "moment2", &var[22]);
    //reader[cattype]->AddVariable( "Cvalue", &var[28]);
    //reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "sj2", &var[139]);  ////for lpos
    reader[cattype]->AddVariable( "sj3", &var[154]);   ////for normal
    //reader[cattype]->AddVariable( "sj4", &var[141]);
    //reader[cattype]->AddVariable( "sj5", &var[142]);
    reader[cattype]->AddVariable( "sj6", &var[157]);   ////for normal
    //reader[cattype]->AddVariable( "sj7", &var[144]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    //reader[cattype]->AddVariable( "cos1", &var[137]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);
    reader[cattype]->AddVariable( "mete", &var[77]);
    reader[cattype]->AddVariable( "h1mass", &var[161]);
    reader[cattype]->AddVariable( "h2mass", &var[163]);
    reader[cattype]->AddVariable( "zmass", &var[162]);
    reader[cattype]->AddVariable( "hmass", &var[147]);
    //reader[cattype]->AddVariable( "z1mass", &var[145]);
    reader[cattype]->BookMVA( "MLP_dilepton_allhad_bbbb_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_dilepton_allhad_bbbb_new.weights.xml" ); 
    break;
  case 2:    //ttbarleptonic
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    //reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    reader[cattype]->AddVariable( "mW1", &var[2]);
    //reader[cattype]->AddVariable( "costhetaZbb", &var[9]);
    //reader[cattype]->AddVariable( "mW2", &var[3]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    //reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    reader[cattype]->AddVariable( "moment2", &var[22]);
    //reader[cattype]->AddVariable( "sphericity", &var[26]);
    //reader[cattype]->AddVariable( "aplanarity", &var[27]);
    //reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "sj8", &var[159]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);
    reader[cattype]->AddVariable( "mete", &var[77]);
    //reader[cattype]->AddVariable( "ctag1", &var[140]);
    //reader[cattype]->AddVariable( "ctag2", &var[141]);
    //reader[cattype]->AddVariable( "ctag3", &var[142]);
    //reader[cattype]->AddVariable( "totctag", &var[137]);
    //reader[cattype]->AddVariable( "sj3", &var[140]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    //reader[cattype]->AddVariable( "zc", &var[131]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "h2mass", &var[163]);
    //reader[cattype]->AddVariable( "zmass", &var[162]);
    reader[cattype]->AddVariable( "hmass", &var[147]);
    //reader[cattype]->AddVariable( "zzhchisq1", &var[143]);
    reader[cattype]->AddVariable( "z1mass", &var[145]);
    reader[cattype]->BookMVA( "MLP_dilepton_allhad_bbcc_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_dilepton_allhad_bbcc_new.weights.xml" ); 
    break;
  case 3:    //ttbarhadronic
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    reader[cattype]->AddVariable( "mZbb", &var[1]);
    //reader[cattype]->AddVariable( "mW1", &var[2]);
    //reader[cattype]->AddVariable( "costhetaZbb", &var[9]);
    //reader[cattype]->AddVariable( "mW2", &var[3]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    //reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    //reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "moment2", &var[22]);
    //reader[cattype]->AddVariable( "sphericity", &var[26]);
    //reader[cattype]->AddVariable( "aplanarity", &var[27]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    //reader[cattype]->AddVariable( "y56", &var[13]);
    //reader[cattype]->AddVariable( "mete", &var[77]);
    reader[cattype]->AddVariable( "sj2", &var[153]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);
    //reader[cattype]->AddVariable( "ctag1", &var[140]);
    //reader[cattype]->AddVariable( "ctag2", &var[141]);
    //reader[cattype]->AddVariable( "ctag3", &var[142]);
    //reader[cattype]->AddVariable( "totctag", &var[137]);
    //reader[cattype]->AddVariable( "sj3", &var[140]);
    //reader[cattype]->AddVariable( "sj8", &var[145]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    //reader[cattype]->AddVariable( "zb", &var[166]);
    reader[cattype]->AddVariable( "hb", &var[167]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "h1mass", &var[162]);
    reader[cattype]->AddVariable( "h2mass", &var[163]);
    reader[cattype]->AddVariable( "zmass", &var[162]);
    reader[cattype]->AddVariable( "hmass", &var[147]);
    //reader[cattype]->AddVariable( "zzhchisq1", &var[143]);
    reader[cattype]->AddVariable( "z1mass", &var[145]);
    reader[cattype]->AddVariable( "z2mass", &var[146]);
    reader[cattype]->AddSpectator("eventtype",&var[198]);

    reader[cattype]->BookMVA( "MLP_ttbar_allhad_bbcc_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_allhad_bbcc_new.weights.xml" ); 
    break;
  case 4:    //ttbarleptonic
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    //reader[cattype]->AddVariable( "mHbb", &var[0]);
    //reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    //reader[cattype]->AddVariable( "costhetaZbb", &var[9]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    //reader[cattype]->AddVariable( "masym1", &var[30]);
    //reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    reader[cattype]->AddVariable( "moment2", &var[22]);
    //reader[cattype]->AddVariable( "Cvalue", &var[28]);
    //reader[cattype]->AddVariable( "Dvalue", &var[29]);
    //reader[cattype]->AddVariable( "sphericity", &var[26]);
    //reader[cattype]->AddVariable( "aplanarity", &var[27]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);
    reader[cattype]->AddVariable( "mete", &var[77]);
    //reader[cattype]->AddVariable( "ctag1", &var[140]);
    //reader[cattype]->AddVariable( "ctag2", &var[141]);
    //reader[cattype]->AddVariable( "ctag3", &var[142]);
    //reader[cattype]->AddVariable( "totctag", &var[143]);
    //reader[cattype]->AddVariable( "sj2", &var[139]);
    //reader[cattype]->AddVariable( "sj3", &var[140]);
    //reader[cattype]->AddVariable( "sj5", &var[142]);
    //reader[cattype]->AddVariable( "sj6", &var[143]);
    //reader[cattype]->AddVariable( "sj8", &var[145]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    //reader[cattype]->AddVariable( "mete", &var[77]);
    //reader[cattype]->AddVariable( "cos1", &var[137]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "h1mass", &var[161]);
    //reader[cattype]->AddVariable( "h2mass", &var[163]);
    //reader[cattype]->AddVariable( "zmass", &var[162]);
    reader[cattype]->AddVariable( "hmass", &var[147]);
    reader[cattype]->AddVariable( "zzhchisq1", &var[143]);
    //reader[cattype]->AddVariable( "z1mass", &var[145]);
    reader[cattype]->BookMVA( "MLP_dilepton_allhad_bbb_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_dilepton_allhad_bbb_new.weights.xml" ); 
    break;
  case 5:    //ttbarhadronic
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    //reader[cattype]->AddVariable( "mHbb", &var[0]);
    //reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    //reader[cattype]->AddVariable( "costhetaZbb", &var[9]);
    //reader[cattype]->AddVariable( "costhetaW1", &var[10]);
    //reader[cattype]->AddVariable( "costhetaW2", &var[11]);
    //reader[cattype]->AddVariable( "mW2", &var[3]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    //reader[cattype]->AddVariable( "m4vect", &var[32]);
    //reader[cattype]->AddVariable( "masym1", &var[30]);
    //reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    //reader[cattype]->AddVariable( "moment1", &var[21]);
    reader[cattype]->AddVariable( "moment2", &var[22]);
    //reader[cattype]->AddVariable( "moment5", &var[25]);
    //reader[cattype]->AddVariable( "sphericity", &var[26]);
    //reader[cattype]->AddVariable( "aplanarity", &var[27]);
    reader[cattype]->AddVariable( "Cvalue", &var[28]);
    reader[cattype]->AddVariable( "Dvalue", &var[29]);
    //reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    //reader[cattype]->AddVariable( "ctag1", &var[140]);
    //reader[cattype]->AddVariable( "ctag2", &var[141]);
    //reader[cattype]->AddVariable( "ctag3", &var[142]);
    //reader[cattype]->AddVariable( "totctag", &var[143]);
    reader[cattype]->AddVariable( "sj1", &var[152]);
    reader[cattype]->AddVariable( "sj2", &var[153]);
    reader[cattype]->AddVariable( "sj3", &var[154]);
    //reader[cattype]->AddVariable( "sj4", &var[141]);
    //reader[cattype]->AddVariable( "sj5", &var[142]);
    //reader[cattype]->AddVariable( "sj6", &var[143]);
    //reader[cattype]->AddVariable( "sj7", &var[144]);
    //reader[cattype]->AddVariable( "sj8", &var[159]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);
    reader[cattype]->AddVariable( "hb", &var[167]);
    //reader[cattype]->AddVariable( "h1mass", &var[161]);
    //reader[cattype]->AddVariable( "h2mass", &var[163]);
    //reader[cattype]->AddVariable( "zmass", &var[162]);
    //reader[cattype]->AddVariable( "hmass", &var[147]);
    //reader[cattype]->AddVariable( "zzhchisq1", &var[143]);
    //reader[cattype]->AddVariable( "z1mass", &var[145]);
    reader[cattype]->BookMVA( "MLP_ttbar_allhad_bbb_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_allhad_bbb_new.weights.xml" ); 
    break;
  case 11:    //lep + jets ttbar cut
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    //reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "sj6", &var[145]);
    reader[cattype]->AddVariable( "totsj", &var[158]);
    reader[cattype]->AddVariable( "topveto1", &var[127]);
    //reader[cattype]->AddVariable( "cos2", &var[139]);
    reader[cattype]->AddVariable( "zb", &var[164]);
    reader[cattype]->AddVariable( "hb", &var[165]);
    reader[cattype]->AddVariable( "h1mass", &var[159]);
    //reader[cattype]->AddVariable( "h2mass", &var[161]);
    //reader[cattype]->AddVariable( "zmass", &var[160]);
    //reader[cattype]->AddVariable( "w1mass2", &var[162]);
    //reader[cattype]->AddVariable( "w2mass2", &var[163]);
    //reader[cattype]->AddVariable( "z2mass", &var[146]);
    //reader[cattype]->AddVariable( "hmass", &var[147]);
   
    reader[cattype]->BookMVA( "MLP_ttbar_lepjets_bbbb_new2", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_lepjets_bbbb_new2.weights.xml" ); 
    break;
  case 12:    //lep + jets dilepton cut
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "sj6", &var[157]);
    reader[cattype]->AddVariable( "totsj", &var[158]);
    //reader[cattype]->AddVariable( "cos2", &var[139]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);
    reader[cattype]->AddVariable( "zb", &var[164]);
    reader[cattype]->AddVariable( "hb", &var[165]);
    reader[cattype]->AddVariable( "h1mass", &var[159]);
    reader[cattype]->AddVariable( "h2mass", &var[161]);
    //reader[cattype]->AddVariable( "zmass", &var[160]);
    reader[cattype]->AddVariable( "w1mass2", &var[162]);
   
    reader[cattype]->BookMVA( "MLP_dilepton_lepjets_bbbb_new2", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_dilepton_lepjets_bbbb_new2.weights.xml" ); 
    break;
  case 13:    //lep + jets
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mZbb", &var[1]);
    reader[cattype]->AddVariable( "mW1", &var[2]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "mtwlnu", &var[79]);
    //reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    //reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    //reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "Cvalue", &var[28]);
    //reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "sj1", &var[140]);
    //reader[cattype]->AddVariable( "sj2", &var[141]);
    //reader[cattype]->AddVariable( "sj5", &var[144]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    //reader[cattype]->AddVariable( "topveto1", &var[127]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);
    reader[cattype]->AddVariable( "mete", &var[77]);
    //reader[cattype]->AddVariable( "h1mass", &var[159]);
    //reader[cattype]->AddVariable( "h2mass", &var[161]);
    //reader[cattype]->AddVariable( "zmass", &var[160]);
    reader[cattype]->AddVariable( "w1mass2", &var[162]);

    reader[cattype]->BookMVA( "MLP_dilepton_lepjets_bbb_new2", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_dilepton_lepjets_bbb_new2.weights.xml" ); 
    break;
  case 14:    //lep + jets hadronic for cc
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mW1", &var[2]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "mtwlnu", &var[79]);
    //reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    //reader[cattype]->AddVariable( "masym2", &var[31]);
    //reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "Cvalue", &var[28]);
    //reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "sj1", &var[140]);
    //reader[cattype]->AddVariable( "sj2", &var[141]);
    //reader[cattype]->AddVariable( "sj5", &var[144]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "topveto1", &var[127]);
    //reader[cattype]->AddVariable( "Evisible", &var[33]);
    reader[cattype]->AddVariable( "mete", &var[77]);
    //reader[cattype]->AddVariable( "h1mass", &var[159]);
    reader[cattype]->AddVariable( "h2mass", &var[161]);
    reader[cattype]->AddVariable( "zmass", &var[160]);
    reader[cattype]->AddVariable( "w1mass2", &var[162]);
  
    reader[cattype]->BookMVA( "MLP_ttbar_lepjets_bbb_new2", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_lepjets_bbb_new2.weights.xml" ); 
    break;
  case 15:    //lep + jets hadronic for cc
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mZbb", &var[1]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    //reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    //reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    reader[cattype]->AddVariable( "Cvalue", &var[28]);
    reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "sj6", &var[145]);
    reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "mete", &var[77]);
    reader[cattype]->AddVariable( "ctag2", &var[141]);
    reader[cattype]->AddVariable( "totctag", &var[143]);
    
    reader[cattype]->BookMVA( "MLP_ttbar_lepjets_cc", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_lepjets_cc.weights.xml" ); 
    break;
  case 21:    //dilepton ttbar
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    reader[cattype]->AddVariable( "mZbb", &var[1]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "topveto1", &var[127]);
    reader[cattype]->AddVariable( "topveto2", &var[128]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);

    reader[cattype]->BookMVA( "MLP_ttbar_dilepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_dilepton.weights.xml" ); 
    break;
  case 22:    //dilepton zzh
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "sj6", &var[145]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->BookMVA( "MLP_zzh_dilepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_zzh_dilepton.weights.xml" ); 
    break;
  case 23:    //dilepton dilepton
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->AddVariable( "Evisible", &var[33]);
    reader[cattype]->BookMVA( "MLP_dilepton_dilepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_dilepton_dilepton.weights.xml" ); 
    break;
  case 24:    //dilepton_blb
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mW1", &var[2]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "Cvalue", &var[28]);
    reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->AddVariable( "sj1", &var[140]);
    reader[cattype]->AddVariable( "sj6", &var[145]);
    reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "cos1", &var[139]);

    reader[cattype]->BookMVA( "MLP_ttbar_dilepton_blb", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_dilepton_blb.weights.xml" ); 
    break;
  case 25:    //dilepton_blb zzh
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "sphericity", &var[26]);
    reader[cattype]->AddVariable( "aplanarity", &var[27]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->AddVariable( "sj3", &var[142]);
    reader[cattype]->AddVariable( "sj4", &var[143]);
    //reader[cattype]->AddVariable( "sj5", &var[144]);
    //reader[cattype]->AddVariable( "sj6", &var[145]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->BookMVA( "MLP_zzh_dilepton_blb", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_zzh_dilepton_blb.weights.xml" ); 
    break;
  case 31:    //trilepton
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mZbb", &var[1]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "topveto1", &var[127]);
    reader[cattype]->AddVariable( "topveto2", &var[128]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->BookMVA( "MLP_ttbar_trilepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_trilepton.weights.xml" ); 
    break;
  case 32:    //trilepton
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "sj4", &var[145]);
    reader[cattype]->AddVariable( "topveto2", &var[128]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->BookMVA( "MLP_zzh_trilepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_zzh_trilepton.weights.xml" ); 
    break;
  case 41:    //nolepton ttbar allhadronic
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    reader[cattype]->AddVariable( "mZbb", &var[1]);
    //reader[cattype]->AddVariable( "costhetaW1", &var[10]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "sphericity", &var[26]);
    //reader[cattype]->AddVariable( "aplanarity", &var[27]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "costhetaZH", &var[81]);
    //reader[cattype]->AddVariable( "sj1", &var[140]);
    //reader[cattype]->AddVariable( "sj2", &var[141]);
    //reader[cattype]->AddVariable( "sj3", &var[142]);
    //reader[cattype]->AddVariable( "sj4", &var[143]);
    //reader[cattype]->AddVariable( "sj5", &var[144]);
    reader[cattype]->AddVariable( "sj6", &var[145]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "mete", &var[77]);
    //reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->AddVariable( "metcos", &var[133]);
    reader[cattype]->BookMVA( "MLP_ttbarhadronic_nolepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbarhadronic_nolepton.weights.xml" ); 
    break;
  case 42:    //nolepton ttbar leptonic
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    //reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    //reader[cattype]->AddVariable( "mZbb", &var[1]);
    //reader[cattype]->AddVariable( "costhetaW1", &var[10]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    //reader[cattype]->AddVariable( "sphericity", &var[26]);
    //reader[cattype]->AddVariable( "aplanarity", &var[27]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "costhetaZH", &var[81]);
    reader[cattype]->AddVariable( "costhetaHWW", &var[83]);
    reader[cattype]->AddVariable( "sj3", &var[142]);
    //reader[cattype]->AddVariable( "sj4", &var[143]);
    //reader[cattype]->AddVariable( "sj1", &var[140]);
    reader[cattype]->AddVariable( "sj6", &var[145]);
    reader[cattype]->AddVariable( "totsj", &var[146]);
    //reader[cattype]->AddVariable( "mete", &var[77]);
    //reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->AddVariable( "metcos", &var[133]);
    reader[cattype]->BookMVA( "MLP_ttbar_nolepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbar_nolepton.weights.xml" ); 
    break;
  }
    return;
}

double expert_ttbar::getMVAvalue(float *par){
  
  for(Int_t i=0;i<200;i++){
    var[i]=(float)par[i];
  }
  
  double value=0.0;
  switch(type){
  case 0:
    value=reader[type]->EvaluateMVA("MLP_ttbar_allhad_bbbb_new");
    break;
  case 1:
    value=reader[type]->EvaluateMVA("MLP_dilepton_allhad_bbbb_new");
    break;
  case 2:
    value=reader[type]->EvaluateMVA("MLP_dilepton_allhad_bbcc_new");
    break;
  case 3:
    value=reader[type]->EvaluateMVA("MLP_ttbar_allhad_bbcc_new");
    break;
  case 4:
    value=reader[type]->EvaluateMVA("MLP_dilepton_allhad_bbb_new");
    break;
  case 5:
    value=reader[type]->EvaluateMVA("MLP_ttbar_allhad_bbb_new");
    break;
  case 11:
    value=reader[type]->EvaluateMVA("MLP_ttbar_lepjets_bbbb_new2");
    break;
  case 12:
    value=reader[type]->EvaluateMVA("MLP_dilepton_lepjets_bbbb_new2");
    break;
  case 13:
    value=reader[type]->EvaluateMVA("MLP_dilepton_lepjets_bbb_new2");
    break;
  case 14:
    value=reader[type]->EvaluateMVA("MLP_ttbar_lepjets_bbb_new2");
    break;
  case 15:
    value=reader[type]->EvaluateMVA("MLP_ttbar_lepjets_cc");
    break;
  case 21:
    value=reader[type]->EvaluateMVA("MLP_ttbar_dilepton");
    break;
  case 22:
    value=reader[type]->EvaluateMVA("MLP_zzh_dilepton");
    break;
  case 23:
    value=reader[type]->EvaluateMVA("MLP_dilepton_dilepton");
    break;
  case 24:
    value=reader[type]->EvaluateMVA("MLP_ttbar_dilepton_blb");
    break;
  case 25:
    value=reader[type]->EvaluateMVA("MLP_zzh_dilepton_blb");
    break;
  case 31:
    value=reader[type]->EvaluateMVA("MLP_ttbar_trilepton");
    break;
  case 32:
    value=reader[type]->EvaluateMVA("MLP_zzh_trilepton");
    break;
  case 41:
    value=reader[type]->EvaluateMVA("MLP_ttbarhadronic_nolepton");
    break;
  case 42:
    value=reader[type]->EvaluateMVA("MLP_ttbar_nolepton");
    break;
  }
  
  return value;
}

expert_ttbar::~expert_ttbar(){
  delete reader[type];
  return;
}

expert_ttbb::expert_ttbb(int cattype){
  
  //cattype indicates the separated category types
  //0: lep+jets, 4btag case
  //1: lep+jets, 3btag case
  //2: lep+jets, 2btag+ctag case
  //3: lep+jets, 2btag case
  //4-10: reserved
  //11: all hadronic 4btag case
  //12: all hadronic 3btag(?) case
  //13: all hadronic 2btag+ctag(?) case
  //14-20: reserved
  //21: dilepton, 2btag case(reserved)
  type=cattype;
  
  switch(cattype){
  case 0:    //all hadronic
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    //reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "moment2", &var[22]);
    reader[cattype]->AddVariable( "Cvalue", &var[28]);
    reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    //reader[cattype]->AddVariable( "y56", &var[13]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "sj1", &var[152]);   ////for lpos
    reader[cattype]->AddVariable( "sj2", &var[153]);   ////forlpos
    //reader[cattype]->AddVariable( "sj3", &var[140]);
    //reader[cattype]->AddVariable( "sj4", &var[141]);
    //reader[cattype]->AddVariable( "sj5", &var[156]);  ////for normal
    //reader[cattype]->AddVariable( "sj6", &var[143]);
    //reader[cattype]->AddVariable( "sj7", &var[144]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "h1mass", &var[161]);
    reader[cattype]->AddVariable( "h2mass", &var[163]);
    reader[cattype]->AddVariable( "zmass", &var[162]);
    reader[cattype]->AddVariable( "hmass", &var[147]);
    //reader[cattype]->AddVariable( "zzhchisq1", &var[143]);
    reader[cattype]->BookMVA( "MLP_ttbb_allhad_bbbb_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbb_allhad_bbbb_new.weights.xml" ); 
    break;
  case 1:    //all hadronic for bbcc
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0]);
    reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    //reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    //reader[cattype]->AddVariable( "m6vect", &var[78]);
    //reader[cattype]->AddVariable( "masym2", &var[31]);
    //reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "moment2", &var[22]);
    reader[cattype]->AddVariable( "Cvalue", &var[28]);
    reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "sj3", &var[140]);
    reader[cattype]->AddVariable( "sj6", &var[157]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    //reader[cattype]->AddVariable( "zb", &var[166]);
    reader[cattype]->AddVariable( "h1mass", &var[161]);
    reader[cattype]->AddVariable( "h2mass", &var[163]);
    reader[cattype]->AddVariable( "z1mass", &var[145]);
    //reader[cattype]->AddVariable( "hmass", &var[147]);

    reader[cattype]->BookMVA( "MLP_ttbb_allhad_bbcc_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbb_allhad_bbcc_new.weights.xml" ); 
    break;
  case 2:    //all hadronic for bbb
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    //reader[cattype]->AddVariable( "mHbb", &var[0]);
    //reader[cattype]->AddVariable( "costhetaHbb", &var[8]);
    reader[cattype]->AddVariable( "mW2", &var[3]);
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    //reader[cattype]->AddVariable( "m6vect", &var[78]);
    //reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    //reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "moment2", &var[22]);
    reader[cattype]->AddVariable( "Cvalue", &var[28]);
    //reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y78", &var[15]);
    //reader[cattype]->AddVariable( "y56", &var[13]);
    //reader[cattype]->AddVariable( "Evisible", &var[33]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    //reader[cattype]->AddVariable( "sj4", &var[141]);
    //reader[cattype]->AddVariable( "sj6", &var[143]);
    //reader[cattype]->AddVariable( "sj5", &var[142]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "h1mass", &var[161]);
    reader[cattype]->AddVariable( "h2mass", &var[163]);
    reader[cattype]->AddVariable( "zmass", &var[162]);
    //reader[cattype]->AddVariable( "hmass", &var[147]);
    //reader[cattype]->AddVariable( "z1mass", &var[145]);
    reader[cattype]->BookMVA( "MLP_ttbb_allhad_bbb_new", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbb_allhad_bbb_new.weights.xml" ); 
    break;
  case 11:    //lep + jets
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mW1", &var[2] );
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    //reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    reader[cattype]->AddVariable( "Cvalue", &var[28]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "sj6", &var[157]);
    //reader[cattype]->AddVariable( "totsj", &var[158]);
    reader[cattype]->AddVariable( "topveto1", &var[127]);
    //reader[cattype]->AddVariable( "h1mass", &var[159]);
    reader[cattype]->AddVariable( "h2mass", &var[161]);
    //reader[cattype]->AddVariable( "zmass", &var[160]);
    reader[cattype]->AddVariable( "w1mass2", &var[162]);
    //reader[cattype]->AddVariable( "w2mass2", &var[163]);

    reader[cattype]->BookMVA( "MLP_ttbb_lepjets_bbbb_new2", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbb_lepjets_bbbb_new2.weights.xml" ); 
    break;
  case 12:    //lep + jets
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mW1", &var[2] );
    //reader[cattype]->AddVariable( "costhetaW1", &var[2] );
    reader[cattype]->AddVariable( "mtwlnu", &var[2] );
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    //reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "Cvalue", &var[28]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    //reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "sj1", &var[152]);
    //reader[cattype]->AddVariable( "sj2", &var[153]);
    //reader[cattype]->AddVariable( "totsj", &var[158]);
    reader[cattype]->AddVariable( "topveto1", &var[127]);
    reader[cattype]->AddVariable( "h1mass", &var[159]);
    reader[cattype]->AddVariable( "h2mass", &var[161]);
    reader[cattype]->AddVariable( "w2mass2", &var[163]);
    reader[cattype]->AddVariable( "sigchisq2", &var[134]);
    reader[cattype]->AddVariable( "w1mass", &var[161]);
    //reader[cattype]->AddVariable( "w2mass2", &var[163]);
    reader[cattype]->BookMVA( "MLP_ttbb_lepjets_bbb_new2", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbb_lepjets_bbb_new2.weights.xml" ); 
    break;
  case 21:    //dilepton
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mZbb", &var[1] );
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "topveto1", &var[127]);
    reader[cattype]->AddVariable( "topveto2", &var[128]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);

    reader[cattype]->BookMVA( "MLP_ttbb_dilepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbb_dilepton.weights.xml" ); 
    break;
  case 22:    //dilepton_blb
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0] );
    reader[cattype]->AddVariable( "costhetaHbb", &var[8] );
    reader[cattype]->AddVariable( "mW2", &var[3] );
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment0", &var[20]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    reader[cattype]->AddVariable( "colzmass", &var[197]);
    //reader[cattype]->AddVariable( "sj1", &var[140]);
    reader[cattype]->AddVariable( "sj2", &var[141]);
    reader[cattype]->AddVariable( "sj3", &var[142]);
    //reader[cattype]->AddVariable( "sj4", &var[143]);
    //reader[cattype]->AddVariable( "totsj", &var[146]);
    reader[cattype]->AddVariable( "cos1", &var[139]);

    reader[cattype]->BookMVA( "MLP_ttbb_dilepton_blb", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbb_dilepton_blb.weights.xml" ); 
    break;
  case 31:    //trilepton
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0] );
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "topveto2", &var[128]);

    reader[cattype]->BookMVA( "MLP_ttz_trilepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttz_trilepton.weights.xml" ); 
    break;
  case 41:    //nolepton
    reader[cattype]=new TMVA::Reader( "!Color:!Silent" );
    reader[cattype]->AddVariable( "mHbb", &var[0] );
    reader[cattype]->AddVariable( "costhetaHbb", &var[8] );
    //reader[cattype]->AddVariable( "mZbb", &var[1] );
    reader[cattype]->AddVariable( "mW2", &var[3] );
    reader[cattype]->AddVariable( "m4vect", &var[32]);
    reader[cattype]->AddVariable( "masym1", &var[30]);
    reader[cattype]->AddVariable( "masym2", &var[31]);
    reader[cattype]->AddVariable( "moment1", &var[21]);
    //reader[cattype]->AddVariable( "Cvalue", &var[28]);
    //reader[cattype]->AddVariable( "Dvalue", &var[29]);
    reader[cattype]->AddVariable( "chisq", &var[12]);
    reader[cattype]->AddVariable( "y56", &var[13]);
    //reader[cattype]->AddVariable( "colzmass", &var[197]);
    reader[cattype]->AddVariable( "costhetaWW", &var[80]);
    reader[cattype]->AddVariable( "sj3", &var[142]);
    reader[cattype]->AddVariable( "sj4", &var[143]);
    //reader[cattype]->AddVariable( "sj5", &var[144]);
    reader[cattype]->AddVariable( "mete", &var[77]);
    reader[cattype]->BookMVA( "MLP_ttbb_nolepton", "/home/ilc/kurata/myAna/myAna/Ana/test/weights/TMVAClassification_MLP_ttbb_nolepton.weights.xml" ); 
    break;
   }
    return;
}

double expert_ttbb::getMVAvalue(float *par){
  
  for(Int_t i=0;i<200;i++){
    var[i]=(float)par[i];
  }
  
  double value=0.0;
  switch(type){
  case 0:
    value=reader[type]->EvaluateMVA("MLP_ttbb_allhad_bbbb_new");
    break;
  case 1:
    value=reader[type]->EvaluateMVA("MLP_ttbb_allhad_bbcc_new");
    break;
  case 2:
    value=reader[type]->EvaluateMVA("MLP_ttbb_allhad_bbb_new");
    break;
  case 11:
    value=reader[type]->EvaluateMVA("MLP_ttbb_lepjets_bbbb_new2");
    break;
  case 12:
    value=reader[type]->EvaluateMVA("MLP_ttbb_lepjets_bbb_new2");
    break;
  case 21:
    value=reader[type]->EvaluateMVA("MLP_ttbb_dilepton");
    break;
  case 22:
    value=reader[type]->EvaluateMVA("MLP_ttbb_dilepton_blb");
    break;
  case 31:
    value=reader[type]->EvaluateMVA("MLP_ttz_trilepton");
    break;
  case 41:
    value=reader[type]->EvaluateMVA("MLP_ttbb_nolepton");
    break;
  }
  
  return value;
}

expert_ttbb::~expert_ttbb(){
  delete reader[type];
  return;
}

softjetFinder::softjetFinder(){
  //clusterNN
  disc_cluster=new TMVA::Reader( "!Color:!Silent" );
  disc_cluster->AddVariable( "trk1", &fvar[0]);
  disc_cluster->AddVariable( "trk2", &fvar[1]);
  disc_cluster->AddVariable( "trk3", &fvar[2]);
  disc_cluster->AddVariable( "trk4", &fvar[3]);
  disc_cluster->AddVariable( "trk5", &fvar[4]);
  disc_cluster->AddVariable( "trk6", &fvar[5]);
  disc_cluster->AddVariable( "trk7", &fvar[6]);
  disc_cluster->AddVariable( "trk8", &fvar[7]);
  disc_cluster->AddVariable( "trk9", &fvar[8]);
  disc_cluster->AddVariable( "trk10", &fvar[9]);
  disc_cluster->AddVariable( "trk11", &fvar[10]);
  disc_cluster->AddVariable( "trk12", &fvar[11]);
  disc_cluster->AddVariable( "trk13", &fvar[12]);
  disc_cluster->AddVariable( "trk14", &fvar[13]);
  disc_cluster->AddVariable( "trk15", &fvar[14]);
  disc_cluster->AddVariable( "trk16", &fvar[15]);
  disc_cluster->AddVariable( "trk17", &fvar[16]);
  disc_cluster->AddVariable( "trk18", &fvar[17]);
  disc_cluster->AddVariable( "trk19", &fvar[18]);
  disc_cluster->AddVariable( "trk20", &fvar[19]);
  // disc_cluster->AddVariable( "trk21", &fvar[20]);
  // disc_cluster->AddVariable( "trk22", &fvar[21]);
  // disc_cluster->AddVariable( "trk23", &fvar[22]);
  // disc_cluster->AddVariable( "trk24", &fvar[23]);
  // disc_cluster->AddVariable( "trk25", &fvar[24]);
  // disc_cluster->AddVariable( "trk26", &fvar[25]);
  // disc_cluster->AddVariable( "trk27", &fvar[26]);
  // disc_cluster->AddVariable( "trk28", &fvar[27]);
  // disc_cluster->AddVariable( "trk29", &fvar[28]);
  // disc_cluster->AddVariable( "trk30", &fvar[29]);
  disc_cluster->BookMVA( "MLP_track_new", "test/weights/TMVAClassification_MLP_track_new.weights.xml" ); 
 
  //trackNN 
  disc_track=new TMVA::Reader( "!Color:!Silent" );
  disc_track->AddVariable( "trk1", &fvar[0]);
  disc_track->AddVariable( "trk2", &fvar[1]);
  disc_track->AddVariable( "trk3", &fvar[2]);
  disc_track->AddVariable( "trk4", &fvar[3]);
  disc_track->AddVariable( "trk5", &fvar[4]);
  disc_track->AddVariable( "trk6", &fvar[5]);
  disc_track->AddVariable( "trk7", &fvar[6]);
  disc_track->AddVariable( "trk8", &fvar[7]);
  disc_track->AddVariable( "trk9", &fvar[8]);
  disc_track->AddVariable( "trk10", &fvar[9]);
  disc_track->AddVariable( "trk11", &fvar[10]);
  disc_track->AddVariable( "trk12", &fvar[11]);
  disc_track->AddVariable( "trk13", &fvar[12]);
  disc_track->AddVariable( "trk14", &fvar[13]);
  disc_track->AddVariable( "trk15", &fvar[14]);
  disc_track->AddVariable( "trk16", &fvar[15]);
  disc_track->AddVariable( "trk17", &fvar[16]);
  disc_track->AddVariable( "trk18", &fvar[17]);
  disc_track->AddVariable( "trk19", &fvar[18]);
  disc_track->AddVariable( "trk20", &fvar[19]);
  // disc_track->AddVariable( "trk21", &fvar[50]);
  // disc_track->AddVariable( "trk22", &fvar[51]);
  // disc_track->AddVariable( "trk23", &fvar[52]);
  // disc_track->AddVariable( "trk24", &fvar[53]);
  // disc_track->AddVariable( "trk25", &fvar[54]);
  // disc_track->AddVariable( "trk26", &fvar[55]);
  // disc_track->AddVariable( "trk27", &fvar[56]);
  // disc_track->AddVariable( "trk28", &fvar[57]);
  // disc_track->AddVariable( "trk29", &fvar[58]);
  // disc_track->AddVariable( "trk30", &fvar[59]);
  disc_track->BookMVA( "MLP_track", "weights/TMVAClassification_MLP_track.weights.xml" ); 

  disc_cluster2=new TMVA::Reader( "!Color:!Silent" );
  disc_cluster2->AddVariable( "trk1", &fvar[0]);
  disc_cluster2->AddVariable( "trk2", &fvar[1]);
  disc_cluster2->AddVariable( "trk3", &fvar[2]);
  disc_cluster2->AddVariable( "trk4", &fvar[3]);
  disc_cluster2->AddVariable( "trk5", &fvar[4]);
  disc_cluster2->AddVariable( "trk6", &fvar[5]);
  disc_cluster2->AddVariable( "trk7", &fvar[6]);
  disc_cluster2->AddVariable( "trk8", &fvar[7]);
  disc_cluster2->AddVariable( "trk9", &fvar[8]);
  disc_cluster2->AddVariable( "trk10", &fvar[9]);
  disc_cluster2->AddVariable( "trk11", &fvar[10]);
  disc_cluster2->AddVariable( "trk12", &fvar[11]);
  disc_cluster2->AddVariable( "trk13", &fvar[12]);
  disc_cluster2->AddVariable( "trk14", &fvar[13]);
  disc_cluster2->AddVariable( "trk15", &fvar[14]);
  disc_cluster2->AddVariable( "trk16", &fvar[15]);
  disc_cluster2->AddVariable( "trk17", &fvar[16]);
  disc_cluster2->AddVariable( "trk18", &fvar[17]);
  disc_cluster2->AddVariable( "trk19", &fvar[18]);
  disc_cluster2->AddVariable( "trk20", &fvar[19]);
  // disc_cluster2->AddVariable( "trk21", &fvar[50]);
  // disc_cluster2->AddVariable( "trk22", &fvar[51]);
  // disc_cluster2->AddVariable( "trk23", &fvar[52]);
  // disc_cluster2->AddVariable( "trk24", &fvar[53]);
  // disc_cluster2->AddVariable( "trk25", &fvar[54]);
  // disc_cluster2->AddVariable( "trk26", &fvar[55]);
  // disc_cluster2->AddVariable( "trk27", &fvar[56]);
  // disc_cluster2->AddVariable( "trk28", &fvar[57]);
  // disc_cluster2->AddVariable( "trk29", &fvar[58]);
  // disc_cluster2->AddVariable( "trk30", &fvar[59]);
  disc_cluster2->BookMVA( "MLP_track2_new", "test/weights/TMVAClassification_MLP_track2_new.weights.xml" ); 

  disc_track2=new TMVA::Reader( "!Color:!Silent" );
  disc_track2->AddVariable( "trk1", &fvar[0]);
  disc_track2->AddVariable( "trk2", &fvar[1]);
  disc_track2->AddVariable( "trk3", &fvar[2]);
  disc_track2->AddVariable( "trk4", &fvar[3]);
  disc_track2->AddVariable( "trk5", &fvar[4]);
  disc_track2->AddVariable( "trk6", &fvar[5]);
  disc_track2->AddVariable( "trk7", &fvar[6]);
  disc_track2->AddVariable( "trk8", &fvar[7]);
  disc_track2->AddVariable( "trk9", &fvar[8]);
  disc_track2->AddVariable( "trk10", &fvar[9]);
  disc_track2->AddVariable( "trk11", &fvar[10]);
  disc_track2->AddVariable( "trk12", &fvar[11]);
  disc_track2->AddVariable( "trk13", &fvar[12]);
  disc_track2->AddVariable( "trk14", &fvar[13]);
  disc_track2->AddVariable( "trk15", &fvar[14]);
  disc_track2->AddVariable( "trk16", &fvar[15]);
  disc_track2->AddVariable( "trk17", &fvar[16]);
  disc_track2->AddVariable( "trk18", &fvar[17]);
  disc_track2->AddVariable( "trk19", &fvar[18]);
  disc_track2->AddVariable( "trk20", &fvar[19]);
  // disc_track2->AddVariable( "trk21", &fvar[50]);
  // disc_track2->AddVariable( "trk22", &fvar[51]);
  // disc_track2->AddVariable( "trk23", &fvar[52]);
  // disc_track2->AddVariable( "trk24", &fvar[53]);
  // disc_track2->AddVariable( "trk25", &fvar[54]);
  // disc_track2->AddVariable( "trk26", &fvar[55]);
  // disc_track2->AddVariable( "trk27", &fvar[56]);
  // disc_track2->AddVariable( "trk28", &fvar[57]);
  // disc_track2->AddVariable( "trk29", &fvar[58]);
  // disc_track2->AddVariable( "trk30", &fvar[59]);
  disc_track2->BookMVA( "MLP_track2", "weights/TMVAClassification_MLP_track2.weights.xml" ); 

  //softjet finder for isolated gluon
  disc_softjet = new TMVA::Reader( "!Color:!Silent" );
  disc_softjet->AddVariable( "jeteta", &sjvar[0]);
  disc_softjet->AddVariable( "hadem", &sjvar[1]);
  disc_softjet->AddVariable( "trkmaxptrel/jete", &sjvar[2]);
  disc_softjet->AddVariable( "neumaxptrel/jete", &sjvar[3]);
  disc_softjet->AddVariable( "trkptrelsum/jete", &sjvar[4]);
  disc_softjet->AddVariable( "neuptrelsum/jete", &sjvar[5]);  
  disc_softjet->AddVariable( "ntrack+nneu", &sjvar[6]);  
  disc_softjet->AddVariable( "cluNN", &sjvar[7]);
  disc_softjet->AddVariable( "Bjet", &sjvar[8] );
  disc_softjet->AddVariable( "Aa/jete", &sjvar[9] );
  //disc_softjet->AddVariable( "twomom", &sjvar[10] );
  disc_softjet->AddVariable( "quartic", &sjvar[11] );
  disc_softjet->AddVariable( "eccentricity", &sjvar[12] );

  disc_softjet->AddSpectator( "jete", &sjvar[14] );
  disc_softjet->BookMVA( "MLP_softjet_new", "test/weights/TMVAClassification_MLP_softjet_new.weights.xml" );

  //softjet finder for collinear jet
  disc_softjet2 = new TMVA::Reader( "!Color:!Silent" );    
  //disc_softjet2->AddVariable( "jeteta", &sjvar[0]);
  disc_softjet2->AddVariable( "hadem", &sjvar[1]);
  disc_softjet2->AddVariable( "trkmaxptrel/jete", &sjvar[2]);
  disc_softjet2->AddVariable( "neumaxptrel/jete", &sjvar[3]);
  disc_softjet2->AddVariable( "trkptrelsum/jete", &sjvar[4]);
  disc_softjet2->AddVariable( "neuptrelsum/jete", &sjvar[5]);  
  disc_softjet2->AddVariable( "ntrack+nneu", &sjvar[6]);  
  disc_softjet2->AddVariable( "cluNN", &sjvar[7]);
  disc_softjet2->AddVariable( "Bjet", &sjvar[8] );
  disc_softjet2->AddVariable( "Aa/jete", &sjvar[9] );
  //disc_softjet2->AddVariable( "twomom", &sjvar[10] );
  //disc_softjet2->AddVariable( "quartic", &sjvar[11] );
  disc_softjet2->AddVariable( "eccentricity", &sjvar[12] );
  disc_softjet2->AddVariable( "planarflow", &sjvar[13] );
  
  disc_softjet2->AddSpectator( "jete", &sjvar[14] );
  disc_softjet2->BookMVA( "MLP_softjet2_new", "test/weights/TMVAClassification_MLP_softjet2_new.weights.xml" );
  
  // string histtxt;
  // for(Int_t i=0;i<2;i++){
  //   histtxt="tmphist" + itos(i+1);
  //   if(i==0) tmphist[i]=new TH1F(histtxt.c_str(),"",50,0.75,1.0);
  //   if(i==1) tmphist[i]=new TH1F(histtxt.c_str(),"",50,0.75,1.0);
  // }
  pullangle2[0]=0.0;
  pullangle2[1]=0.0;
  pullangle[0]=0.0;
  pullangle[1]=0.0;
  
  return;
}

softjetFinder::~softjetFinder(){
  //delete disc_cluster;
  //delete disc_track;
  //delete disc_softjet;

  return;
}

void softjetFinder::get_nearjetPull(int jetid, jetdata data){  //please call it before get_MVAvalue
  get_Pullvalue((int)data.nearjetid[jetid], jetid, data, 2);
  return;
}

void softjetFinder::get_jetPull(int jetid, jetdata data){  //please call it before get_MVAvalue
  get_Pullvalue(jetid, (int)data.nearjetid[jetid], data, 1);
  return;
}

void softjetFinder::get_Pullvalue(int jetid1, int jetid2, jetdata data, int type){
  TLorentzVector jt(data.jetpx[jetid1],
		    data.jetpy[jetid1],
		    data.jetpz[jetid1],
		    data.jete[jetid1]);

  TLorentzVector jt2(data.jetpx[jetid2],
		     data.jetpy[jetid2],
		     data.jetpz[jetid2],
		     data.jete[jetid2]);
  double pvars[2]={0.0,0.0};

  //Pull
  double tt[2]={0.0,0.0};
  TLorentzVector tr,tr2;
  double delphi,delphi2,mag;
  if(data.tr_npart>50) data.tr_npart=50;
  if(data.nu_npart>50) data.nu_npart=50;

  for(int n=0;n<data.tr_npart;n++){
    if(data.tr_e[n]<sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0)))
      data.tr_e[n]=sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0));
    tr.SetPxPyPzE(data.tr_px[n],
		  data.tr_py[n],
		  data.tr_pz[n],
		  data.tr_e[n]);
    delphi=tr.DeltaPhi(jt);
    mag=sqrt(pow(tr.Rapidity()-jt.Rapidity(),2.0)+pow(delphi,2.0));
    //cout << "check1 " << tr.Rapidity() << endl;
    tt[0]+=tr.Pt()*mag*(tr.Rapidity()-jt.Rapidity());
    tt[1]+=tr.Pt()*mag*delphi;  //tr.DeltaPhi(jt);
  }
  for(int n=0;n<data.nu_npart;n++){
    if(data.nu_e[n]<sqrt(pow(data.nu_px[n],2.0)+pow(data.nu_py[n],2.0)+pow(data.nu_pz[n],2.0)))
      data.nu_e[n]=sqrt(pow(data.nu_px[n],2.0)+pow(data.nu_py[n],2.0)+pow(data.nu_pz[n],2.0));
    tr2.SetPxPyPzE(data.nu_px[n],
		   data.nu_py[n],
		   data.nu_pz[n],
		   data.nu_e[n]);
    delphi2=tr2.DeltaPhi(jt);
    mag=sqrt(pow(tr2.Rapidity()-jt.Rapidity(),2.0)+pow(delphi2,2.0));

    //cout << "check2 " << tr2.Rapidity() << endl;
    tt[0]+=tr2.Pt()*mag*(tr2.Rapidity()-jt.Rapidity());
    tt[1]+=tr2.Pt()*mag*delphi2;  //tr.DeltaPhi(jt);
  }
  tt[0]=tt[0]/sqrt(pow(data.jetpx[jetid1],2.0)+pow(data.jetpy[jetid1],2.0));
  tt[1]=tt[1]/sqrt(pow(data.jetpx[jetid1],2.0)+pow(data.jetpy[jetid1],2.0));
  
  TVector2 ttt1(tt[0],tt[1]);
  TVector2 ttt2(jt2.Rapidity()-jt.Rapidity(),jt2.DeltaPhi(jt));
  pvars[0]=ttt1.DeltaPhi(ttt2);
  pvars[1]=(ttt1.Px()*ttt2.Px()+ttt1.Py()*ttt2.Py())
    /sqrt(ttt2.Px()*ttt2.Px()+ttt2.Py()*ttt2.Py());

  switch(type){
  case 1:
    pullangle[1]=pvars[1];
    pullangle[0]=pvars[0];
    break;
  case 2:
    pullangle2[1]=pvars[1];
    pullangle2[0]=pvars[0];
    //cout << pullangle2[0] << " " << pullangle2[1] << endl;
    break;
  }

  return;
}

void softjetFinder::get_Bjet(int jetid, jetdata data){
  //cout << "come Bjet" << endl;

  TLorentzVector jt(data.jetpx[jetid],
		    data.jetpy[jetid],
		    data.jetpz[jetid],
		    data.jete[jetid]);

  //jet broadening
  Bjet=0.0;
  double totp=0.0;
  TVector3 nn=jt.Vect().Unit();
  TVector3 tr,tr2,cr,cr2;
  if(data.tr_npart>50) data.tr_npart=50;
  if(data.nu_npart>50) data.nu_npart=50;

  for(int n=0;n<data.tr_npart;n++){
    if(data.tr_e[n]<0.0) continue;
    tr.SetXYZ(data.tr_px[n],
	      data.tr_py[n],
	      data.tr_pz[n]);
    cr=tr.Cross(nn);
    //cout << "cr " << cr.Mag()<< endl;
    Bjet+=cr.Mag();
    totp+=tr.Mag();
  }
  //cout << "Bjet:" << Bjet << " " << totp << endl;
  
  for(int n=0;n<data.nu_npart;n++){
    if(data.nu_e[n]<0.0) continue;
    tr2.SetXYZ(data.nu_px[n],
	       data.nu_py[n],
	       data.nu_pz[n]);
    cr2=tr2.Cross(nn);
    Bjet+=cr2.Mag();
    totp+=tr2.Mag();
  }
  //cout << "Bjet2:" << Bjet << " " << totp << endl;

  Bjet=Bjet/totp;

  //cout << "Bjet " << Bjet << endl;
  return;
}

void softjetFinder::get_Aa(int jetid, jetdata data){
  //cout << "come Bjet" << endl;

  TLorentzVector jt(data.jetpx[jetid],
		    data.jetpy[jetid],
		    data.jetpz[jetid],
		    data.jete[jetid]);

  //jet angularities
  Aa=0.0;
  double aaa=0.99;
  Double_t RR=1.5;   //edge of jet cone
  TLorentzVector tr;
  double tmptheta,ftheta;
  if(data.tr_npart>50) data.tr_npart=50;
  if(data.nu_npart>50) data.nu_npart=50;

    for(int n=0;n<data.tr_npart;n++){
    if(data.tr_e[n]<0.0) continue;
    if(data.tr_e[n]<sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0)))
      data.tr_e[n]=sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0));
    tr.SetPxPyPzE(data.tr_px[n],
		  data.tr_py[n],
		  data.tr_pz[n],
		  data.tr_e[n]);
    if(jt.DeltaR(tr)>RR) continue;
    tmptheta=TMath::Pi()*jt.DeltaR(tr)/(2*RR);
    ftheta=pow(TMath::Sin(tmptheta),aaa)*pow(1-TMath::Cos(tmptheta),1-aaa);
    //cout << "check " << tmptheta << " " << ftheta << endl;
    Aa+=tr.E()*ftheta;
  }
  for(int n=0;n<data.nu_npart;n++){
    if(data.nu_e[n]<0.0) continue;
    if(data.nu_e[n]<sqrt(pow(data.nu_px[n],2.0)+pow(data.nu_py[n],2.0)+pow(data.nu_pz[n],2.0)))
      data.nu_e[n]=sqrt(pow(data.nu_px[n],2.0)+pow(data.nu_py[n],2.0)+pow(data.nu_pz[n],2.0));
    tr.SetPxPyPzE(data.nu_px[n],
		  data.nu_py[n],
		  data.nu_pz[n],
		  data.nu_e[n]);
    if(jt.DeltaR(tr)>RR) continue;
    tmptheta=TMath::Pi()*jt.DeltaR(tr)/(2*RR);
    ftheta=pow(TMath::Sin(tmptheta),aaa)*pow(1-TMath::Cos(tmptheta),1-aaa);
    //cout << "check " << tmptheta << " " << ftheta << endl;
    Aa+=tr.E()*ftheta;
  }
  
  //cout << "Aa " << Aa << endl;
  return;
}

void softjetFinder::get_twomom(int jetid, jetdata data){
  Tbeta=0.0;

  TLorentzVector tra,trb;
  if(data.tr_npart>50) data.tr_npart=50;
  if(data.nu_npart>50) data.nu_npart=50;

  for(int n=0;n<data.tr_npart+data.nu_npart;n++){
    if(n<data.tr_npart){
      if(data.tr_e[n]<sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0)))
	data.tr_e[n]=sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0));
      
      tra.SetPxPyPzE(data.tr_px[n],
		     data.tr_py[n],
		     data.tr_pz[n],
		     data.tr_e[n]);
    }else{
      if(data.nu_e[n-data.tr_npart]<sqrt(pow(data.nu_px[n-data.tr_npart],2.0)+pow(data.nu_py[n-data.tr_npart],2.0)+pow(data.nu_pz[n-data.tr_npart],2.0)))
	data.nu_e[n-data.tr_npart]=sqrt(pow(data.nu_px[n-data.tr_npart],2.0)+pow(data.nu_py[n-data.tr_npart],2.0)+pow(data.nu_pz[n-data.tr_npart],2.0));
      tra.SetPxPyPzE(data.nu_px[n-data.tr_npart],
		     data.nu_py[n-data.tr_npart],
		     data.nu_pz[n-data.tr_npart],
		     data.nu_e[n-data.tr_npart]);
    }
    
    for(int m=0;m<data.tr_npart+data.nu_npart;m++){
      if(m<data.tr_npart){
	if(data.tr_e[m]<sqrt(pow(data.tr_px[m],2.0)+pow(data.tr_py[m],2.0)+pow(data.tr_pz[m],2.0)))
	  data.tr_e[m]=sqrt(pow(data.tr_px[m],2.0)+pow(data.tr_py[m],2.0)+pow(data.tr_pz[m],2.0));
	
	trb.SetPxPyPzE(data.tr_px[m],
		       data.tr_py[m],
		       data.tr_pz[m],
		       data.tr_e[m]);
      }else{
	if(data.nu_e[m-data.tr_npart]<sqrt(pow(data.nu_px[m-data.tr_npart],2.0)+pow(data.nu_py[m-data.tr_npart],2.0)+pow(data.nu_pz[m-data.tr_npart],2.0)))
	  data.nu_e[m-data.tr_npart]=sqrt(pow(data.nu_px[m-data.tr_npart],2.0)+pow(data.nu_py[m-data.tr_npart],2.0)+pow(data.nu_pz[m-data.tr_npart],2.0));
	
	trb.SetPxPyPzE(data.nu_px[m-data.tr_npart],
		       data.nu_py[m-data.tr_npart],
		       data.nu_pz[m-data.tr_npart],
		       data.nu_e[m-data.tr_npart]);
      } 
      //Tbeta+=tra.Pt()*trb.Pt()*pow(tra.DeltaR(trb),2.0);
      Tbeta+=tra.P()*trb.P()*pow(TMath::Cos(tra.Angle(trb.Vect())),8.0);
    }
  }
  Tbeta=Tbeta/(pow(data.jetpx[jetid],2.0)+pow(data.jetpy[jetid],2.0)+pow(data.jetpz[jetid],2.0));

  //cout << "Tbeta " << Tbeta << endl;
  return;
}

void softjetFinder::get_geomom(int jetid, jetdata data){
  //cout << "come Bjet" << endl;

  TLorentzVector jt(data.jetpx[jetid],
		    data.jetpy[jetid],
		    data.jetpz[jetid],
		    data.jete[jetid]);
  
  //2d geometric moments
  TMatrixD cc(2,2);
  cc[0][0]=0.0;
  cc[0][1]=0.0;
  cc[1][0]=0.0;
  cc[1][1]=0.0;
  
  //making matrix
  TLorentzVector tra;
  double mag;
  if(data.tr_npart>50) data.tr_npart=50;
  if(data.nu_npart>50) data.nu_npart=50;

  for(int n=0;n<data.tr_npart+data.nu_npart;n++){
    if(n<data.tr_npart){
      if(data.tr_e[n]<sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0)))
	data.tr_e[n]=sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0));
      
      tra.SetPxPyPzE(data.tr_px[n],
		     data.tr_py[n],
		     data.tr_pz[n],
		     data.tr_e[n]);
    }else{
      if(data.nu_e[n-data.tr_npart]<sqrt(pow(data.nu_px[n-data.tr_npart],2.0)+pow(data.nu_py[n-data.tr_npart],2.0)+pow(data.nu_pz[n-data.tr_npart],2.0)))
	data.nu_e[n-data.tr_npart]=sqrt(pow(data.nu_px[n-data.tr_npart],2.0)+pow(data.nu_py[n-data.tr_npart],2.0)+pow(data.nu_pz[n-data.tr_npart],2.0));
    
      tra.SetPxPyPzE(data.nu_px[n-data.tr_npart],
		     data.nu_py[n-data.tr_npart],
		     data.nu_pz[n-data.tr_npart],
		     data.nu_e[n-data.tr_npart]);
    }
    mag=sqrt(pow(tra.Rapidity()-jt.Rapidity(),2.0)+pow(tra.Phi()-jt.Phi(),2.0));
    cc[0][0]+=tra.Pt()*mag*pow(tra.Rapidity()-jt.Rapidity(),2.0);
    cc[0][1]+=tra.Pt()*mag*(tra.Phi()-jt.Phi())*(tra.Rapidity()-jt.Rapidity());
    cc[1][0]+=tra.Pt()*mag*(tra.Phi()-jt.Phi())*(tra.Rapidity()-jt.Rapidity());
    cc[1][1]+=tra.Pt()*mag*pow(tra.Phi()-jt.Phi(),2.0);
  }
  //divide by jet pt
  cc[0][0]=cc[0][0]/jt.Pt();
  cc[0][1]=cc[0][1]/jt.Pt();
  cc[1][0]=cc[1][0]/jt.Pt();
  cc[1][1]=cc[1][1]/jt.Pt();
  //get eigenvalues
  TMatrixDEigen dd(cc);
  TVectorD rr=dd.GetEigenValuesRe();
  if(rr[0]!=0 && rr[0]!=0){
    quartic=sqrt(rr[0]*rr[0]+rr[1]*rr[1]);
    eccentricity=sqrt((rr[0]*rr[0]-rr[1]*rr[1])/rr[0]);
    pflow=4*rr[0]*rr[1]/pow(rr[0]+rr[1],2.0);
  }else{
    quartic=0.0;
    eccentricity=0.0;
    pflow=0.0;
  }
  return;
}

void softjetFinder::get_trkNN(int jetid, jetdata data){
  //tmphist[1]->Reset();
  
  TLorentzVector jt(data.jetpx[jetid],
		    data.jetpy[jetid],
		    data.jetpz[jetid],
		    data.jete[jetid]);

  //integrated jet shape
  //cout << data.tr_npart << " " << data.nu_npart << endl;
  double mom[50][3];
  for(int i=0;i<50;i++){
    for(int j=0;j<3;j++){
      mom[i][j]=0.0;
    }
  }

  TLorentzVector tr,tr2;
  if(data.tr_npart>50) data.tr_npart=50;
  if(data.nu_npart>50) data.nu_npart=50;

  for(int n=0;n<data.tr_npart;n++){
    if(data.tr_e[n]<0.0) continue;
    if(data.tr_e[n]<sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0)))
      data.tr_e[n]=sqrt(pow(data.tr_px[n],2.0)+pow(data.tr_py[n],2.0)+pow(data.tr_pz[n],2.0));
    
    tr.SetPxPyPzE(data.tr_px[n],
		  data.tr_py[n],
		  data.tr_pz[n],
		  data.tr_e[n]);
    for(int m=0;m<50;m++){
      if(TMath::Cos(jt.Angle(tr.Vect()))>1.0-(m+1)*0.01){
	mom[m][0]+=tr.Px();
	mom[m][1]+=tr.Py();
	mom[m][2]+=tr.Pz();
      } 
    }
  }
  
  for(int n=0;n<data.nu_npart;n++){
    if(data.nu_e[n]<0.0) continue;
    if(data.nu_e[n]<sqrt(pow(data.nu_px[n],2.0)+pow(data.nu_py[n],2.0)+pow(data.nu_pz[n],2.0)))
      data.nu_e[n]=sqrt(pow(data.nu_px[n],2.0)+pow(data.nu_py[n],2.0)+pow(data.nu_pz[n],2.0));
    tr2.SetPxPyPzE(data.nu_px[n],
		   data.nu_py[n],
		   data.nu_pz[n],
		   data.nu_e[n]);
    for(int m=0;m<50;m++){
      if(TMath::Cos(jt.Angle(tr2.Vect()))>1.0-(m+1)*0.01){
	mom[m][0]+=tr2.Px();
	mom[m][1]+=tr2.Py();
	mom[m][2]+=tr2.Pz();
      } 
    }
  }

  //for(int m=0;m<50;m++){
  //  tmphist[1]->SetBinContent(50-m, sqrt(pow(mom[m][0],2.0)+pow(mom[m][1],2.0)+pow(mom[m][2],2.0))/jt.P());
  //}

  //get bin content
  for(int m=0;m<20;m++){
    fvar[m]=sqrt(pow(mom[m][0],2.0)+pow(mom[m][1],2.0)+pow(mom[m][2],2.0))/jt.P();
    //tmphist[1]->GetBinContent(50-m);
  }

  return;
}

double softjetFinder::get_MVAvalue(int jetid, jetdata data){
  TLorentzVector j1,j2,j1j2;
  //cout << "start MVAvalue" << endl;
  //cout << "LorentzVector done" << endl;
  j1.SetPxPyPzE(data.jetpx[jetid],
		data.jetpy[jetid],
		data.jetpz[jetid],
		data.jete[jetid]);
  j2.SetPxPyPzE(data.jetpx[data.nearjetid[jetid]],
		data.jetpy[data.nearjetid[jetid]],
		data.jetpz[data.nearjetid[jetid]],
		data.jete[data.nearjetid[jetid]]);
  j1j2=j1+j2;

  //checking new method!!!
  //input to soft jet NN
  get_trkNN(jetid, data);
  //for(int i=0;i<20;i++) cout << "fvar: " << i << fvar[i] << endl;      
  sjvar[0]=j1.Eta();
  sjvar[1]=0.0;
  if(data.emcal+data.hadcal!=0.0) sjvar[1]=data.emcal/(data.emcal+data.hadcal);
  //other vals
  sjvar[2]=data.trkmaxptrel[jetid]/data.jete[jetid];
  sjvar[3]=data.neumaxptrel[jetid]/data.jete[jetid];
  sjvar[4]=data.trkptrelsum[jetid]/data.jete[jetid];
  sjvar[5]=data.neuptrelsum[jetid]/data.jete[jetid];
  sjvar[6]=data.ntrack[jetid]+data.nneu[jetid];
  
  if(TMath::Cos(j1.Vect().Angle(j2.Vect()))<0.75){
    sjvar[7]=disc_cluster2->EvaluateMVA("MLP_track2_new");
  }else{
    sjvar[7]=disc_cluster->EvaluateMVA("MLP_track_new");
  }
  get_Bjet(jetid, data);
  sjvar[8]=Bjet;
  get_Aa(jetid, data);
  sjvar[9]=Aa/data.jete[jetid];
  
  get_twomom(jetid, data);
  sjvar[10]=Tbeta; 
  get_geomom(jetid, data);
  sjvar[11]=quartic;
  sjvar[12]=eccentricity;
  sjvar[13]=pflow;
  
  sjvar[14]=j1.E();
  
  double tmpdisc=0.0;
  if(TMath::Cos(j1.Vect().Angle(j2.Vect()))<0.75){
    tmpdisc=disc_softjet2->EvaluateMVA("MLP_softjet2_new");
  }else{
    tmpdisc=disc_softjet->EvaluateMVA("MLP_softjet_new");
  }

  return tmpdisc;
}

double calWeight(int leptype, int proctype,double evtnum,double sigcs){
  double weight=1.0;
  double lum=500.0; //500 fb-1  1 atbarn(unit fb-1)
  double cs=1.0; //cross section(unit fb)
 
  switch(leptype){
  case 0:  //case 0 is reserved for test
    weight=1.0/evtnum;  //just normalize to 1
    break;
  case 1: //normalize their cs and luminosity
    switch(proctype){
    case 0:  //signal - qqhh
      cs=sigcs;   //0.12502408;
      break;
    case 1:  //signal - qqhh
      cs=sigcs;   //0.12502408;
      break;
    case 2:  //signal - qqhh
      cs=sigcs;   //0.12502408;
      break;

    //eLpR
    case 3: //l+jets
      cs=50.296;   //bbn1e1du
      break;
    case 4: //l+jets
      cs=49.8838;   //bbn2e2du
      break;
    case 5: //l+jets
      cs=49.8330;   //bbn3e3du
      break;
    case 6: //l+jets 
      cs=50.1043;    //bbn1e1sc
      break;
    case 7: //l+jets
      cs=49.8201;   //bbn2e2sc
      break;
    case 8: //l+jets
      cs=49.8330;   //bbn3e3sc
      break;
    case 9: //all hadronic
      cs=126.0276;  //bbuddu
      break;
    case 10: //all hadronic
      cs=125.3797;   //bbcsdu
      break;
    case 11: //all hadronic
      cs=126.1838;    //bbcssc
      break;
    case 12: //ttbar
      cs=126.0276;  //bbuddu2
      break;
    case 13: //ttbar
      cs=125.3797;   //bbcsdu2
      break;
    case 14: //ttbar
      cs=125.3797;   //bbudsc
      break;
    case 15: //ttbar
      cs=126.1838;    //bbcssc2
      break;
    case 16: //ttbar
      cs=126.0276;  //bbuddu3
      break;
    case 17: //ttbar
      cs=125.3797;   //bbcsdu3
      break;
    case 18: //ttbar
      cs=125.3797;   //bbudsc2
      break;
    case 19: //ttbar
      cs=126.1838;    //bbcssc3
      break;
    case 20: //bbnlln
      cs=167.1098;
      break;
    case 21: //bbnlln
      cs=167.1098;
      break;
    case 22: //bbnlln
      cs=167.1098;
      break;
    case 23: //bbnlln
      cs=167.1098;
      break;
    case 24: //ttbar +bb
      cs=0.26902;
      break;
    case 25: //ttcc
      cs=0.2923;
      break;
    case 26: //ttqq
      cs=0.8278;
      break;
    case 27: //ttz
      cs=1.6535;
      break;
    case 28: //ttz
      cs=1.6535;
      break;
    case 29: //ttbarh
      cs=0.2396;
      break;
    case 30: //zzh
      cs=1.2873;
      break;
    case 31: //zzh
      cs=1.2873;
      break;
    case 32: //zzz
      cs=3.0734;
      break;
  
    //eRpL
    case 33: //l+jets
      cs=20.03;   //bbn1e1du
      break;
    case 34: //l+jets
      cs=20.13;   //bbn2e2du
      break;
    case 35: //l+jets
      cs=20.13;   //bbn3e3du
      break;
    case 36: //l+jets 
      cs=20.03;    //bbn1e1sc
      break;
    case 37: //l+jets
      cs=20.00;   //bbn2e2sc
      break;
    case 38: //l+jets
      cs=20.04;   //bbn3e3sc
      break;
    case 39: //bbnlln
      cs=71.0;
      break;
    case 40: //bbnlln
      cs=71.0;
      break;
    case 41: 
      cs=0.4639;   //ttz
      break;
    case 42:  
      cs=0.09259;    //tth
      break;
    case 43:
      cs=0.5275;   //zzh
      break;
    case 44: //zzz
      cs=0.9808;
      break;
    case 99: //qqWW
      cs=20.1741;
      break;
    }
    
    weight=cs*lum/evtnum;
    break;
  }

  return weight;
}

double calWeight_sig(int leptype, int proctype,double evtnum,double sigcs){
  double weight=1.0;
  double lum=500.0; //500 fb-1  1 atbarn(unit fb-1)
  double cs=1.0; //cross section(unit fb)
 
  switch(leptype){
  case 0:  //case 0 is reserved for test
    weight=1.0/evtnum;  //just normalize to 1
    break;
  case 1: //normalize their cs and luminosity
    switch(proctype){
      //----------- for eLpR ----------
    case 0:  //signal - qqhh
      cs=0.125;   //0.12502408;
      break;
    case 1:  //signal - qqhh
      cs=0.125;   //0.12502408;
      break;
    case 2:  //signal - qqhh
      cs=0.125;   //0.12502408;
      break;
    case 3: //signal - bbhh
      cs=0.045;   //bbn1e1du
      break;
    case 4: //signal - bbhh
      cs=0.045;   //bbn2e2du
      break;
    case 5: //signal - bbhh
      cs=0.045;   //bbn3e3du
      break;
    case 6: //signal - cchh
      cs=0.035;    //bbn1e1sc
      break;
    case 7: //signal - cchh
      cs=0.035;   //bbn2e2sc
      break;
    case 8: //signal - cchh
      cs=0.035;   //bbn3e3sc
      break;
    case 9: //signal - llhh
      cs=0.03071;  //bbuddu
      break;
      //for eRpL
    case 10: //signal - qqhh
      cs=0.08019;   //bbcsdu
      break;
    case 11: //signal - bbhh
      cs=0.02732;    //bbcssc
      break;
    case 12: //signal - cchh
      cs=0.02249;
      break;
    case 13: //signal - llhh
      cs=0.01961;
      break;
    case 14: //signal - n1n1hh
      cs=0.01798;    //bbcssc
      break;
    case 15: //signal - n2n2hh
      cs=0.01336;
      break;
    case 16: //signal - n3n3hh
      cs=0.01336;
      break;
    }
    
    weight=cs*lum/evtnum;
    break;
  }

  return weight;
}

TH1F* make_SNplot(int proctype, TH1F *s, TH1F *b){
  int nbins=s->GetNbinsX();
  double interval=s->GetBinWidth(1);  
  double min=s->GetBinLowEdge(1);
  double max=s->GetBinLowEdge(nbins)+interval;
  
  double error=interval/2.0;
  
  //cout << min << " " << max<< endl;
  
  //cal s/sqrt(s+n) for each bin
  double xval[500],yval[500];
  double maxyval=-10.0;
  for(Int_t i=1;i<=nbins;i++){
    double tmps=s->Integral(i,nbins+1,"");
    double tmpb=b->Integral(i,nbins+1,"");
    yval[i-1]=0.0;
    if(tmps+tmpb!=0.0) yval[i-1]=tmps/sqrt(tmps+tmpb);
    if(maxyval<yval[i-1]) maxyval=yval[i-1];
    xval[i-1]=min+(i-1)*interval+error;
  }
  
  //cal. error
  double xerr[500],yerr[500];
  for(Int_t i=1;i<=nbins;i++){
    xerr[i-1]=0.0;  //error;
    yerr[i-1]=0.0;
  }
  
  //make histogram (graph is not good for making good graph)
  string tmpstr="ret_" + itos(proctype);  
  TH1F* sn=new TH1F(tmpstr.c_str(),"",nbins,min,max);
  for(Int_t i=1;i<=nbins;i++){
    sn->SetBinContent(i,yval[i-1]);
  }
  
  sn->SetMarkerStyle(20);
  sn->SetMarkerSize(1.8);
  sn->SetLineWidth(2);
  sn->GetXaxis()->SetTitle(s->GetXaxis()->GetTitle());
  sn->GetYaxis()->SetTitle("s/#sqrt{s+b}");
  sn->SetTitle("");
  sn->SetMaximum(1.80*maxyval);
  sn->SetMinimum(0.0);
  
  //using color
  switch(proctype){
  case 0:
    sn->SetMarkerColor(52);
    sn->SetLineColor(52);
    break;
  case 1:
    sn->SetMarkerColor(54);
    sn->SetLineColor(54);
    break;
  case 2:
    sn->SetMarkerColor(56);
    sn->SetLineColor(56);
    break;
  case 3:
    sn->SetMarkerColor(58);
    sn->SetLineColor(58);
    break;
  case 4:
    sn->SetMarkerColor(60);
    sn->SetLineColor(60);
    break;
  case 5:
    sn->SetMarkerColor(62);
    sn->SetLineColor(62);
    break;
  case 6:
    sn->SetMarkerColor(64);
    sn->SetLineColor(64);
    break;
  case 7:
    sn->SetMarkerColor(66);
    sn->SetLineColor(66);
    break;
  case 8:
    sn->SetMarkerColor(68);
    sn->SetLineColor(68);
    break;
  case 9:
    sn->SetMarkerColor(70);
    sn->SetLineColor(70);
    break;
  default:
    sn->SetMarkerColor(51);
    sn->SetLineColor(51);
    break;
  }
  
  return sn;
}

TH1F* make_Splot(int proctype, TH1F *s){
  int nbins=s->GetNbinsX();
  double interval=s->GetBinWidth(1);  
  double min=s->GetBinLowEdge(1);
  double max=s->GetBinLowEdge(nbins)+interval;
  
  double error=interval/2.0;
  
  //cout << min << " " << max<< endl;
  
  //cal s/sqrt(s+n) for each bin
  double xval[500],yval[500];
  double totval=s->Integral();
  for(Int_t i=1;i<=nbins;i++){
    double tmps=s->Integral(i,nbins+1,"");
    yval[i-1]=0.0;
    if(totval!=0.0) yval[i-1]=tmps/totval;
    xval[i-1]=min+(i-1)*interval+error;
  }
  
  //cal. error
  double xerr[500],yerr[500];
  for(Int_t i=1;i<=nbins;i++){
    xerr[i-1]=0.0;   //error;
    yerr[i-1]=0.0;
  }
  
  //make graph
  string tmpstr="ret_" + itos(proctype);  
  TH1F* sn=new TH1F(tmpstr.c_str(),"",nbins,min,max);
  for(Int_t i=1;i<=nbins;i++){
    sn->SetBinContent(i,yval[i-1]);
  }
  sn->SetMarkerStyle(22);
  sn->SetMarkerSize(1.8);
  sn->SetLineWidth(2);
  sn->GetXaxis()->SetTitle(s->GetXaxis()->GetTitle());
  sn->GetYaxis()->SetTitle("survival rate");
  sn->SetTitle("");
  sn->SetMaximum(1.10);
  sn->SetMinimum(0.0);

  //using color
  switch(proctype){
  case 0:
    sn->SetMarkerColor(2);
    sn->SetLineColor(2);
    break;
  case 1:
    sn->SetMarkerColor(3);
    sn->SetLineColor(3);
    break;
  case 2:
    sn->SetMarkerColor(4);
    sn->SetLineColor(4);
    break;
  case 3:
    sn->SetMarkerColor(6);
    sn->SetLineColor(6);
    break;
  case 4:
    sn->SetMarkerColor(7);
    sn->SetLineColor(7);
    break;
  case 5:
    sn->SetMarkerColor(8);
    sn->SetLineColor(8);
    break;
  case 6:
    sn->SetMarkerColor(9);
    sn->SetLineColor(9);
    break;
  case 7:
    sn->SetMarkerColor(11);
    sn->SetLineColor(11);
    break;
  case 8:
    sn->SetMarkerColor(12);
    sn->SetLineColor(12);
    break;
  case 9:
    sn->SetMarkerColor(13);
    sn->SetLineColor(13);
    break;
  default:
    sn->SetMarkerColor(1);
    sn->SetLineColor(1);
    break;
  }
  
  return sn;
}

TH1F* make_SNplotRev(int proctype, TH1F *s, TH1F *b){
  int nbins=s->GetNbinsX();
  double interval=s->GetBinWidth(1);  
  double min=s->GetBinLowEdge(1);
  double max=s->GetBinLowEdge(nbins)+interval;
  
  double error=interval/2.0;
  
  //cout << min << " " << max<< endl;
  
  //cal s/sqrt(s+n) for each bin
  double xval[500],yval[500];
  double maxyval=-10.0;
  for(Int_t i=1;i<=nbins;i++){
    double tmps=s->Integral(0,i,"");
    double tmpb=b->Integral(0,i,"");
    yval[i-1]=0.0;
    if(tmps+tmpb!=0.0) yval[i-1]=tmps/sqrt(tmps+tmpb);
    if(maxyval<yval[i-1]) maxyval=yval[i-1];
    xval[i-1]=min+(i-1)*interval+error;
  }
  
  //cal. error
  double xerr[500],yerr[500];
  for(Int_t i=1;i<=nbins;i++){
    xerr[i-1]=0.0;   //error;
    yerr[i-1]=0.0;
  }
  
  //make graph
  string tmpstr="ret_" + itos(proctype);  
  TH1F* sn=new TH1F(tmpstr.c_str(),"",nbins,min,max);
  for(Int_t i=1;i<=nbins;i++){
    sn->SetBinContent(i,yval[i-1]);
  }
  sn->SetMarkerStyle(20);
  sn->SetMarkerSize(1.8);
  sn->SetLineWidth(2);
  sn->GetXaxis()->SetTitle(s->GetXaxis()->GetTitle());
  sn->GetYaxis()->SetTitle("s/#sqrt{s+b}");
  sn->SetTitle("");
  sn->SetMaximum(1.80*maxyval);
  sn->SetMinimum(0.0);
  
  //using color
  switch(proctype){
  case 0:
    sn->SetMarkerColor(52);
    sn->SetLineColor(52);
    break;
  case 1:
    sn->SetMarkerColor(53);
    sn->SetLineColor(53);
    break;
  case 2:
    sn->SetMarkerColor(54);
    sn->SetLineColor(54);
    break;
  case 3:
    sn->SetMarkerColor(56);
    sn->SetLineColor(56);
    break;
  case 4:
    sn->SetMarkerColor(57);
    sn->SetLineColor(57);
    break;
  case 5:
    sn->SetMarkerColor(58);
    sn->SetLineColor(58);
    break;
  case 6:
    sn->SetMarkerColor(59);
    sn->SetLineColor(59);
    break;
  case 7:
    sn->SetMarkerColor(60);
    sn->SetLineColor(60);
    break;
  case 8:
    sn->SetMarkerColor(61);
    sn->SetLineColor(61);
    break;
  case 9:
    sn->SetMarkerColor(62);
    sn->SetLineColor(62);
    break;
  default:
    sn->SetMarkerColor(51);
    sn->SetLineColor(51);
    break;
  }
  
  return sn;
}

TH1F* make_SplotRev(int proctype, TH1F *s){
  int nbins=s->GetNbinsX();
  double interval=s->GetBinWidth(1);  
  double min=s->GetBinLowEdge(1);
  double max=s->GetBinLowEdge(nbins)+interval;
  
  double error=interval/2.0;
  
  //cout << min << " " << max<< endl;
  
  //cal s/sqrt(s+n) for each bin
  double xval[500],yval[500];
  double totval=s->Integral();
  for(Int_t i=1;i<=nbins;i++){
    double tmps=s->Integral(0,i,"");
    yval[i-1]=0.0;
    if(totval!=0.0) yval[i-1]=tmps/totval;
    xval[i-1]=max-(i-1)*interval+error;
  }
  
  //cal. error
  double xerr[500],yerr[500];
  for(Int_t i=1;i<=nbins;i++){
    xerr[i-1]=0.0;   //error;
    yerr[i-1]=0.0;
  }
  
  

  //make graph
  string tmpstr="ret_" + itos(proctype);  
  TH1F* sn=new TH1F(tmpstr.c_str(),"",nbins,min,max);
  for(Int_t i=1;i<=nbins;i++){
    sn->SetBinContent(i,yval[i-1]);
  }
  sn->SetMarkerStyle(22);
  sn->SetMarkerSize(1.8);
  sn->SetLineWidth(2);
  sn->GetXaxis()->SetTitle(s->GetXaxis()->GetTitle());
  sn->GetYaxis()->SetTitle("survival rate");
  sn->SetTitle("");
  sn->SetMaximum(1.80);
  sn->SetMinimum(0.0);

  //using color
  switch(proctype){
  case 0:
    sn->SetMarkerColor(2);
    sn->SetLineColor(2);
    break;
  case 1:
    sn->SetMarkerColor(3);
    sn->SetLineColor(3);
    break;
  case 2:
    sn->SetMarkerColor(4);
    sn->SetLineColor(4);
    break;
  case 3:
    sn->SetMarkerColor(6);
    sn->SetLineColor(6);
    break;
  case 4:
    sn->SetMarkerColor(7);
    sn->SetLineColor(7);
    break;
  case 5:
    sn->SetMarkerColor(8);
    sn->SetLineColor(8);
    break;
  case 6:
    sn->SetMarkerColor(9);
    sn->SetLineColor(9);
    break;
  case 7:
    sn->SetMarkerColor(11);
    sn->SetLineColor(11);
    break;
  case 8:
    sn->SetMarkerColor(12);
    sn->SetLineColor(12);
    break;
  case 9:
    sn->SetMarkerColor(13);
    sn->SetLineColor(13);
    break;
  default:
    sn->SetMarkerColor(1);
    sn->SetLineColor(1);
    break;
  }
  
  return sn;
}

TH1F* make_SNplotAbs(int proctype, TH1F *s, TH1F *b){
  int nbins=(int)s->GetNbinsX()/2;
  double interval=s->GetBinWidth(1);  
  double min=0.0;
  double max=s->GetBinLowEdge(nbins*2)+interval;
  
  double error=interval/2.0;
  
  //cout << min << " " << max<< endl;
  
  //cal s/sqrt(s+n) for each bin
  double xval[500],yval[500];
  double maxyval=-10.0;
  for(Int_t i=0;i<=nbins;i++){
    double tmps=s->Integral(i,nbins*2+1-i,"");
    double tmpb=b->Integral(i,nbins*2+1-i,"");
    yval[i]=0.0;
    if(tmps+tmpb!=0.0) yval[i]=tmps/sqrt(tmps+tmpb);
    if(maxyval<yval[i]) maxyval=yval[i];
    xval[i]=max-i*interval+error;
  }
  
  //cal. error
  double xerr[500],yerr[500];
  for(Int_t i=0;i<=nbins;i++){
    xerr[i]=0.0;   //error;
    yerr[i]=0.0;
  }
  
  //make graph
  string tmpstr="ret_" + itos(proctype);  
  TH1F* sn=new TH1F(tmpstr.c_str(),"",nbins,min,max);
  for(Int_t i=1;i<=nbins;i++){
    sn->SetBinContent(i,yval[nbins-i+1]);
  }
  sn->SetMarkerStyle(20);
  sn->SetMarkerSize(1.8);
  sn->SetLineWidth(2);
  sn->GetXaxis()->SetTitle(s->GetXaxis()->GetTitle());
  sn->GetYaxis()->SetTitle("s/#sqrt{s+b}");
  sn->SetTitle("");
  sn->SetMaximum(1.8*maxyval);
  sn->SetMinimum(0.0);
  
  //using color
  switch(proctype){
  case 0:
    sn->SetMarkerColor(52);
    sn->SetLineColor(52);
    break;
  case 1:
    sn->SetMarkerColor(53);
    sn->SetLineColor(53);
    break;
  case 2:
    sn->SetMarkerColor(54);
    sn->SetLineColor(54);
    break;
  case 3:
    sn->SetMarkerColor(56);
    sn->SetLineColor(56);
    break;
  case 4:
    sn->SetMarkerColor(57);
    sn->SetLineColor(57);
    break;
  case 5:
    sn->SetMarkerColor(58);
    sn->SetLineColor(58);
    break;
  case 6:
    sn->SetMarkerColor(60);
    sn->SetLineColor(60);
    break;
  case 7:
    sn->SetMarkerColor(62);
    sn->SetLineColor(62);
    break;
  case 8:
    sn->SetMarkerColor(64);
    sn->SetLineColor(64);
    break;
  case 9:
    sn->SetMarkerColor(66);
    sn->SetLineColor(66);
    break;
  default:
    sn->SetMarkerColor(51);
    sn->SetLineColor(51);
    break;
  }
  
  return sn;
}

TH1F* make_SplotAbs(int proctype, TH1F *s){
  int nbins=(int)s->GetNbinsX()/2;
  double interval=s->GetBinWidth(1);  
  double min=0.0;
  double max=s->GetBinLowEdge(nbins*2)+interval;
  
  double error=interval/2.0;
  
  //cout << min << " " << max<< endl;
  
  //cal s/sqrt(s+n) for each bin
  double xval[500],yval[500];
  double totval=s->Integral();
  for(Int_t i=0;i<=nbins;i++){
    double tmps=s->Integral(i,nbins*2+1-i,"");
    yval[i]=0.0;
    if(totval!=0.0) yval[i]=tmps/totval;
    xval[i]=max-i*interval+error;
  }
  
  //cal. error
  double xerr[500],yerr[500];
  for(Int_t i=0;i<=nbins;i++){
    xerr[i]=0.0;   //error;
    yerr[i]=0.0;
  }
  
  //make graph
  string tmpstr="ret_" + itos(proctype);  
  TH1F* sn=new TH1F(tmpstr.c_str(),"",nbins,min,max);
  for(Int_t i=1;i<=nbins;i++){
    sn->SetBinContent(i,yval[nbins-i+1]);
  }
  sn->SetMarkerStyle(22);
  sn->SetMarkerSize(1.8);
  sn->SetLineWidth(2);
  sn->GetXaxis()->SetTitle(s->GetXaxis()->GetTitle());
  sn->GetYaxis()->SetTitle("survival rate");
  sn->SetTitle("");
  sn->SetMaximum(1.80);
  sn->SetMinimum(0.0);

  //using color
  switch(proctype){
  case 0:
    sn->SetMarkerColor(2);
    sn->SetLineColor(2);
    break;
  case 1:
    sn->SetMarkerColor(3);
    sn->SetLineColor(3);
    break;
  case 2:
    sn->SetMarkerColor(4);
    sn->SetLineColor(4);
    break;
  case 3:
    sn->SetMarkerColor(6);
    sn->SetLineColor(6);
    break;
  case 4:
    sn->SetMarkerColor(7);
    sn->SetLineColor(7);
    break;
  case 5:
    sn->SetMarkerColor(8);
    sn->SetLineColor(8);
    break;
  case 6:
    sn->SetMarkerColor(9);
    sn->SetLineColor(9);
    break;
  case 7:
    sn->SetMarkerColor(11);
    sn->SetLineColor(11);
    break;
  case 8:
    sn->SetMarkerColor(12);
    sn->SetLineColor(12);
    break;
  case 9:
    sn->SetMarkerColor(13);
    sn->SetLineColor(13);
    break;
  default:
    sn->SetMarkerColor(1);
    sn->SetLineColor(1);
    break;
  }
  
  return sn;
}

void setColor(TH1F *h, int type){

  switch(type){
  case 0:  //signal
    h->SetLineColor(2);
    break;
  case 1:  //ttbar - lepton+jets
    h->SetLineColor(3);
    break;
  case 2:  //ttbar - all hadronic
    h->SetLineColor(4);
    break;
  case 3:  //ttbar - dilepton
    h->SetLineColor(6);
    break;
  case 4:  //ttbar + bb
    h->SetLineColor(7);
    break;
  case 5: //ttz
    h->SetLineColor(8);
    break;
  case 6: //zzh
    h->SetLineColor(9);
    break;
  case 7: //e2e2bbbb
    h->SetLineColor(11);
    break;
  case 8:
    h->SetLineColor(12);
    break;
  case 9:
    h->SetLineColor(1);
    break;
  case 10:
    h->SetLineColor(1);
    break;
  case 11:
    h->SetLineColor(1);
    break;
  case 12:
    h->SetLineColor(1);
    break;
  }
  h->GetYaxis()->SetLabelOffset(0.001);

  return;
}

string itos(int i)  
{
  stringstream s;
  s << i;
  return s.str();
}

string ftos(float f)  
{
  stringstream s;
  s << f;
  return s.str();
}

double calDelPhi(double phi1, double phi2){
  double a=0.0;

  //first, correct phi
  if(phi1<0.0) phi1=2*3.141592-fabs(phi1);
  if(phi2<0.0) phi2=2*3.141592-fabs(phi2);
  a=fabs(phi1-phi2);
  if(a>3.141592) a=2*3.141592-a;

  return a;
}

double cal_moment(jetdata data, double Evis, int order){
  double lgd=0.0,h=0.0,h0=0.0;
  
  for(int j=0;j<data.njets;j++){
    for(int l=0;l<data.njets;l++){
      //set tlorentzvector of 2 particles
      TLorentzVector jet1(data.jetpx[j],
			  data.jetpy[j],
			  data.jetpz[j],
			  data.jete[j]);
      
      TLorentzVector jet2(data.jetpx[l],
			  data.jetpy[l],
			  data.jetpz[l],
			  data.jete[l]);
      
      //cal angle
      double tangle=TMath::Cos(jet1.Angle(jet2.Vect()));
      double tprod=jet1.P()*jet2.P();
      
      //for normalize
      h0+=tprod/pow(Evis,2.0);

      //legendre polynominals
      switch(order){
      case 0:
	lgd=1.0;
	break;
      case 1:
	lgd=tangle;
	break;
      case 2:
	lgd=0.5*(3*pow(tangle,2.0)-1);
	break;
      case 3:
	lgd=0.5*(5*pow(tangle,3.0)-3*tangle);
	break;
      case 4:
	lgd=0.125*(35*pow(tangle,4.0)-30*pow(tangle,2.0)+3);
	break;
      case 5:
	lgd=0.125*(63*pow(tangle,5.0)-70*pow(tangle,3.0)+15*tangle);
	break;
      }
      
      h+=tprod*lgd/pow(Evis,2.0);
    }
  }

  //if(order!=0) h=h/h0;

  return h;
}

double cal_sphericity(jetdata data, int type){
  TMatrixD sab(3,3);
  for(int j=0;j<3;j++){
    for(int l=0;l<3;l++){
      sab[j][l]=0.0;
    }
  }
  
  double norm=0.0;
  for(int j=0;j<data.njets;j++){
    //set tlorentzvector of 2 particles
    TLorentzVector jet1(data.jetpx[j],
			data.jetpy[j],
			data.jetpz[j],
			data.jete[j]);
    
    sab[0][0]+=jet1.Px()*jet1.Px();
    sab[0][1]+=jet1.Px()*jet1.Py();
    sab[0][2]+=jet1.Px()*jet1.Pz();
    sab[1][0]+=jet1.Py()*jet1.Px();
    sab[1][1]+=jet1.Py()*jet1.Py();
    sab[1][2]+=jet1.Py()*jet1.Pz();
    sab[2][0]+=jet1.Pz()*jet1.Px();
    sab[2][1]+=jet1.Pz()*jet1.Py();
    sab[2][2]+=jet1.Pz()*jet1.Pz();
    
    norm+=jet1.P()*jet1.P();
  }
  
  //divide
  for(int j=0;j<3;j++){
    for(int l=0;l<3;l++){
      sab[j][l]=sab[j][l]/norm;
    }
  }
  
  TMatrixDEigen lambdas(sab);
  TVectorD lambda(3);
  lambda=lambdas.GetEigenValuesRe();
 
  if(type==0) return 1.5*(lambda[1]+lambda[2]);
  if(type==1) return 1.5*lambda[2];
  if(type==2) return 3.0*(lambda[0]*lambda[1]+lambda[1]*lambda[2]+lambda[2]*lambda[0]);
  if(type==3) return 27.0*lambda[0]*lambda[1]*lambda[2];

  return 0.0;
}

void jetclustering_2jets(int snjets,int cjets, float *var, double *jets){
  //start jet clustering until 2 jets
  vector<TLorentzVector> jetvec;
  for(int i=0;i<snjets;i++){
    TLorentzVector jt(var[34+4*i],
		      var[35+4*i],
		      var[36+4*i],
		      var[37+4*i]);
    jetvec.push_back(jt);
  }
  //cout << jetvec.size() << endl;

  TLorentzVector j1,j2;
  while(jetvec.size()>cjets){
    double okyval=1.0e20;
    int oki=0,okj=0;
    for(int i=0;i<jetvec.size();i++){
      for(int j=i+1;j<jetvec.size();j++){
	j1.SetPxPyPzE(jetvec[i].Px(),
		      jetvec[i].Py(),
		      jetvec[i].Pz(),
		      jetvec[i].E());
	j2.SetPxPyPzE(jetvec[j].Px(),
		      jetvec[j].Py(),
		      jetvec[j].Pz(),
		      jetvec[j].E());
	
	//cal. yvalue
	double tmpyval=TMath::Min(pow(j1.E(),2.0),pow(j2.E(),2.0))*(1-TMath::Cos(j1.Angle(j2.Vect())))
	  /pow(500.0,2.0);
	
	if(tmpyval<okyval){
	  okyval=tmpyval;
	  oki=i;
	  okj=j;
	}
      }
    }

    //start to add min. y value jets
    j1.SetPxPyPzE(jetvec[oki].Px(),
		  jetvec[oki].Py(),
		  jetvec[oki].Pz(),
		  jetvec[oki].E());
    j2.SetPxPyPzE(jetvec[okj].Px(),
		  jetvec[okj].Py(),
		  jetvec[okj].Pz(),
		  jetvec[okj].E());
    TLorentzVector jj=j1+j2;

    vector<TLorentzVector> tmpjetvec;
    for(int i=0;i<jetvec.size();i++){
      if(i==oki || i==okj) continue;
      tmpjetvec.push_back(jetvec[i]);
    }
    tmpjetvec.push_back(jj);
    
    jetvec.clear();
    
    for(int i=0;i<tmpjetvec.size();i++){
      jetvec.push_back(tmpjetvec[i]);
    }
  }

  //get vector
  if(cjets==2){
    if(jetvec[0].E()>jetvec[1].E()){
      jets[0]=jetvec[0].Px();
      jets[1]=jetvec[0].Py();
      jets[2]=jetvec[0].Pz();
      jets[3]=jetvec[0].E();
      jets[4]=jetvec[1].Px();
      jets[5]=jetvec[1].Py();
      jets[6]=jetvec[1].Pz();
      jets[7]=jetvec[1].E();
    }else{
      jets[0]=jetvec[1].Px();
      jets[1]=jetvec[1].Py();
      jets[2]=jetvec[1].Pz();
      jets[3]=jetvec[1].E();
      jets[4]=jetvec[0].Px();
      jets[5]=jetvec[0].Py();
      jets[6]=jetvec[0].Pz();
      jets[7]=jetvec[0].E();
    }
  }else if(cjets==3){
    if(jetvec[0].E()>TMath::Max(jetvec[1].E(),jetvec[2].E())){
      jets[0]=jetvec[0].Px();
      jets[1]=jetvec[0].Py();
      jets[2]=jetvec[0].Pz();
      jets[3]=jetvec[0].E();
      jets[4]=jetvec[1].Px();
      jets[5]=jetvec[1].Py();
      jets[6]=jetvec[1].Pz();
      jets[7]=jetvec[1].E();
      jets[8]=jetvec[2].Px();
      jets[9]=jetvec[2].Py();
      jets[10]=jetvec[2].Pz();
      jets[11]=jetvec[2].E();
    }else if(jetvec[1].E()>TMath::Max(jetvec[0].E(),jetvec[2].E())){
      jets[0]=jetvec[1].Px();
      jets[1]=jetvec[1].Py();
      jets[2]=jetvec[1].Pz();
      jets[3]=jetvec[1].E();
      jets[4]=jetvec[0].Px();
      jets[5]=jetvec[0].Py();
      jets[6]=jetvec[0].Pz();
      jets[7]=jetvec[0].E();
      jets[8]=jetvec[2].Px();
      jets[9]=jetvec[2].Py();
      jets[10]=jetvec[2].Pz();
      jets[11]=jetvec[2].E();
    }else{
      jets[0]=jetvec[2].Px();
      jets[1]=jetvec[2].Py();
      jets[2]=jetvec[2].Pz();
      jets[3]=jetvec[2].E();
      jets[4]=jetvec[0].Px();
      jets[5]=jetvec[0].Py();
      jets[6]=jetvec[0].Pz();
      jets[7]=jetvec[0].E();
      jets[8]=jetvec[1].Px();
      jets[9]=jetvec[1].Py();
      jets[10]=jetvec[1].Pz();
      jets[11]=jetvec[1].E();
    }
  }
  return;
}

double calThrust(int njets, double *var){

  vector<TVector3> _partMom;

  for(int j=0;j<8;j++){
    int tmpj=j;
    if(njets==6 && j==7) continue;

    double px = var[34+tmpj*4];
    double py = var[35+tmpj*4];
    double pz = var[36+tmpj*4];
    double pt = sqrt(px*px+py*py);
    double mom = sqrt(px*px+py*py+pz*pz);
    double cosTheta = pz/mom;

    if (pt<0.1) continue;
    if (pt<0.5) continue;
    if (cosTheta>0.997) continue;

    TVector3 a( px, py, pz );
    _partMom.push_back(a);
    
  }
  
  const int nwork=11,iFastMax = 4,iGood=2;  
  const float dConv=0.0001; // 0.0001
  int sgn;
  double theta=0,phi=0;
  double thp,thps,tds,tmax,dOblateness;
  vector<TVector3> TAxes(3),Fast(iFastMax+1),Workv(nwork);
  vector<double> Workf(nwork),dThrust(3);
  TVector3 tdi,tpr,mytest;
  
  tmax = 0;
  for ( unsigned int i=0; i < _partMom.size(); i++)
    tmax += _partMom[i].Mag();
  
  // pass = 0: find thrust axis
  // pass = 1: find major axis
  for ( int pass=0; pass <= 1; pass++ ) 
    {
      if ( pass == 1 )
	{
	  phi   = TAxes[0].Phi();
	  theta = TAxes[0].Theta();
	  for ( unsigned  int i = 0;i < _partMom.size(); i++)
	    {
	      _partMom[i].RotateZ(-phi);
	      _partMom[i].RotateY(-theta);
	    }
	  TAxes[0].SetXYZ(0,0,1);
	} // if pass == 1
      
      // Find the ifast highest momentum particles and
      // put the highest in Fast[0], next in Fast[1],....Fast[iFast-1].
      // Fast[iFast] is just a workspace.
      
      for ( unsigned  int i = 0; i < Fast.size(); i++ )
	Fast[i].SetXYZ(0,0,0);
      
      for ( unsigned int i = 0; i < _partMom.size(); i++ )
	{
	  for ( int ifast = iFastMax -1; ifast >= 0 ; ifast-- ) 
	    {
	      if (_partMom[i].Mag2() > Fast[ifast].Mag2() ) 
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
      
      for ( int iw = 0; iw < Workv.size(); iw++ ) 
	{
	  Workf[iw] = 0.;
	}
      int p = (int) min( iFastMax, (int)_partMom.size() ) - 1 ;
      int nc = 1 << p;
      for ( int n = 0; n < nc; n++ )
	{
	  tdi.SetXYZ(0,0,0);
	  for (unsigned int i = 0; i < min(iFastMax,nc) ; i++)
	    {
	      if ( (1 << (i+1)) * ( (n + (1<<i)) / (1<<(i+1)) ) >= n+1) //i+1 
		{ sgn = -1;} else {sgn = 1;}
	      tdi += sgn*Fast[i];
	      if (pass==1) tdi.SetZ(0);
	    } // for i 
	  tds = tdi.Mag2(); 
	  for (int iw = (int) min(n,9); iw >= 0; iw-- )
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
      for (int iw = 0; iw < min(nc,10) && nagree < iGood; iw++ )
	{
	  thp = 0;
	  thps = -99999.;
	  while ( thp > thps + dConv )
	    {
	      thps = thp;
	      if ( thp <= 1E-10 )
		{ tdi = Workv[iw]; } else { tdi=tpr; }
	      tpr.SetXYZ(0,0,0);
	      for ( unsigned int i = 0; i < _partMom.size(); i++ )
		{
		  sgn = (int) TMath::Sign(1.0,tdi.Dot(_partMom[i]));
		  tpr += sgn*_partMom[i];
		  if (pass == 1) { tpr.SetZ(0); } // ###
		} // for i 
	      thp = tpr.Mag()/tmax;
	    } // while 
	  // Save good axis. Try new initial axis until enough
	  // tries agree.
	  if ( thp < dThrust[pass] - dConv ) continue;
	  if ( thp > dThrust[pass] + dConv )
	    {
	      nagree = 0;
	      // 	      if (myrnd.flat() > 0.49999)
	      // 		{sgn = 1;} else {sgn=-1;}
	      sgn = 1; 
	      TAxes[pass] = (sgn/(tmax*thp))*tpr;
	      dThrust[pass] = thp;
	    } // if thp
	  nagree++;
	} // for iw (2)
    } // for pass ...
  
  // Find minor axis and value by orthogonality.
  //make random generator
  TRandom myrnd(12345);
  if (myrnd.Uniform(1.0) > 0.49999)
    {sgn = 1;} else {sgn=-1;}
  TAxes[2].SetXYZ( -sgn*TAxes[1].y(), sgn*TAxes[1].x(), 0);
  thp = 0.;
  for ( unsigned int i = 0; i < _partMom.size(); i++ )
    {
      thp += fabs(TAxes[2].Dot(_partMom[i]) );
    } // for i 
  dThrust[2] = thp/tmax;
  
  // Rotate back to original coordinate system.
  for ( unsigned int i = 0;i < TAxes.size(); i++)
    {
      TAxes[i].RotateY(theta); 
      TAxes[i].RotateZ(phi); 
    }
  dOblateness = dThrust[1] - dThrust[2];
  
  //_principleThrustValue = dThrust[0];
  //_majorThrustValue     = dThrust[1];
  //_minorThrustValue     = dThrust[2];
  //_principleThrustAxis  =   TAxes[0];
  //_majorThrustAxis      =   TAxes[1];
  //_minorThrustAxis      =   TAxes[2];
  
  return dThrust[0];
  
}

double corrEvtWeight(double evtweight, int proctype){
  //correct the event rate for merging eLpR and eRpL
  //(0.8,-0.3)             (1.0,-1.0)      (-1.0 1.0)
  //--- signal ---
  //bbhh   0.02732         0.045           0.02732  ??
  //cchh   0.02131         0.035           0.02249  ??
  //qqhh   0.07595         0.125           0.08019  ??
  //llhh   0.01873         0.03071         0.01961  ??
  //--- ttbar_lepjets ---
  //bbn1e1du  30.1402      50.296          20.03
  //bbn2e2du  29.9388      49.8838         20.13
  //bbn3e3du  29.8948      49.8330         20.13
  //bbn1e1sc  30.0770      50.1043         20.03
  //bbn2e2sc  29.8671      49.8201         20.00
  //bbn3e3sc  29.8692      49.8330         20.04
  //--- ttbar_allhad ---
  //bbuddu    75.2970      126.0276        na
  //bbcsdu    75.6460      125.3797        na
  //bbudsc    75.4263      125.3797        na
  //bbcssc    75.3600      126.1838        na
  //--- ttbar_dilepton ---
  //bbnlln    100.2875     167.1098        71.0
  //--- ttbar + gamma ---
  //ttbb     0.1604        0.26902         na
  //ttcc     0.1747        0.2923          na
  //ttqq     0.4935        0.8278          na
  //--- ttbar + Z --- 
  //ttz      0.9830        1.6535          0.4639
  //--- ttbar + H ---
  //tth      0.1433        0.2396          0.09259
  //--- ZZ + H ---
  //zzh      0.7710        1.2873          0.5275
  //--- ZZZ ---
  //zzz      1.830         3.0734          0.9808

  //eLpR values
  double sigma1[100];
  //signal
  sigma1[0]=0.045;
  sigma1[1]=0.035;
  sigma1[2]=0.125;
  sigma1[3]=0.03071;
  //ttbar_lepjets
  sigma1[4]=50.296;
  sigma1[5]=49.8838;
  sigma1[6]=49.8330;
  sigma1[7]=50.1043;
  sigma1[8]=49.8201;
  sigma1[9]=49.8330;
  //ttbar_allhad
  sigma1[10]=126.0276;
  sigma1[11]=125.3797;
  sigma1[12]=125.3797;
  sigma1[13]=126.1838;
  //ttbar_dilepton
  sigma1[14]=167.1098;
  //ttbar+gamma
  sigma1[15]=0.26902;
  sigma1[16]=0.2923;
  sigma1[17]=0.8278;
  //ttbar+Z
  sigma1[18]=1.6535;
  //ttbar+H
  sigma1[19]=0.2796;
  //ZZ+H
  sigma1[20]=1.2873;
  //ZZ+Z
  sigma1[21]=3.0734;

  //eRpL values
  double sigma2[100];
  //signal
  sigma2[0]=0.02732;
  sigma2[1]=0.02249;
  sigma2[2]=0.08019;
  sigma2[3]=0.01961;
  //ttbar_lepjets
  sigma2[4]=20.03;
  sigma2[5]=20.13;
  sigma2[6]=20.13;
  sigma2[7]=20.03;
  sigma2[8]=20.00;
  sigma2[9]=20.04;
  //ttbar_allhad
  sigma2[10]=0.0;
  sigma2[11]=0.0;
  sigma2[12]=0.0;
  sigma2[13]=0.0;
  //ttbar_dilepton
  sigma2[14]=71.00;
  //ttbar+gamma
  sigma2[15]=0.0;
  sigma2[16]=0.0;
  sigma2[17]=0.0;
  //ttbar+Z
  sigma2[18]=0.4639;
  //ttbar+H
  sigma2[19]=0.09259;
  //ZZ+H
  sigma2[20]=0.5275;
  //ZZ+Z
  sigma2[21]=0.9808;

  //merged value
  double sigma3[100];
  //signal
  sigma3[0]=0.03297;
  sigma3[1]=0.02649;
  sigma3[2]=0.09453;
  sigma3[3]=0.02316;
  //ttbar_lepjets
  sigma3[4]=30.1402;
  sigma3[5]=29.9388;
  sigma3[6]=29.8948;
  sigma3[7]=30.0770;
  sigma3[8]=29.8671;
  sigma3[9]=29.8692;
  //ttbar_allhad
  sigma3[10]=75.2970;
  sigma3[11]=75.6460;
  sigma3[12]=75.4263;
  sigma3[13]=75.3600;
  //ttbar_dilepton
  sigma3[14]=100.2875;
  //ttbar+gamma
  sigma3[15]=0.1604;
  sigma3[16]=0.1747;
  sigma3[17]=0.4935;
  //ttbar+Z
  sigma3[18]=0.9830;
  //ttbar+H
  sigma3[19]=0.1433;
  //ZZ+H
  sigma3[20]=0.7710;
  //ZZ+Z
  sigma3[21]=1.830;

  //final rate
  double rate[100];
  for(int i=0;i<22;i++){
    rate[i]=0.0;
    if(sigma1[i]-sigma2[i]!=0.0) rate[i]=(sigma3[i]-sigma2[i])/(sigma1[i]-sigma2[i]);
  }
  
  //cout << evtweight << endl;
  switch(proctype){
  case 0:   //signal
    //for 6 jets
    //qqhh
    if(fabs(evtweight-0.000209553)<1.0e-8) return rate[2];  //eLpR//
    if(fabs(evtweight-0.000242634)<1.0e-8) return 1.0-rate[2];  //eRpL  //
    //bbhh
    if(fabs(evtweight-7.82554e-05)<1.0e-8) return rate[0];   //eLpR//
    if(fabs(evtweight-9.07901e-05)<1.0e-8) return 1.0-rate[0];   //eLpR//
    //cchh
    if(fabs(evtweight-5.84926e-05)<1.0e-8) return rate[1];   //eLpR//
    if(fabs(evtweight-6.97936e-05)<1.0e-8) return 1.0-rate[1];    //eRpL//
    //llhh
    if(fabs(evtweight-6.19023e-05)<1.0e-8) return rate[3];   //eLpR//
    if(fabs(evtweight-5.47435e-05)<1.0e-8) return 1.0-rate[3];    //eRpL //
    //nnhh
    if(fabs(evtweight-9.00443e-05)<1.0e-8) return 1.0;  //eLpR//
    if(fabs(evtweight-7.12344e-05)<1.0e-8) return 1.0;  //eLpR//   
    if(fabs(evtweight-0.00017975)<1.0e-8) return 1.0;  //eLpR//   
    if(fabs(evtweight-0.000128429)<1.0e-8) return 1.0;  //eLpR//
    if(fabs(evtweight-9.59756e-05)<1.0e-8) return 1.0;  //eLpR//   
    if(fabs(evtweight-9.58118e-05)<1.0e-8) return 1.0;  //eLpR//  

    //for 4 jets
    //qqhh
    if(fabs(evtweight-0.000209553)<1.0e-8) return rate[2];  //eLpR//
    if(fabs(evtweight-0.000242634)<1.0e-8) return 1.0-rate[2];  //eRpL  //
    //bbhh
    if(fabs(evtweight-7.82554e-05)<1.0e-8) return rate[0];   //eLpR//
    if(fabs(evtweight-9.07901e-05)<1.0e-8) return 1.0-rate[0];   //eLpR//
    //cchh
    if(fabs(evtweight-5.84926e-05)<1.0e-8) return rate[1];   //eLpR//
    if(fabs(evtweight-6.97936e-05)<1.0e-8) return 1.0-rate[1];    //eRpL//
    //llhh
    if(fabs(evtweight-6.19023e-05)<1.0e-8) return rate[3];   //eLpR//
    if(fabs(evtweight-5.47435e-05)<1.0e-8) return 1.0-rate[3];    //eRpL //
 
    break;
  case 1:   //ttbar_lepjets
    //for 6 jets
    //bbn1e1du  3   
    if(fabs(evtweight-0.0487284)<1.0e-6) return rate[4]; //eLpR
    //bbn2e2du   4   
    if(fabs(evtweight-0.0489073)<1.0e-6) return rate[5];   //eLpR
    //bbn3e3du   5   
    if(fabs(evtweight-0.34153)<1.0e-6) return rate[6];   //eLpR
    //bbn1e1sc   6  
    if(fabs(evtweight-0.0699755)<1.0e-6) return rate[7];  //eLpR
    //bbn2e2sc   7 
    if(fabs(evtweight-0.0490384)<1.0e-6) return rate[8];   //eLpR
    //bbn3e3sc   8
    if(fabs(evtweight-0.0563704)<1.0e-6) return rate[9];   //eLpR
    //bbn1e1du   33  
    if(fabs(evtweight-0.054411)<1.0e-6) return 1.0-rate[4]; //eRpL
    //bbn2e2du   34   
    if(fabs(evtweight-0.050413)<1.0e-6) return 1.0-rate[5];   //eRpL
    //bbn3e3du   35   
    if(fabs(evtweight-0.0517973)<1.0e-6) return 1.0-rate[6];   //eRpL
    //bbn1e1sc   36  
    if(fabs(evtweight-0.0505627)<1.0e-6) return 1.0-rate[7];  //eRpL
    //bbn2e2sc   37 
    if(fabs(evtweight-0.0759267)<1.0e-6) return 1.0-rate[8];   //eRpL
    //bbn3e3sc   38
    if(fabs(evtweight-0.0519836)<1.0e-6) return 1.0-rate[9];   //eRpL
 
    //for 4 jets
    //bbn1e1du  3   
    if(fabs(evtweight-0.0486594)<1.0e-6) return rate[4]; //eLpR
    //bbn2e2du   4   
    if(fabs(evtweight-0.0489073)<1.0e-6) return rate[5];   //eLpR
    //bbn3e3du   5   
    if(fabs(evtweight-0.341453)<1.0e-6) return rate[6];   //eLpR
    //bbn1e1sc   6  
    if(fabs(evtweight-0.0699755)<1.0e-6) return rate[7];  //eLpR
    //bbn2e2sc   7 
    if(fabs(evtweight-0.0490384)<1.0e-6) return rate[8];   //eLpR
    //bbn3e3sc   8
    if(fabs(evtweight-0.0584459)<1.0e-6) return rate[9];   //eLpR
    //bbn1e1du   33  
    if(fabs(evtweight-0.054411)<1.0e-6) return 1.0-rate[4]; //eRpL
    //bbn2e2du   34   
    if(fabs(evtweight-0.050413)<1.0e-6) return 1.0-rate[5];   //eRpL
    //bbn3e3du   35   
    if(fabs(evtweight-0.0517973)<1.0e-6) return 1.0-rate[6];   //eRpL
    //bbn1e1sc   36  
    if(fabs(evtweight-0.0505627)<1.0e-6) return 1.0-rate[7];  //eRpL
    //bbn2e2sc   37 
    if(fabs(evtweight-0.0759267)<1.0e-6) return 1.0-rate[8];   //eRpL
    //bbn3e3sc   38
    if(fabs(evtweight-0.0519836)<1.0e-6) return 1.0-rate[9];   //eRpL

    //for 8 jets
    //bbn2e2du   4   
    if(fabs(evtweight-0.0490516)<1.0e-6) return rate[5];   //eLpR
    break;
  case 2:   //ttbar_allhad
    //for 6 jets
    //bbuddu    9
    if(fabs(evtweight-0.0509658)<1.0e-6) return rate[10];  //eLpR
    //bbcsdu   10
    if(fabs(evtweight-0.0384459)<1.0e-6) return rate[11];  //eLpR
    //bbcssc   11
    if(fabs(evtweight-0.0549975)<1.0e-6) return rate[13];   //eLpR
    //bbudsc   14
    if(fabs(evtweight-0.0384459)<1.0e-6) return rate[12];   //eLpR

     //for 4 jets
    //bbuddu    9
    if(fabs(evtweight-0.0510001)<1.0e-6) return rate[10];  //eLpR
    //bbcsdu   10
    if(fabs(evtweight-0.0389798)<1.0e-6) return rate[11];  //eLpR
    //bbcssc   11
    if(fabs(evtweight-0.0549926)<1.0e-6) return rate[13];   //eLpR
    //bbudsc   14
    if(fabs(evtweight-0.0389798)<1.0e-6) return rate[12];   //eLpR

     //for 8 jets
    //bbuddu    9
    if(fabs(evtweight-0.0510084)<1.0e-6) return rate[10];  //eLpR
    //bbcsdu   10
    if(fabs(evtweight-0.0390188)<1.0e-6) return rate[11];  //eLpR
    //bbcssc   11
    if(fabs(evtweight-0.0473653)<1.0e-6) return rate[13];   //eLpR
    //bbudsc   14
    if(fabs(evtweight-0.0390188)<1.0e-6) return rate[12];   //eLpR
   break;
  case 3:   //ttbar_dilepton
    //for 6 jets
    //bbnlln    20
    if(fabs(evtweight-0.0885308)<1.0e-5)  return rate[14];  //eLpR
    if(fabs(evtweight-0.0885777)<1.0e-5)  return rate[14];  //eLpR
    //bbnlln    39
    if(fabs(evtweight-0.130493)<1.0e-5)  return 1.0-rate[14];  //eRpL
    break;
  case 4:   //ttbar+gamma
    //for 6 jets
    //ttbarbb    24
    if(fabs(evtweight-0.00134789)<1.0e-8) return rate[15];   //eLpR
    //ttbarcc    25
    if(fabs(evtweight-0.00129176)<1.0e-8) return rate[16];   //eLpR
    //ttbarqq    26
    if(fabs(evtweight-0.00414377)<1.0e-8) return rate[17];   //eLpR
    break;
  case 5:   //tt + Z
    //for 6 jets
    //tt+Z  27    41
    if(fabs(evtweight-0.00506017)<1.0e-8)  return rate[18];  //eLpR  18
    if(fabs(evtweight-0.00232326)<1.0e-8) return 1.0-rate[18];     //eRpL    41
    break;
  case 6:  //ttbar + H   
    //for 6 jets
    //ttbar+H     29   42
    if(fabs(evtweight-0.00128075)<1.0e-7)  return rate[19];  //eLpR   32
    if(fabs(evtweight-0.000464348)<1.0e-8) return 1.0-rate[19];      //eRpL    41
    break;
  case 7:   //ZZ + H
    //for 6 jets
    //ZZ+H     30     43
    if(fabs(evtweight-0.00216898)<1.0e-8)  return rate[20];  //eLpR   19
    if(fabs(evtweight-0.00217191)<1.0e-8)  return rate[20];  //eLpR   19
    if(fabs(evtweight-0.0026375)<1.0e-8) return 1.0-rate[20];      //eRpL  43
    if(fabs(evtweight-0.00264003)<1.0e-8) return 1.0-rate[20];      //eRpL  43
    break;
  case 8:   //ZZZ
    //for 6 jets
    //ZZZ     32     44
    if(fabs(evtweight-0.0144831)<1.0e-8) return rate[21];   //eLpR   34
    if(fabs(evtweight-0.00492879)<1.0e-8) return 1.0-rate[21];      //eRpL
    break;
  }

  return 0.0;
}
