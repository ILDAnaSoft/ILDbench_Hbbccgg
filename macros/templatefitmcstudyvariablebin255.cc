#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <string>
#include <sstream>
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TH2.h"
#include "TH3.h"
#include "TRandom.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooPlot.h"
#include "RooPoisson.h"
#include "RooHistPdf.h"
#include "RooMCStudy.h"
#include "RooExtendPdf.h"
#include "RooDataSet.h"
#include "RooHist.h"
#include "RooFitResult.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TF1.h"

using namespace RooFit;

string itos(int i);

int templatefitmcstudyvariablebin255(string pol){
  //gROOT->Reset();
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  //gStyle->SetEndErrorSize(0.0);

  const int cats = 9;
  string fname = "template_flavortag_" + pol + ".root";

  TFile * tf = new TFile(fname.c_str());
  TH3F *temp[cats];   //from rootfile
  TH3F *tempbb, *tempcc, *tempgg, *tempall, *tempbg;
  TH3F *temps[2];

  //histograms of parameters
  TH1F *pars[3];
  pars[0] = new TH1F("rbb","rbb",100,0.0,1.0e+6);
  pars[1] = new TH1F("rcc","rcc",100,0.0,1.0e+6);
  pars[2] = new TH1F("rgg","rgg",100,0.0,1.0e+6);

  string hname,hname2;
  for(int i=0;i<cats;i++){
    hname ="template_" + itos(i);
    hname2 ="template_" + itos(i) +"_2";
    temp[i] = (TH3F*) tf->Get(hname.c_str())->Clone(hname2.c_str());
    //temp[i]->Scale(500.0/1600.0);
  }

  //create variable bi n tenplates
  float xbin[3]={0.00, 0.50, 1.00};    //divide into 10 bins  //b-likeliness
  //float ybin[6]={0.00, 0.10, 0.20, 0.50, 0.90, 1.00};    //divide into 10 bins  //c-likeliness
  float ybin[6]={0.00, 0.30, 0.50, 0.80, 0.90, 1.00};    //divide into 10 bins  //c-likeliness
  float zbin[6]={0.00, 0.50, 0.70, 0.80, 0.90, 1.00};    //divide into 10 bins  //bc-likeliness

  //float xbin[11]={0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00};    //divide into 10 bins  //b-likeliness
  //float ybin[11]={0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00};    //divide into 10 bins  //b-likeliness
  //float zbin[11]={0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00};    //divide into 10 bins  //b-likeliness

  //float ybin[4]={0.00, 0.20, 0.90, 1.00};    //divide into 10 bins  //c-likeliness      
  //float zbin[4]={0.00, 0.70, 0.80, 1.00};    //divide into 10 bins  //bc-likeliness       

  const int nbinsx = 2;  
  const int nbins = 5;
  tempbb =  new TH3F("template_bb","",nbinsx,xbin,nbins,ybin,nbins,zbin);
  tempcc =  new TH3F("template_cc","",nbinsx,xbin,nbins,ybin,nbins,zbin);
  tempgg =  new TH3F("template_gg","",nbinsx,xbin,nbins,ybin,nbins,zbin);

  tempall = new TH3F("template_all","",nbinsx,xbin,nbins,ybin,nbins,zbin);
  tempbg = new TH3F("template_bg","",nbinsx,xbin,nbins,ybin,nbins,zbin);

  temps[0] = new TH3F("template_higgsbg","",nbinsx,xbin,nbins,ybin,nbins,zbin);
  temps[1] = new TH3F("template_smbg","",nbinsx,xbin,nbins,ybin,nbins,zbin);

  for(int cat=0;cat<cats;cat++){  //loop for templates
    for(int i1=1;i1<=nbinsx;i1++){
      for(int j1=1;j1<=nbins;j1++){
	for(int k1=1;k1<=nbins;k1++){
	  
	  double tmpval =0.0;
	  for(int i2=1;i2<=50;i2++){
	    for(int j2=1;j2<=50;j2++){
	      for(int k2=1;k2<=50;k2++){
		//cout << (i2-1)*0.02 << " " << i2*0.02 << endl;
		if((xbin[i1-1] <= (float)((i2-1)*0.02) && xbin[i1] > (float)((i2-1)*0.02))
		   && (ybin[j1-1] <= (float)((j2-1)*0.02) && ybin[j1] > (float)((j2-1)*0.02))
		   && (zbin[k1-1] <= (float)((k2-1)*0.02) && zbin[k1] > (float)((k2-1)*0.02)))
		  tmpval += temp[cat]->GetBinContent(i2,j2,k2);
	      }
	    }
	  }
	  
	  //signal
	  if(cat==0)
	    tempbb->SetBinContent(i1,j1,k1, tmpval);
	  //signal
	  if(cat==1)
	    tempcc->SetBinContent(i1,j1,k1, tmpval);
	  //signal
	  if(cat==2)
	    tempgg->SetBinContent(i1,j1,k1, tmpval);
	  if(cat>2)
	    tempbg->SetBinContent(i1,j1,k1,tempbg->GetBinContent(i1,j1,k1)
				  + tmpval);
	  tempall->SetBinContent(i1,j1,k1,tempall->GetBinContent(i1,j1,k1)
				 + tmpval);

	  if(cat>2 && cat<8)
	    temps[1]->SetBinContent(i1,j1,k1,temps[1]->GetBinContent(i1,j1,k1)
				 + tmpval);
	  if(cat==8)
	    temps[0]->SetBinContent(i1,j1,k1,temps[0]->GetBinContent(i1,j1,k1)
				 + tmpval);
	}
      }
    }
  }

  //crete 2D histo
  TH2F *temp2d[6];
  temp2d[0] = (TH2F*) tempbb->Project3D("yx")->Clone("bb2");
  temp2d[1] = (TH2F*) tempcc->Project3D("yx")->Clone("cc2");
  temp2d[2] = (TH2F*) tempgg->Project3D("yx")->Clone("gg2");
  temp2d[3] = (TH2F*) temps[0]->Project3D("yx")->Clone("o2");
  temp2d[4] = (TH2F*) temps[1]->Project3D("yx")->Clone("sm2");
  temp2d[5] = (TH2F*) tempall->Project3D("yx")->Clone("all2");
  for(int i=0;i<6;i++){
    temp2d[i]->SetTitle("");
    temp2d[i]->GetXaxis()->SetTitle("b-Likeliness");
    temp2d[i]->GetYaxis()->SetTitle("c-Likeliness");
  }

  //start pseudo-experiment
  //parameters
  double rbb=0.0, rcc=0.0, rgg=0.0;
  const int times = 10000;   
    
  //define roofit
  //x,y,z
  RooRealVar b("b","b",0.0,1.0);
  RooRealVar c("c","c",0.0,1.0);
  RooRealVar bc("bc","bc",0.0,1.0);
  
  //datahist
  RooDataHist dh("dh","data with b,c,bc",RooArgList(b,c,bc),tempall);
  RooDataHist bh("bh","bb with b,c,bc",RooArgList(b,c,bc),tempbb);
  RooDataHist ch("ch","cc with b,c,bc",RooArgList(b,c,bc),tempcc);
  RooDataHist gh("gh","gg with b,c,bc",RooArgList(b,c,bc),tempgg);
  RooDataHist bgh("bgh","background with b,c,bc",RooArgList(b,c,bc),tempbg);
  
  //histpdf
  RooHistPdf pdf0("pdfall","pdf of all the events",RooArgSet(b,c,bc),dh);
  RooHistPdf pdf1("pdfbb","pdf of bb",RooArgSet(b,c,bc),bh);
  RooHistPdf pdf2("pdfcc","pdf of cc",RooArgSet(b,c,bc),ch);
  RooHistPdf pdf3("pdfgg","pdf of gg",RooArgSet(b,c,bc),gh);
  RooHistPdf pdf4("pdfbg","pdf of bg",RooArgSet(b,c,bc),bgh);

  //extendpdf
  RooRealVar nbb("nbb","constant",tempbb->Integral());
  RooRealVar ncc("ncc","constant",tempcc->Integral());
  RooRealVar ngg("ngg","constant",tempgg->Integral());
  RooRealVar nbg("nbg","constant",tempbg->Integral());
  
  RooExtendPdf pdf21("pdf2bb", "scaled pdf", pdf1, nbb);
  RooExtendPdf pdf22("pdf2cc", "scaled pdf", pdf2, ncc);
  RooExtendPdf pdf23("pdf2gg", "scaled pdf", pdf3, ngg);
  RooExtendPdf pdf24("pdf2bg", "scaled pdf", pdf4, nbg);
  
  RooRealVar r0("r0","rbb",tempbb->Integral(),0.0,1.0e+7);
  RooRealVar r1("r1","rcc",tempcc->Integral(),0.0,1.0e+6);
  RooRealVar r2("r2","rgg",tempgg->Integral(),0.0,1.0e+6);
  RooRealVar r3("r3","constant",tempbg->Integral());
  
  RooAddPdf model("model","fitting model", RooArgList(pdf1,pdf2,pdf3,pdf4), RooArgList(r0,r1,r2,r3));
  //model.fitTo(dh);
  
    //under investigation
  RooMCStudy* mcs = new RooMCStudy(pdf0,RooArgSet(b,c,bc),FitModel(model),Binned(),Extended(kTRUE),  Silence(),
    				   FitOptions(Save(kTRUE),Extended(kTRUE),PrintEvalErrors(-1)));
  
  mcs->generateAndFit((int)(times*1.05),tempall->Integral());
  //mcs->fit(10,"test.dat");  //, 10);   //tempall->Integral());
  //mcs->calcPulls();

  //create hists
  TH1F *rs[3], *pulls[3];
  rs[0] = new TH1F("rs0","",100,0.55,1.45);
  rs[0]->SetLineColor(2);
  rs[0]->GetYaxis()->SetTitle("Number of Toy MC events");
  rs[0]->GetXaxis()->SetTitle("r_{bb}");
  rs[1] = new TH1F("rs1","",100,0.7,1.3);
  rs[1]->SetLineColor(3);
  rs[1]->GetYaxis()->SetTitle("Number of Toy MC events");
  rs[1]->GetXaxis()->SetTitle("r_{cc}");
  rs[2] = new TH1F("rs2","",100,0.0,2.0);
  rs[2]->SetLineColor(4);
  rs[2]->GetYaxis()->SetTitle("Number of Toy MC events");
  rs[2]->GetXaxis()->SetTitle("r_{gg}");
  pulls[0] = new TH1F("pulls0","",100,-5.0,5.0);
  pulls[0]->SetLineColor(2);
  pulls[0]->GetYaxis()->SetTitle("Number of Toy MC events");
  pulls[0]->GetXaxis()->SetTitle("Pull_{rbb}");
  pulls[1] = new TH1F("pulls1","",100,-5.0,5.0);
  pulls[1]->SetLineColor(3);
  pulls[1]->GetYaxis()->SetTitle("Number of Toy MC events");
  pulls[1]->GetXaxis()->SetTitle("Pull_{rcc}");
  pulls[2] = new TH1F("pulls2","",100,-5.0,5.0);
  pulls[2]->SetLineColor(4);
  pulls[2]->GetYaxis()->SetTitle("Number of Toy MC events");
  pulls[2]->GetXaxis()->SetTitle("Pull_{rgg}");
  
  //gaussian for fit
  TF1 *gg[3], *g0;
  gg[0] = new TF1("g0","gaus",0.5,+1.5);
  gg[0]->SetLineColor(1);
  gg[1] = new TF1("g1","gaus",0.0,+2.0);
  gg[1]->SetLineColor(1);
  gg[2] = new TF1("g2","gaus",0.0,+2.0);
  gg[2]->SetLineColor(1);
  g0 = new TF1("gg","gaus",-5.0,+5.0);
  g0->SetLineColor(1);

  for(int i=0;i<times;i++){
    rs[0]->Fill(mcs->fitParams(i)->getRealValue("r0") / tempbb->Integral());
    rs[1]->Fill(mcs->fitParams(i)->getRealValue("r1") / tempcc->Integral());
    rs[2]->Fill(mcs->fitParams(i)->getRealValue("r2") / tempgg->Integral());

    //pulls
    pulls[0]->Fill((mcs->fitParams(i)->getRealValue("r0") - tempbb->Integral())
		   /sqrt(mcs->fitResult(i)->covarianceMatrix()[0][0]));
    pulls[1]->Fill((mcs->fitParams(i)->getRealValue("r1") - tempcc->Integral())
		   /sqrt(mcs->fitResult(i)->covarianceMatrix()[1][1]));
    pulls[2]->Fill((mcs->fitParams(i)->getRealValue("r2") - tempgg->Integral())
		   /sqrt(mcs->fitResult(i)->covarianceMatrix()[2][2]));
  }

  string epsname, pdfname;

  //create visualize stuff                                                                      
  TPaveText *titletxt=new TPaveText(0.16,0.928,0.46,1.008,"NDC");
  titletxt->SetFillColor(0);
  titletxt->SetBorderSize(0);
  //titletxt->AddText("Higgs Self Coupling Analysis");
  titletxt->AddText("Higgs Branching Ratio Analysis");
  titletxt->AddText("e^{+}e^{-}#rightarrow#nu#nuH");
  //"ILD            "
  TPaveLabel *lumtxt=new TPaveLabel(0.40,0.930,1.0,0.985,"ILD        ","NDC");  //   #intL=9.4fb^{-1}","NDC");
  lumtxt->SetFillColor(0);
  lumtxt->SetBorderSize(0);
  lumtxt->SetTextFont(62);
  lumtxt->SetTextSize(1.1);
  TPaveLabel *lumtxt2=new TPaveLabel(0.67,0.925,0.82,0.97,"Preliminary","NDC");  //   #intL=9.4fb^{-1}","NDC");
  lumtxt2->SetFillColor(0);
  lumtxt2->SetBorderSize(0);
  lumtxt2->SetTextFont(62);
  lumtxt2->SetTextSize(0.8);

  TCanvas *c1 = new TCanvas("c1","Dynamic Filling Example",300,30,1600,1200);
  c1->GetFrame()->SetBorderSize(6);
  c1->GetFrame()->SetBorderMode(-1);
  c1->Divide(2,2);

  c1->cd(1);
  rs[0]->Draw("Hist");
  rs[0]->Fit("g0");
  gg[0]->Draw("same");
  lumtxt->Draw();
  titletxt->Draw();
 
  c1->cd(2);
  rs[1]->Draw("Hist");
  rs[1]->Fit("g1");
  gg[1]->Draw("same");
  lumtxt->Draw();
  titletxt->Draw();
 
  c1->cd(3);
  rs[2]->Draw("Hist");
  rs[2]->Fit("g2");
  gg[2]->Draw("same");
  lumtxt->Draw();
  titletxt->Draw();
 
  c1->cd(0);

  epsname = "plots/pseudoexperiment_" + pol + ".eps";
  pdfname = "plots/pseudoexperiment_" + pol + ".pdf";

  c1->SaveAs(epsname.c_str());
  c1->SaveAs(pdfname.c_str());

  TCanvas *c2 = new TCanvas("c2","Dynamic Filling Example",300,30,1600,1200);
  c2->GetFrame()->SetBorderSize(6);
  c2->GetFrame()->SetBorderMode(-1);
  c2->Divide(2,2);

  c2->cd(1);
  pulls[0]->Draw("Hist");
  pulls[0]->Fit("gg");
  g0->Draw("same");
  lumtxt->Draw();
  titletxt->Draw();
 
  c2->cd(2);
  pulls[1]->Draw("Hist");
  pulls[1]->Fit("gg");
  g0->Draw("same");
  lumtxt->Draw();
  titletxt->Draw();
 
  c2->cd(3);
  pulls[2]->Draw("Hist");
  pulls[2]->Fit("gg");
  g0->Draw("same");
  lumtxt->Draw();
  titletxt->Draw();

  epsname = "plots/pseudoexperiment_pull_" + pol + ".eps";
  pdfname = "plots/pseudoexperiment_pull_" + pol + ".pdf";
 
  c2->cd(0);
  c2->SaveAs(epsname.c_str());
  c2->SaveAs(pdfname.c_str());

  TCanvas *c3 = new TCanvas("c3","Dynamic Filling Example",300,30,2200,1200);
  c3->GetFrame()->SetBorderSize(6);
  c3->GetFrame()->SetBorderMode(-1);
  c3->Divide(3,2);

  //"ILD            "
  TPaveLabel *lumtxt3=new TPaveLabel(0.48,0.930,0.80,0.985,"ILD        ","NDC");  //   #intL=9.4fb^{-1}","NDC");
  lumtxt3->SetFillColor(0);
  lumtxt3->SetBorderSize(0);
  lumtxt3->SetTextFont(62);
  lumtxt3->SetTextSize(1.1);
  TPaveLabel *lumtxt4=new TPaveLabel(0.67,0.925,0.70,0.97,"Preliminary","NDC");  //   #intL=9.4fb^{-1}","NDC");
  lumtxt4->SetFillColor(0);
  lumtxt4->SetBorderSize(0);
  lumtxt4->SetTextFont(62);
  lumtxt4->SetTextSize(0.8);

  TPaveText *cc1=new TPaveText(0.53,0.84,0.78,0.78,"NDC");
  cc1->SetFillColor(0);
  cc1->SetBorderSize(0);
  //titletxt->AddText("Higgs Self Coupling Analysis");
  //cc1->AddText("For Visualization");
  //cc1->AddText("IDR-L, P(-0.8, +0.3)");
  cc1->AddText("H#rightarrowb#bar{b}");

  TPaveText *cc2=new TPaveText(0.53,0.84,0.78,0.78,"NDC");
  cc2->SetFillColor(0);
  cc2->SetBorderSize(0);
  //titletxt->AddText("Higgs Self Coupling Analysis");
  //cc2->AddText("For Visualization");
  //cc2->AddText("IDR-L, P(-0.8, +0.3)");
  cc2->AddText("H#rightarrowc#bar{c}");

  TPaveText *cc3=new TPaveText(0.53,0.84,0.78,0.78,"NDC");
  cc3->SetFillColor(0);
  cc3->SetBorderSize(0);
  //titletxt->AddText("Higgs Self Coupling Analysis");
  //cc3->AddText("For Visualization");
  //cc3->AddText("IDR-L, P(-0.8, +0.3)");
  cc3->AddText("H#rightarrowgg");

  TPaveText *cc4=new TPaveText(0.53,0.84,0.78,0.78,"NDC");
  cc4->SetFillColor(0);
  cc4->SetBorderSize(0);
  //titletxt->AddText("Higgs Self Coupling Analysis");
  //cc4->AddText("For Visualization");
  //cc4->AddText("IDR-L, P(-0.8, +0.3)");
  cc4->AddText("Other Higgs");

  TPaveText *cc5=new TPaveText(0.53,0.84,0.78,0.78,"NDC");
  cc5->SetFillColor(0);
  cc5->SetBorderSize(0);
  //titletxt->AddText("Higgs Self Coupling Analysis");
  //cc5->AddText("For Visualization");
  //cc5->AddText("IDR-L, P(-0.8, +0.3)");
  cc5->AddText("SM background");

  TPaveText *cc6=new TPaveText(0.53,0.84,0.78,0.78,"NDC");
  cc6->SetFillColor(0);
  cc6->SetBorderSize(0);
  //titletxt->AddText("Higgs Self Coupling Analysis");
  //cc6->AddText("For Visualization");
  //cc6->AddText("IDR-L, P(-0.8, +0.3)");
  cc6->AddText("Data");

  //common texts
  TPaveText *cc7=new TPaveText(0.40,0.06,1.00,0.0,"NDC");
  cc7->SetFillColor(0);
  cc7->SetBorderSize(0);
  cc7->AddText("IDR-L, P(-0.8, +0.3)");
  //cc7->AddText("IDR-L, P(-0.8, +0.3)");
  //cc7->AddText("Data");

  c3->cd(1);
  gPad->SetLogz(0);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.20);
  temp2d[0]->GetYaxis()->SetTitleOffset(1.0);
  temp2d[0]->GetXaxis()->SetTitleOffset(0.8);
  temp2d[0]->SetLineWidth(1.0);
  temp2d[0]->Draw("colz");
  lumtxt3->Draw();
  lumtxt4->Draw();
  titletxt->Draw();
  //ilcleg->Draw();
  cc1->Draw();

  c3->cd(2);
  gPad->SetLogz(0);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.20);
  temp2d[1]->GetYaxis()->SetTitleOffset(1.0);
  temp2d[1]->GetXaxis()->SetTitleOffset(0.8);
  temp2d[1]->SetLineWidth(1.0);
  temp2d[1]->Draw("colz");
  lumtxt3->Draw();
  lumtxt4->Draw();
  titletxt->Draw();
  //ilcleg->Draw();
  cc2->Draw();

  c3->cd(3);
  gPad->SetLogz(0);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.20);
  temp2d[2]->GetYaxis()->SetTitleOffset(1.0);
  temp2d[2]->GetXaxis()->SetTitleOffset(0.8);
  temp2d[2]->SetLineWidth(1.0);
  temp2d[2]->Draw("colz");
  lumtxt3->Draw();
  lumtxt4->Draw();
  titletxt->Draw();
  //ilcleg->Draw();
  cc3->Draw();

  c3->cd(4);
  gPad->SetLogz(0);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.20);
  temp2d[3]->GetYaxis()->SetTitleOffset(1.0);
  temp2d[3]->GetXaxis()->SetTitleOffset(0.8);
  temp2d[3]->SetLineWidth(1.0);
  temp2d[3]->Draw("colz");
  lumtxt3->Draw();
  lumtxt4->Draw();
  titletxt->Draw();
  //ilcleg->Draw();
  cc4->Draw();

  c3->cd(5);
  gPad->SetLogz(0);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.20);
  temp2d[4]->GetYaxis()->SetTitleOffset(1.0);
  temp2d[4]->GetXaxis()->SetTitleOffset(0.8);
  temp2d[4]->SetLineWidth(1.0);
  temp2d[4]->Draw("colz");
  lumtxt3->Draw();
  lumtxt4->Draw();
  titletxt->Draw();
  //ilcleg->Draw();
  cc5->Draw();

  c3->cd(6);
  gPad->SetLogz(0);
  gPad->SetLeftMargin(0.20);
  gPad->SetRightMargin(0.20);
  temp2d[5]->GetYaxis()->SetTitleOffset(1.0);
  temp2d[5]->GetXaxis()->SetTitleOffset(0.8);
  temp2d[5]->SetLineWidth(1.0);
  temp2d[5]->Draw("colz");
  lumtxt3->Draw();
  lumtxt4->Draw();
  titletxt->Draw();
  //ilcleg->Draw();
  cc6->Draw();
  cc7->Draw();

  c3->cd(0);
  c3->SaveAs("./template_variablebin.eps");
  c3->SaveAs("./template_variablebin.eps");


  TFile *tout = new TFile("./pseudoexperiment.root","RECREATE");
  rs[0]->Write();  
  rs[1]->Write();  
  rs[2]->Write();  
  pulls[0]->Write();
  pulls[1]->Write();
  pulls[2]->Write();
  tout->Close();

  return 0;
}

string itos(int i)
{
  stringstream s;
  s << i;
  return s.str();
}
