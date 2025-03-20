//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

#include <time.h>

using namespace std;

TString getObjectName(TString namein, Int_t i);

Int_t plots(){

  TString fileN;
  TString objectName;
  //fileN = "./hist.root";
  fileN = "./histSingle.root";

  TFile *f1 = new TFile(fileN.Data());

  TH1D *h1_1 = (TH1D*)f1->Get("h1_dummyHisto");

  const Int_t nChannels = 10;
  TH1D *h1_Arr[nChannels];
  TGraph *gr_Arr[nChannels];
  Int_t i = 0;

  const Int_t colorArr[nChannels] = {    kBlack,     kRed,  kGreen + 3, kCyan - 3, kBlue + 1,
				        kGray+1, kMagenta, kYellow + 2,  kRed - 2, kBlue - 6};

  const Int_t markerArr[nChannels] = { 20, 20, 20, 20, 20,
				       21, 21, 21, 21, 21};

  for(i = 0;i<nChannels;i++){
    objectName=getObjectName("h1_dummyHistoArr",i);
    //cout<<objectName<<endl;
    h1_Arr[i] = (TH1D*)f1->Get(objectName.Data());
    objectName=getObjectName("gr_dummyGraph",i);
    gr_Arr[i] = (TGraph*)f1->Get(objectName.Data());
  }
  
  //h1_1->Rebin(4);
  //h1_1->SetTitle("");

  //TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,800,800);
  TCanvas *c1 = new TCanvas("c1",fileN.Data(),10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
 
  h1_1->SetLineColor(kBlack);
  h1_1->SetLineWidth(3.0);

  //h1_1->Draw();
  //h1_1->GetXaxis()->SetTitle("Value, Unit");
  //h1_1->GetXaxis()->SetRangeUser(-0.13,-0.12);

  /*
  for(i = 0;i<nChannels;i++){
    h1_Arr[i]->SetLineColor(colorArr[i]);
    h1_Arr[i]->SetLineWidth(3.0);
    if(i == 0){
      h1_Arr[i]->Draw();
      h1_Arr[i]->GetXaxis()->SetTitle("Value, Unit");
      //h1_1->GetXaxis()->SetRangeUser(-0.13,-0.12);
    }
    else
      h1_Arr[i]->Draw("same");
  }
  
  //h1_1->Fit("gaus");
  TString legInfo;
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  for(i = 0;i<nChannels;i++){
    legInfo = "ch ";legInfo += i;
    leg->AddEntry(h1_Arr[i], legInfo.Data(), "l");
  }
  leg->Draw();
  */

  TMultiGraph *mg = new TMultiGraph();
  
  for(i = 0;i<nChannels;i++){
    gr_Arr[i]->SetLineColor(colorArr[i]);
    gr_Arr[i]->SetLineWidth(3.0);
    gr_Arr[i]->SetMarkerColor(colorArr[i]);
    gr_Arr[i]->SetMarkerStyle(markerArr[i]);
    mg->Add(gr_Arr[i]);
  }

  mg->Draw("apl");
  
  mg->GetXaxis()->SetTitle("ValueX, Unit");
  mg->GetYaxis()->SetTitle("ValueY, Unit");
  
  TString legInfo;
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  for(i = 0;i<nChannels;i++){
    legInfo = "ch ";legInfo += i;
    leg->AddEntry(gr_Arr[i], legInfo.Data(), "apl");
  }
  leg->Draw();
  

  return 0;
}

TString getObjectName(TString namein, Int_t i){
  TString h1nameh = namein; h1nameh += "_"; h1nameh += "ch_"; h1nameh += i;
  return h1nameh;
}
