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

Int_t plots_sim_vs_data(){

  TString fileN01;
  TString fileN02;
  //
  fileN01 = "./hist.root";
  fileN02 = "/home/burmist/home2/work/CTA/muons_info/histOut_plot_df_short.root";
  
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());

  TH1D *h1_1 = (TH1D*)f01->Get("h1_ring_center_r_cut_d");
  TH1D *h1_2 = (TH1D*)f02->Get("h1_ring_center_r_cut_d");


  TCanvas *c1 = new TCanvas("c1","c1",10,10,600,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
 
  h1_1->SetLineColor(kBlack);
  h1_1->SetLineWidth(3.0);
  h1_2->SetLineColor(kRed);
  h1_2->SetLineWidth(3.0);

  h1_1->Draw();
  h1_2->Draw("same");
  
  //TString legInfo;
  TLegend *leg = new TLegend(0.6,0.6,0.9,0.9,"","brNDC");
  leg->AddEntry(h1_1, "Simulation", "apl");
  leg->AddEntry(h1_2, "Data", "apl");
  leg->Draw();
  

  return 0;
}
