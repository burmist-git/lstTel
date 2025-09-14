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

Int_t plots( TString fileN01, TString fileN02);

Int_t plots_mu_ring(Int_t key = 1){
  Int_t nn = 8;
  TString fileN01[nn];
  TString fileN02[nn];
  //
  fileN01[0] = "histSingle_lstTel_180theta_0phi_7mx_0my.root";
  fileN01[1] = "histSingle_lstTel_180theta_0phi_5mx_5my.root";
  fileN01[2] = "histSingle_lstTel_180theta_0phi_0mx_7my.root";
  fileN01[3] = "histSingle_lstTel_180theta_0phi_m5mx_5my.root";
  fileN01[4] = "histSingle_lstTel_180theta_0phi_m7mx_0my.root";
  fileN01[5] = "histSingle_lstTel_180theta_0phi_m5mx_m5my.root";
  fileN01[6] = "histSingle_lstTel_180theta_0phi_0mx_m7my.root";
  fileN01[7] = "histSingle_lstTel_180theta_0phi_5mx_m5my.root";
  //
  if(key == 1){
    for(Int_t i = 0;i<nn;i++){
      fileN02[i] = fileN01[i];
      plots( fileN01[i], fileN02[i]);
    }
    return 0;
  }
  return plots(  fileN01[1],  fileN02[1]);
}

Int_t plots( TString fileN01, TString fileN02){
  TFile *f01 = new TFile(fileN01.Data());
  TFile *f02 = new TFile(fileN02.Data());
  //
  TH2D *h2_PosY_vs_PosX = (TH2D*)f01->Get("h2_PosY_vs_PosX");
  TH2D *h2_PosY_vs_PosX_clean = (TH2D*)f01->Get("h2_PosY_vs_PosX_clean");
  TGraph *gr_recoring = (TGraph*)f01->Get("gr_recoring");
  TGraph *gr_center = (TGraph*)f01->Get("gr_center");
  TH1D *h1_phi_deg = (TH1D*)f01->Get("h1_phi_deg");
  TH1D *h1_phi_reco_deg = (TH1D*)f01->Get("h1_phi_predicted_deg");
  TH1D *h1_phi_reco_cam_deg = (TH1D*)f01->Get("h1_phi_predicted_camera_deg");
  //
  TH2D *h2_PosY_vs_PosX_02 = (TH2D*)f02->Get("h2_PosY_vs_PosX");
  TH2D *h2_PosY_vs_PosX_clean_02 = (TH2D*)f02->Get("h2_PosY_vs_PosX_clean");
  TGraph *gr_recoring_02 = (TGraph*)f02->Get("gr_recoring");
  TGraph *gr_center_02 = (TGraph*)f02->Get("gr_center");
  TH1D *h1_phi_deg_02 = (TH1D*)f02->Get("h1_phi_deg");
  //
  //
  //
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1500,500);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  c1->Divide(3,1);
  //
  /*
  c1->cd(1);
  h2_PosY_vs_PosX->SetMinimum(20);
  h2_PosY_vs_PosX->SetMaximum(250);
  h2_PosY_vs_PosX->Draw("ZCOLOR");
  h2_PosY_vs_PosX->GetXaxis()->SetTitle("x, mm");
  h2_PosY_vs_PosX->GetYaxis()->SetTitle("y, mm");
  c1->cd(2);
  h2_PosY_vs_PosX_clean->Draw("ZCOLOR");
  gr_recoring->Draw("P SAME");
  gr_center->SetMarkerStyle(43);
  gr_recoring->Draw("Pl SAME");
  gr_recoring->SetLineWidth(2.0);
  gr_recoring->SetLineColor(kMagenta + 2);
  gr_recoring->SetMarkerColor(kMagenta + 2);
  gr_center->Draw("P SAME");
  h2_PosY_vs_PosX_clean->GetXaxis()->SetTitle("x, mm");
  h2_PosY_vs_PosX_clean->GetYaxis()->SetTitle("y, mm");
  c1->cd(3);
  h1_phi_deg->SetLineWidth(2.0);
  h1_phi_deg->SetLineColor(kBlack);
  h1_phi_deg->Draw();
  h1_phi_deg->GetXaxis()->SetTitle("#phi, deg");
  //
  */
  c1->cd(1);
  h2_PosY_vs_PosX_02->SetMinimum(20);
  h2_PosY_vs_PosX_02->SetMaximum(250);
  h2_PosY_vs_PosX_02->Draw("ZCOLOR");
  h2_PosY_vs_PosX_02->GetXaxis()->SetTitle("x, mm");
  h2_PosY_vs_PosX_02->GetYaxis()->SetTitle("y, mm");
  c1->cd(2);
  h2_PosY_vs_PosX_clean_02->Draw("ZCOLOR");
  gr_recoring_02->Draw("P SAME");
  gr_center_02->SetMarkerStyle(43);
  gr_recoring_02->Draw("Pl SAME");
  gr_recoring_02->SetLineWidth(2.0);
  gr_recoring_02->SetLineColor(kMagenta + 2);
  gr_recoring_02->SetMarkerColor(kMagenta + 2);
  gr_center_02->Draw("P SAME");
  h2_PosY_vs_PosX_clean_02->GetXaxis()->SetTitle("x, mm");
  h2_PosY_vs_PosX_clean_02->GetYaxis()->SetTitle("y, mm");
  c1->cd(3);
  h1_phi_deg_02->SetLineWidth(2.0);
  h1_phi_deg_02->SetLineColor(kBlack);
  h1_phi_reco_deg->SetLineWidth(2.0);
  h1_phi_reco_deg->SetLineColor(kRed+2.0);
  h1_phi_reco_cam_deg->SetLineWidth(2.0);
  h1_phi_reco_cam_deg->SetLineColor(kMagenta+2.0);  
  //h1_phi_deg_02->SetMaximum(1000);
  h1_phi_deg_02->SetMaximum(300);
  h1_phi_deg_02->SetMinimum(0);
  h1_phi_deg_02->Draw();
  h1_phi_reco_deg->Draw("sames");
  h1_phi_reco_cam_deg->Draw("sames");
  h1_phi_deg_02->GetXaxis()->SetTitle("#phi, deg");

  
  TString fileout = fileN01;
  fileout += ".pdf";
  c1->SaveAs(fileout.Data());   
  
  /*
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
  */  

  return 0;
}
