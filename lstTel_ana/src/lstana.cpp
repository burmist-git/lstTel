//my
#include "lstana.hh"

//root
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TVector3.h>

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>
#include <bits/stdc++.h>

using namespace std;

void lstana::Loop(TString histOut){
  //
  //
  _rnd = new TRandom3(1231231);
  //
  //
  Int_t i = 0;
  ///
  TH1D *h1_primPosX = new TH1D("h1_primPosX","h1_primPosX",2000,-100,100);
  TH1D *h1_primPosY = new TH1D("h1_primPosY","h1_primPosY",2000,-100,100);
  TH1D *h1_primPosZ = new TH1D("h1_primPosZ","h1_primPosZ",2000,500,1000);
  TH1D *h1_nPhot = new TH1D("h1_nPhot","h1_nPhot",1000,0,100000);
  TH1D *h1_nPhot_cut = new TH1D("h1_nPhot_cut","h1_nPhot_cut",1000,0,100000);
  ///
  TH1D *h1_trk_mirror_impact_point_X = new TH1D("h1_trk_mirror_impact_point_X","h1_trk_mirror_impact_point_X",2000,-100,100);
  TH1D *h1_trk_mirror_impact_point_Y = new TH1D("h1_trk_mirror_impact_point_Y","h1_trk_mirror_impact_point_Y",2000,-100,100);
  ///
  TH1D *h1_trk_theta_deg = new TH1D("h1_trk_theta_deg","h1_trk_theta_deg", 400, 0.0, 5.0);
  TH1D *h1_trk_phi_deg = new TH1D("h1_trk_phi_deg","h1_trk_phi_deg", 400, -200.0, 200.0);
  //
  TH1D *h1_Wavelength = new TH1D("h1_Wavelength","h1_Wavelength", 400, 250.0, 700.0);  
  //  
  TH1D *h1_trk_mom_GeV = new TH1D("h1_trk_mom_GeV","h1_trk_mom_GeV",400, 0.0, 250.0);
  //
  Double_t trk_mom_GeV;
  Double_t trk_theta;
  Double_t trk_phi;
  Double_t trk_theta_deg;
  Double_t trk_phi_deg;
  ///
  Double_t trk_mirror_impact_point_X;
  Double_t trk_mirror_impact_point_Y;
  Double_t trk_mirror_impact_point_Z;
  Double_t trk_mirror_impact_point_R;
  ///
  vector<TH2D*> h2_cam_v;
  vector<TH1D*> h1_phi_v;
  vector<TH1D*> h1_true_info_v;
  ///
  Int_t event_counter = 0;
  ///
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    if(jentry%1000 == 0)
      cout<<jentry<<endl;
    /////
    trk_mom_GeV = get_trk_mom_GeV();
    h1_trk_mom_GeV->Fill(trk_mom_GeV);
    /////
    h1_primPosX->Fill(primPosX/1000.0);
    h1_primPosY->Fill(primPosY/1000.0);
    h1_primPosZ->Fill(primPosZ/1000.0);
    h1_nPhot->Fill(nPhot);
    //
    get_trk_theta(trk_theta,trk_phi);
    trk_theta_deg = trk_theta*180.0/TMath::Pi();
    trk_phi_deg = trk_phi*180.0/TMath::Pi();
    //
    h1_trk_theta_deg->Fill(trk_theta_deg);
    h1_trk_phi_deg->Fill(trk_phi_deg);
    //
    for(Int_t i = 0;i<nPhot;i++){
      h1_Wavelength->Fill(Wavelength[i]);
    }
    //
    trk_mirror_impact_point_Z = -749.8; // m
    if(get_trk_mirror_impact_point(trk_mirror_impact_point_X, trk_mirror_impact_point_Y, trk_mirror_impact_point_Z)){
      ////
      trk_mirror_impact_point_R = TMath::Sqrt(trk_mirror_impact_point_X*trk_mirror_impact_point_X + trk_mirror_impact_point_Y*trk_mirror_impact_point_Y);
      h1_trk_mirror_impact_point_X->Fill(trk_mirror_impact_point_X);
      h1_trk_mirror_impact_point_Y->Fill(trk_mirror_impact_point_Y);
      //	  
      if(trk_theta_deg<1.2){
	if(trk_mirror_impact_point_R<10){
	  if(trk_mom_GeV>10){
	    h1_nPhot_cut->Fill(nPhot);
	    TString h2_tmp_name = "h2_cam_";
	    h2_tmp_name += event_counter;
	    h2_tmp_name += "_ev";
	    TString h1_phi_tmp_name = "h1_phi_";
	    h1_phi_tmp_name += event_counter;
	    h1_phi_tmp_name += "_ev";
	    TH2D *h2_tmp = new TH2D(h2_tmp_name.Data(),h2_tmp_name.Data(), 70,-1.2,1.2, 70,-1.2,1.2);
	    TH1D *h1_phi_tmp = new TH1D(h1_phi_tmp_name.Data(),h1_phi_tmp_name.Data(), 40,-200.0,200.0);
	    TString h1_true_tmp_name = "h1_true_";
	    h1_true_tmp_name += event_counter;
	    h1_true_tmp_name += "_ev";
	    TH1D *h1_true_tmp = new TH1D(h1_true_tmp_name.Data(),h1_true_tmp_name.Data(), 100,0.0,100);
	    ///
	    h1_true_tmp->SetBinContent(1,trk_mirror_impact_point_X);
	    h1_true_tmp->SetBinContent(2,trk_mirror_impact_point_Y);
	    h1_true_tmp->SetBinContent(3,trk_mirror_impact_point_Z);
	    h1_true_tmp->SetBinContent(4,trk_theta_deg);
	    h1_true_tmp->SetBinContent(5,trk_phi_deg);
	    ///
	    simulate_optical_systems(h2_tmp, h1_phi_tmp);
	    h2_cam_v.push_back(h2_tmp);
	    h1_phi_v.push_back(h1_phi_tmp);
	    h1_true_info_v.push_back(h1_true_tmp);
	    //
	    //save_to_csv(trk_mirror_impact_point_X,trk_mirror_impact_point_Y,trk_mirror_impact_point_Z,0.1,event_counter);
	  }
	}
      }
    }
    ////
    //
    event_counter++;    
    /////
  }
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //////
  TDirectory *camimage = rootFile->mkdir("camimage");
  camimage->cd();
  for(unsigned int ii = 0;ii<h2_cam_v.size();ii++){
    h2_cam_v.at(ii)->Write();
    h1_phi_v.at(ii)->Write();
    h1_true_info_v.at(ii)->Write();
  }
  rootFile->cd();
  //////
  h1_trk_mom_GeV->Write();
  h1_primPosX->Write();
  h1_primPosY->Write();
  h1_primPosZ->Write();
  h1_trk_mirror_impact_point_X->Write();
  h1_trk_mirror_impact_point_Y->Write();
  h1_nPhot->Write();
  h1_nPhot_cut->Write();
  //////
  h1_trk_theta_deg->Write();
  h1_trk_phi_deg->Write();
  //////
  h1_Wavelength->Write();
  //////
  rootFile->Close();
}

Double_t lstana::get_trk_mom_GeV(){
  TVector3 v(primMomX, primMomY, primMomZ);
  return v.Mag()/1000.0;
}

bool lstana::get_trk_mirror_impact_point(Double_t &trk_mirror_impact_point_X, Double_t &trk_mirror_impact_point_Y, Double_t trk_mirror_impact_point_Z){
  TVector3 trk_r0( primPosX/1000, primPosY/1000, primPosZ/1000);
  TVector3 trk_v( primMomX, primMomY, primMomZ);
  TVector3 plane_r0(0.0, 0.0, trk_mirror_impact_point_Z);
  TVector3 plane_n(0.0, 0.0, 1.0);
  Double_t div = trk_v.Dot(plane_n); 
  if(div == 0){
    if(plane_r0.Dot(plane_n) == trk_r0.Dot(plane_n)){ 
      trk_mirror_impact_point_X = primPosX;
      trk_mirror_impact_point_Y = primPosY;
    }
    else{
      trk_mirror_impact_point_X = -999.0;
      trk_mirror_impact_point_Y = -999.0;
    }
    return false;
  }
  Double_t t = (plane_r0.Dot(plane_n) - trk_r0.Dot(plane_n))/div;
  TVector3 trk_gen = trk_v*t + trk_r0;
  trk_mirror_impact_point_X = trk_gen.x();
  trk_mirror_impact_point_Y = trk_gen.y();
  //
  return true;
}

void lstana::get_trk_theta(Double_t &trk_theta, Double_t &trk_phi){
  TVector3 v( primMomX, primMomY, primMomZ);
  trk_theta = TMath::Pi() - v.Theta();
  //trk_phi = get_canonical_phi_from_Vphi(v.Phi());
  trk_phi = v.Phi();
}

Double_t lstana::get_canonical_phi_from_Vphi(Double_t phiv){
  if( phiv >= 0.0 && phiv <= TMath::Pi())
    return phiv;
  return TMath::Pi()-phiv;
}

void lstana::simulate_optical_systems(TH2D *h2, TH1D *h1_phi){
  Double_t phot_theta;
  Double_t phot_phi;
  Double_t phot_x_cam;
  Double_t phot_y_cam;
  Double_t phot_r_mirror;
  for(Int_t i = 0;i<nPhot;i++){
    if(Wavelength[i]>250.0 && Wavelength[i]<700.0){
      if(_rnd->Uniform()<=1.0){
	phot_r_mirror = TMath::Sqrt(PosX[i]*PosX[i] + PosY[i]*PosY[i]);
	//cout<<phot_r_mirror<<endl;
	if(phot_r_mirror/1000.0<_lst_Rmax_m){
	  TVector3 v(MomX[i],MomY[i],MomZ[i]);
	  phot_theta = TMath::Pi() - v.Theta();
	  phot_phi = v.Phi();
	  TVector2 v_cam;
	  v_cam.SetMagPhi(_lst_effective_focal_length_m*TMath::Tan(phot_theta), phot_phi);
	  phot_x_cam = v_cam.X();
	  phot_y_cam = v_cam.Y();
	  h2->Fill(phot_x_cam,phot_y_cam);
	  //cout<<phot_x_cam<<endl;
	  h1_phi->Fill(phot_phi*180.0/TMath::Pi());
	}
      }
    }
  }
}

void lstana::save_to_csv( Double_t trk_mir_impact_X, Double_t trk_mir_impact_Y, Double_t trk_mir_impact_Z,
			  Double_t save_percentage, Int_t event_counter){
  ofstream csvfile;
  TString csvname = "photon_info_"; 
  csvname += event_counter;
  csvname += "ev.csv";
  //
  csvfile.open(csvname.Data());
  csvfile<<"primMomX,"
	 <<"primMomY,"
	 <<"primMomZ,"
	 <<"primPosX,"
	 <<"primPosY,"
	 <<"primPosZ,"
	 <<"trk_mir_impact_X,"
	 <<"trk_mir_impact_Y,"
	 <<"trk_mir_impact_Z,"
	 <<"MomX,"
	 <<"MomY,"
	 <<"MomZ"<<endl;
  for(Int_t i = 0;i<nPhot;i++){
    if(Wavelength[i]>250.0 && Wavelength[i]<700.0){
      if(_rnd->Uniform()<save_percentage){
        csvfile<<primMomX<<","
	       <<primMomY<<","
	       <<primMomZ<<","
	       <<primPosX/1000<<","
	       <<primPosY/1000<<","
	       <<primPosZ/1000<<","
	       <<trk_mir_impact_X<<","
	       <<trk_mir_impact_Y<<","
	       <<trk_mir_impact_Z<<","
	       <<MomX[i]<<","
	       <<MomY[i]<<","
	       <<MomZ[i]<<endl;
      }
    }
  }
  //for(unsigned int i = 0; i <simp_hist_crop_v.size() ; i++){
  //for(unsigned int j = 0; j < _n_pixels; j++){
  //  if(j == (_n_pixels-1))
  //    csvfile<<simp_hist_crop_v.at(i)->GetBinContent(j+1);
  //  else
  //    csvfile<<simp_hist_crop_v.at(i)->GetBinContent(j+1)<<" ";
  //}
  //csvfile<<endl;
  //}
  csvfile.close();  
}
