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

void lstana::Loop01(TString histOut){
  ///
  cout<<"lstana::Loop01"<<endl;
  ///
  _rnd = new TRandom3(1231231);
  ///
  TH1D *h1_trk_mom_GeV = new TH1D("h1_trk_mom_GeV","h1_trk_mom_GeV",400, 0.0, 250.0);
  TH1D *h1_nPhot = new TH1D("h1_nPhot","h1_nPhot",1000,0,100000);
  TH1D *h1_PosX = new TH1D("h1_PosX","h1_PosX",1000,-2000,2000);
  TH1D *h1_PosY = new TH1D("h1_PosY","h1_PosY",1000,-2000,2000);
  TH1D *h1_PosZ = new TH1D("h1_PosZ","h1_PosZ",1000,-3000000,3000000);
  TH1D *h1_ParentID = new TH1D("h1_ParentID","h1_ParentID",200,0,100);
  TH1D *h1_photPathLen = new TH1D("h1_photPathLen","h1_photPathLen",10000,0,10000000);
  //
  //TH2D *h2_PosY_vs_PosX = new TH2D("h2_PosY_vs_PosX","h2_PosY_vs_PosX",300,-2000,2000,300,-2000,2000);
  TH2D *h2_PosY_vs_PosX = new TH2D("h2_PosY_vs_PosX","h2_PosY_vs_PosX",100,-1200,1200,100,-1200,1200);
  TH2D *h2_PosY_vs_PosX_clean = new TH2D("h2_PosY_vs_PosX_clean","h2_PosY_vs_PosX_clean",100,-1200,1200,100,-1200,1200);
  //
  //
  //
  TH2D *h2_PosY_vs_PosX_arr[nChannels];
  h2D2Init(h2_PosY_vs_PosX_arr, "h2_PosY_vs_PosX", "h2_PosY_vs_PosX",
	   100,-1200,1200,100,-1200,1200);
   
  
  //
  TH1D *h1_phi_predicted_deg = new TH1D("h1_phi_predicted_deg","h1_phi_predicted_deg",32,-12.0, 372.0);
  TH1D *h1_phi_predicted_camera_deg = new TH1D("h1_phi_predicted_camera_deg","h1_phi_predicted_camera_deg",32,-12.0, 372.0);
  //
  Long64_t nentries = fChain->GetEntriesFast();
  cout<<"nentries = "<<nentries<<endl;
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  //for (Long64_t jentry=0; jentry<1;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //if(jentry%1000 == 0)
    //cout<<jentry<<endl;
    cout<<"jentry = "<<jentry<<endl;
    /////
    //
    h1_trk_mom_GeV->Fill(get_trk_mom_GeV());
    h1_nPhot->Fill(nPhot);
    //
    //
    for(Int_t i = 0;i<nPhot;i++){
      h1_ParentID->Fill(ParentID[i]);
      h1_photPathLen->Fill(photPathLen[i]);
      if(TMath::Sqrt(PosX[i]*PosX[i] + PosY[i]*PosY[i])<1500.0){
	h1_PosX->Fill(PosX[i]);
	h1_PosY->Fill(PosY[i]);
	h1_PosZ->Fill(PosZ[i]);
	h2_PosY_vs_PosX->Fill(PosX[i], PosY[i]);
	if((Int_t)jentry<nChannels)
	  h2_PosY_vs_PosX_arr[(Int_t)jentry]->Fill(PosX[i], PosY[i]);
      }
    }
    //
    //
    //
    //
  }
  //
  //
  TGraph *gr_data = new TGraph();
  TGraph *gr_weight = new TGraph();
  TGraph *gr_center = new TGraph();
  TGraph *gr_test = new TGraph();
  gr_data->SetNameTitle("gr_data","gr_data");
  gr_weight->SetNameTitle("gr_weight","gr_weight");
  gr_center->SetNameTitle("gr_center","gr_center");
  gr_test->SetNameTitle("gr_test","gr_test");
  //
  //
  get_weight_data_from_TH2D( h2_PosY_vs_PosX, gr_data, gr_weight, 0);
  Double_t cx_reco,  cy_reco,  r_reco;
  get_optimum_circular_fit_Chaudhuri(gr_data, gr_weight, cx_reco,  cy_reco,  r_reco);
  TH1D *h1_ring_radius = new TH1D("h1_ring_radius","h1_ring_radius",1000, 0.0, 1200);
  TH1D *h1_delta_ring_radius = new TH1D("h1_delta_ring_radius","h1_delta_ring_radius",1000, -1200, 1200);
  //get_ring_radius(h1_ring_radius, h1_delta_ring_radius, gr_data, gr_weight, cx_reco,  cy_reco, r_reco);
  get_partial_ring_radius(h1_ring_radius, h1_delta_ring_radius, gr_test, gr_data, gr_weight, cx_reco,  cy_reco, r_reco);
  
  cout<<"cx_reco "<<cx_reco<<endl
      <<"cy_reco "<<cy_reco<<endl
      <<"r_reco  "<<r_reco<<endl;
  //
  TGraph *gr_recoring = new TGraph();
  gr_recoring->SetNameTitle("gr_recoring","gr_recoring");
  gen_ring(gr_recoring, 1000, cx_reco, cy_reco, r_reco);
  gr_center->SetPoint(0,cx_reco, cy_reco);
  //
  //
  //TH1D *h1_phi_deg = new TH1D("h1_phi_deg","h1_phi_deg",380, -10.0, 370.0);
  TH1D *h1_phi_deg = new TH1D("h1_phi_deg","h1_phi_deg",32,-12.0, 372.0);
  get_phi_dist_and_clean_ring(h1_phi_deg, h2_PosY_vs_PosX_clean, (Float_t)cx_reco, (Float_t)cy_reco, (Float_t)r_reco, (Float_t)r_reco*0.1);
  predict_phi_dist(h1_phi_predicted_deg, h1_PosZ, 10);
  predict_phi_dist_camera(h1_phi_predicted_camera_deg, h1_PosZ, 10);
  //
  //
  //
  //
  //
  //
  TFile* rootFile = new TFile(histOut.Data(), "RECREATE", " Histograms", 1);
  rootFile->cd();
  if (rootFile->IsZombie()){
    cout<<"  ERROR ---> file "<<histOut.Data()<<" is zombi"<<endl;
    assert(0);
  }
  else
    cout<<"  Output Histos file ---> "<<histOut.Data()<<endl;
  //////
  rootFile->cd();
  h1_trk_mom_GeV->Write();
  h1_nPhot->Write();
  h1_PosX->Write();
  h1_PosY->Write();
  h1_PosZ->Write();
  h2_PosY_vs_PosX->Write();
  h2_PosY_vs_PosX_clean->Write();
  h1_ParentID->Write();
  h1_photPathLen->Write();
  gr_data->Write();
  gr_weight->Write();
  gr_recoring->Write();
  gr_center->Write();
  h1_phi_deg->Write();
  h1_phi_predicted_deg->Write();
  h1_phi_predicted_camera_deg->Write();
  //
  h1_ring_radius->Write();
  h1_delta_ring_radius->Write();
  gr_test->Write();
  //
  for(Int_t ii = 0; ii<nChannels; ii++)
    h2_PosY_vs_PosX_arr[ii]->Write();
  //
  //
  rootFile->Close();  
}

void lstana::Loop(TString histOut){
  //
  //
  _rnd = new TRandom3(1231231);
  //
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
  TH1D *h1_ringCompleteness = new TH1D("h1_ringCompleteness","h1_ringCompleteness",50, 0.0, 1.1);
  //
  //
  TH1D *h1_ring_center_r_cut = new TH1D("h1_ring_center_r_cut","h1_ring_center_r_cut",50, 0.0, 3.0);
  TH1D *h1_ring_center_r_cut_w = new TH1D("h1_ring_center_r_cut_w","h1_ring_center_r_cut_w",50, 0.0, 3.0);
  TH1D *h1_ring_center_r_cut_d = new TH1D("h1_ring_center_r_cut_d","h1_ring_center_r_cut_d",50, 0.0, 3.0);
  //
  //
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
  Double_t ringCompleteness;
  ///
  vector<TH2D*> h2_cam_v;
  vector<TH1D*> h1_phi_v;
  vector<TH1D*> h1_true_info_v;
  ///
  Int_t event_counter = 0;
  ///
  TH2D *h2_impact_r_array_vs_ring_center_r = new TH2D("h2_impact_r_array_vs_ring_center_r","h2_impact_r_array_vs_ring_center_r",30, 0.0, 3.0,60, 0.0, 24.0);
  TH2D *h2_impact_r_array_vs_ring_center_r_w = new TH2D("h2_impact_r_array_vs_ring_center_r_w","h2_impact_r_array_vs_ring_center_r_w",30, 0.0, 3.0,60, 0.0, 24.0);
  TH2D *h2_impact_r_array_vs_ring_center_r_d = new TH2D("h2_impact_r_array_vs_ring_center_r_d","h2_impact_r_array_vs_ring_center_r_d",30, 0.0, 3.0,60, 0.0, 24.0);    
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
    trk_mirror_impact_point_Z = -749.8; // m
    get_trk_mirror_impact_point(trk_mirror_impact_point_X, trk_mirror_impact_point_Y, trk_mirror_impact_point_Z);
    trk_mirror_impact_point_R = TMath::Sqrt(trk_mirror_impact_point_X*trk_mirror_impact_point_X + trk_mirror_impact_point_Y*trk_mirror_impact_point_Y);
    //if(trk_theta_deg<0.8){
    //if(trk_mirror_impact_point_R<8){
    //if(trk_mom_GeV>12){
    if(trk_theta_deg<10.0){
      if(trk_mirror_impact_point_R<100){
	if(trk_mom_GeV>10){
	  //
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
	  for(Int_t i = 0;i<nPhot;i++)
	    h1_Wavelength->Fill(Wavelength[i]);
	  //
	  h1_trk_mirror_impact_point_X->Fill(trk_mirror_impact_point_X);
	  h1_trk_mirror_impact_point_Y->Fill(trk_mirror_impact_point_Y);
	  //
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
	  ringCompleteness = get_ring_Completeness(h1_phi_tmp, 50.0);
	  h1_ringCompleteness->Fill(ringCompleteness);
	  h2_cam_v.push_back(h2_tmp);
	  h1_phi_v.push_back(h1_phi_tmp);
	  h1_true_info_v.push_back(h1_true_tmp);
	  //	  
	  h2_impact_r_array_vs_ring_center_r->Fill(trk_theta_deg,trk_mirror_impact_point_R );
	  h2_impact_r_array_vs_ring_center_r_w->Fill(trk_theta_deg,trk_mirror_impact_point_R, ringCompleteness);
	  //
	  if(trk_mirror_impact_point_R<5){
	    h1_ring_center_r_cut->Fill( trk_theta_deg);
	    h1_ring_center_r_cut_w->Fill( trk_theta_deg, ringCompleteness);
	  }
	  //
	  //save_to_csv(trk_mirror_impact_point_X,trk_mirror_impact_point_Y,trk_mirror_impact_point_Z,0.1,event_counter);
	}
      }
    }
    ////
    //
    event_counter++;    
    /////
  }
  //
  TH2D_divide(h2_impact_r_array_vs_ring_center_r_w,h2_impact_r_array_vs_ring_center_r,h2_impact_r_array_vs_ring_center_r_d);
  TH1D_divide(h1_ring_center_r_cut_w,h1_ring_center_r_cut,h1_ring_center_r_cut_d);
  //
  h2_impact_r_array_vs_ring_center_r_d->SetMinimum(0.0);
  h2_impact_r_array_vs_ring_center_r_d->SetMaximum(1.0);
  //
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
  h1_ringCompleteness->Write();
  //
  h2_impact_r_array_vs_ring_center_r->Write();
  h2_impact_r_array_vs_ring_center_r_w->Write();
  h2_impact_r_array_vs_ring_center_r_d->Write();
  //
  //
  h1_ring_center_r_cut_w->Write();
  h1_ring_center_r_cut->Write();
  h1_ring_center_r_cut_d->Write();
  //
  //
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

Double_t lstana::get_ring_Completeness(TH1D *h1, Double_t threshold){
  Double_t counter_yes = 0.0;
  Double_t norm = 0.0;
  for(Int_t i = 1;i<=h1->GetNbinsX();i++){
    if(h1->GetBinCenter(i)>=-180.0 && h1->GetBinCenter(i)<=180.0)
      norm += 1.0;
    if(h1->GetBinContent(i)>=threshold)
      counter_yes += 1.0;
  }
  if(norm > 0.0)
    return counter_yes/norm;
  else
    return -999.0;
}

void lstana::TH2D_divide( TH2D *h2_w, TH2D *h2, TH2D *h2_norm){
  Double_t val;
  Double_t norm;
  Double_t val_norm;
  for(Int_t i = 1;i<=h2_w->GetNbinsX();i++){
    for(Int_t j = 1;j<=h2_w->GetNbinsY();j++){
      val = h2_w->GetBinContent(i,j);
      norm = h2->GetBinContent(i,j);
      if(norm>0)
        val_norm = val/norm;
      else
        val_norm = 0.0;
      h2_norm->SetBinContent(i,j,val_norm);
    }
  }
}

void lstana::TH1D_divide( TH1D *h1_w, TH1D *h1, TH1D *h1_norm){
  Double_t val;
  Double_t norm;
  Double_t val_norm;
  for(Int_t i = 1;i<=h1_w->GetNbinsX();i++){
    val = h1_w->GetBinContent(i);
    norm = h1->GetBinContent(i);
    if(norm>0)
      val_norm = val/norm;
    else
      val_norm = 0.0;
    h1_norm->SetBinContent(i,val_norm);
  }
}

void lstana::get_optimum_circular_fit_Chaudhuri(TGraph *gr_data, TGraph *gr_weight, Double_t &cx_reco,  Double_t &cy_reco,  Double_t &r_reco){
  Double_t x, y;
  Double_t M_tot = 0.0;
  Double_t M_w_tot = 0.0;
  Double_t x_mean = 0.0;
  Double_t y_mean = 0.0;
  Double_t w = 0.0;
  Double_t A = 0.0;
  Double_t A_ = 0.0;
  Double_t B = 0.0;
  Double_t B_ = 0.0;
  Double_t C = 0.0;
  Double_t C_ = 0.0;
  for(Int_t i = 0; i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    gr_weight->GetPoint(i,w,w);
    x_mean += x*w;
    y_mean += y*w;
    M_w_tot += w;
    M_tot++;
  }
  x_mean /= M_w_tot;
  y_mean /= M_w_tot; 
  //
  for(Int_t i = 0; i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    gr_weight->GetPoint(i,w,w);
    A  += w*(x - x_mean)*x;
    A_ += w*(y - y_mean)*x;
    B  += w*(x - x_mean)*y;
    B_ += w*(y - y_mean)*y;
    C  += w*(x - x_mean)*(x*x + y*y);
    C_ += w*(y - y_mean)*(x*x + y*y);
  }
  C  /= 2;
  C_ /= 2;
  //
  cx_reco = (B_*C - B*C_) / (A*B_ - A_*B);
  cy_reco = (A_*C - A*C_) / (A_*B - A*B_);
  //
  r_reco = 0.0;
  for(Int_t i = 0; i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    gr_weight->GetPoint(i,w,w);
    r_reco += w*((x - cx_reco)*(x - cx_reco) + (y - cy_reco)*(y - cy_reco));
  }  
  r_reco /= M_w_tot;
  r_reco = TMath::Sqrt(r_reco);
  //cx_reco = x_mean;
  //cy_reco = y_mean;
  //gr_data->GetPoint(0,x,y);
  //r_reco = TMath::Sqrt((x - cx_reco)*(x - cx_reco) + (y - cy_reco)*(y - cy_reco));    
}

void lstana::get_weight_data_from_TH2D( TH2D *h2, TGraph *gr_data, TGraph *gr_weight, Double_t treshold){
  Double_t x;
  Double_t y;
  Double_t w;
  for(Int_t i = 1; i <= h2->GetNbinsX(); i++){
    for(Int_t j = 1; j<= h2->GetNbinsY(); j++){
      x = h2->GetXaxis()->GetBinCenter(i);
      y = h2->GetYaxis()->GetBinCenter(j);
      w = h2->GetBinContent(i,j);
      if(w>treshold){
	gr_data->SetPoint(gr_data->GetN(),x,y);
	gr_weight->SetPoint(gr_weight->GetN(),w,w);
      }
      //cout<<"x "<<x<<endl
      //  <<"y "<<y<<endl
      //  <<"w "<<w<<endl;
    }
  }
  //cout<<"h2->GetNbinsY() "<<h2->GetNbinsY()<<endl;
}

void lstana::gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R){
  Double_t phi = 0.0;
  TVector2 rc(x0,y0);
  for(Int_t i = 0;i<np;i++){
    TVector2 p;
    p.SetMagPhi(R,2*TMath::Pi()/(np-1)*i);
    TVector2 pt = rc + p;
    gr->SetPoint( i, pt.X(), pt.Y());
  }
}

void lstana::get_phi_dist_and_clean_ring(TH1D *h1_deg, TH2D *h2_PosY_vs_PosX_clean, Float_t cx_reco, Float_t cy_reco, Float_t r_reco, Float_t d_r_reco){
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  //for (Long64_t jentry=0; jentry<nentries;jentry++) {
  for (Long64_t jentry=0; jentry<1;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //if(jentry%1000 == 0)
    //cout<<jentry<<endl;
    //
    for(Int_t i = 0;i<nPhot;i++){
      if(_rnd->Uniform()<=0.2){
	TVector2 v( PosX[i] - cx_reco, PosY[i] - cy_reco);
	if( (v.Mod()>r_reco - d_r_reco) && (v.Mod()<r_reco + d_r_reco) ){
	  h1_deg->Fill(v.Phi()*180.0/TMath::Pi());
	  h2_PosY_vs_PosX_clean->Fill(PosX[i],PosY[i]);
	}
      }
    }
  }
}

void lstana::predict_phi_dist_camera(TH1D *h1_phi_deg, TH1D *h1_PosZ, Double_t n_norm){
  Double_t parabolicMirror_f  = 28.0; //m
  Double_t sensitive_sizeX = 3.140;   //m
  Double_t sensitive_sizeY = 3.140;   //m
  Double_t sensitive_eff_R = TMath::Sqrt(sensitive_sizeX*sensitive_sizeY/TMath::Pi());   //m
  //cout<<sensitive_eff_R<<endl;
  Double_t sensitive_dZ0 = 1.0/(1.0/parabolicMirror_f - 1.0/(12000)); //m
  Double_t photon_camera_pos_m = h1_PosZ->GetMean()/1000.0;
  Double_t parabolicMirror_R = 12.393;
  Double_t phi;
  //
  //
  //cout<<"parabolicMirror_f   "<<parabolicMirror_f<<endl
  //    <<"sensitive_sizeX     "<<sensitive_sizeX<<endl
  //   <<"sensitive_sizeY     "<<sensitive_sizeY<<endl
  //   <<"sensitive_dZ0       "<<sensitive_dZ0<<endl
  //   <<"photon_camera_pos_m "<<photon_camera_pos_m<<endl;
  //
  //
  //
  Double_t impact_x_m;
  Double_t impact_y_m;
  Double_t impact_z_m;
  //
  reco_muon_impact_point(photon_camera_pos_m,sensitive_dZ0, impact_x_m, impact_y_m, impact_z_m);
  for(Int_t i = 1;i<=h1_phi_deg->GetNbinsX();i++){
    phi = h1_phi_deg->GetBinCenter(i)/180.0*TMath::Pi();
    h1_phi_deg->SetBinContent(i,d_zero( impact_x_m, impact_y_m, sensitive_eff_R, n_norm, phi));
  }
}


void lstana::predict_phi_dist(TH1D *h1_phi_deg, TH1D *h1_PosZ, Double_t n_norm){
  Double_t parabolicMirror_f  = 28.0; //m
  Double_t sensitive_sizeX = 3.140;   //m
  Double_t sensitive_sizeY = 3.140;   //m
  Double_t sensitive_dZ0 = 1.0/(1.0/parabolicMirror_f - 1.0/(12000)); //m
  Double_t photon_camera_pos_m = h1_PosZ->GetMean()/1000.0;
  Double_t parabolicMirror_R = 12.393;
  Double_t phi;
  //
  //
  //cout<<"parabolicMirror_f   "<<parabolicMirror_f<<endl
  //    <<"sensitive_sizeX     "<<sensitive_sizeX<<endl
  //   <<"sensitive_sizeY     "<<sensitive_sizeY<<endl
  //   <<"sensitive_dZ0       "<<sensitive_dZ0<<endl
  //   <<"photon_camera_pos_m "<<photon_camera_pos_m<<endl;
  //
  //
  //
  Double_t impact_x_m;
  Double_t impact_y_m;
  Double_t impact_z_m;
  //
  reco_muon_impact_point(photon_camera_pos_m,sensitive_dZ0, impact_x_m, impact_y_m, impact_z_m);
  for(Int_t i = 1;i<=h1_phi_deg->GetNbinsX();i++){
    phi = h1_phi_deg->GetBinCenter(i)/180.0*TMath::Pi();
    h1_phi_deg->SetBinContent(i,d_zero( impact_x_m, impact_y_m, parabolicMirror_R, n_norm, phi));
  }
  //
  cout<<"lstana primPosX_m "<<primPosX/1000.0<<endl
      <<"lstana primPosY_m "<<primPosY/1000.0<<endl
      <<"lstana primPosZ_m "<<primPosZ/1000.0<<endl;
  //
  //
  cout<<"lstana impact_x_m "<<impact_x_m<<endl
      <<"lstana impact_y_m "<<impact_y_m<<endl
      <<"lstana impact_z_m "<<impact_z_m<<endl;
}

void lstana::reco_muon_impact_point( Double_t phot_z_hit, Double_t sensitive_dZ, Double_t &impact_x, Double_t &impact_y, Double_t &impact_z){
  //
  impact_z = phot_z_hit - sensitive_dZ;
  //
  Double_t dz = primPosZ/1000.0 - impact_z;
  //cout<<"dz "<<dz<<endl;
  TVector3 mu_v(primMomX, primMomY, primMomZ);
  //cout<<"primMomX "<<primMomX<<endl
  //    <<"primMomY "<<primMomY<<endl
  //    <<"primMomZ "<<primMomZ<<endl;
  TVector3 mu_gen(primPosX/1000.0, primPosY/1000.0, primPosZ/1000.0);
  //TVector3 mu_int = (-dz/primMomZ)*mu_v + mu_gen;
  //
  //cout<<"mu_v.Theta() "<<mu_v.Theta()*180.0/TMath::Pi()<<endl
  //  <<"mu_v.Phi()   "<<mu_v.Phi()*180.0/TMath::Pi()<<endl;
  //
  impact_x = dz*TMath::Tan(TMath::Pi() - mu_v.Theta())*TMath::Cos(mu_v.Phi()) + mu_gen.X();
  impact_y = dz*TMath::Tan(TMath::Pi() - mu_v.Theta())*TMath::Sin(mu_v.Phi()) + mu_gen.Y();
  //impact_z = mu_int.Z();
}

Double_t lstana::d_zero( Double_t impact_x, Double_t impact_y, Double_t R_mirror, Double_t norm, Double_t phi){
  TVector2 v(impact_x, impact_y);
  Double_t rho = v.Mod();
  Double_t phi0 = v.Phi() + TMath::Pi();
  Double_t rho_R = rho/R_mirror;
  Double_t d = 0.0;
  Double_t sqrt_sin_2 = 0.0;
  Double_t phi_max;
  if(rho_R > 1.0){
    sqrt_sin_2 = rho_R*rho_R*TMath::Power(TMath::Sin(phi-phi0),2.0);
    phi_max = TMath::ASin(1.0/rho_R);
    if(sqrt_sin_2<=1.0){
      if(TMath::Abs(phi-phi0) <= phi_max || phi<=phi_max){
	d = 2.0*norm*R_mirror*TMath::Sqrt(1.0 - sqrt_sin_2);
      }
      else{
	d = 0.0;
      }
    }
    else{
      d = 0.0;
    }
  }
  else{
    sqrt_sin_2 = rho_R*rho_R*TMath::Power(TMath::Sin(phi-phi0),2.0);
    if(sqrt_sin_2<=1.0){
      d = norm*R_mirror*TMath::Sqrt(1.0 - sqrt_sin_2);
      d += norm*R_mirror*rho_R*TMath::Cos(phi-phi0);
    }
    else{
      d = 0.0;
    }
  }
  return d;
}

void lstana::get_ring_radius(TH1D *h1_ring_radius, TH1D *h1_delta_ring_radius, TGraph *gr_data, TGraph *gr_weight, Double_t cx_reco, Double_t cy_reco, Double_t r_reco){
  Double_t x,y,w;
  Double_t r;
  for(Int_t i = 0;i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    gr_weight->GetPoint(i,w,w);
    r = TMath::Sqrt( (x - cx_reco) * (x - cx_reco) + (y - cy_reco) * (y - cy_reco) );
    h1_ring_radius->Fill(r,w);
    h1_delta_ring_radius->Fill(r - r_reco, w);
  }
}

void lstana::get_partial_ring_radius(TH1D *h1_ring_radius, TH1D *h1_delta_ring_radius, TGraph *gr_test, TGraph *gr_data, TGraph *gr_weight, Double_t cx_reco, Double_t cy_reco, Double_t r_reco){
  Double_t x,y,w;
  Double_t r;
  TVector2 vc(cx_reco, cy_reco);
  TVector2 v;
  TVector2 v_shift;
  TVector2 v_shift_zero;
  Double_t d_angle;
  Double_t d_angle_min = 100;
  Double_t v_shift_phi_min;
  for(Int_t i = 0;i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    v.SetX(x);
    v.SetY(y);
    v_shift.SetX(x - cx_reco);
    v_shift.SetY(y - cy_reco);
    if( (v.Mod() > 0.0) && (vc.Mod() > 0.0)){
      d_angle = TMath::ACos(v*vc/v.Mod()/vc.Mod());
      d_angle = TMath::Abs(d_angle)*180.0/TMath::Pi();
    }
    if(d_angle_min>d_angle){
      d_angle_min = d_angle;
      v_shift_phi_min = v_shift.Phi();
      v_shift_zero.SetX(x - cx_reco);
      v_shift_zero.SetX(y - cy_reco);
      cout<<"d_angle_min = "<<d_angle_min<<endl;
      //cout<<"dotprod = "<<dotprod<<endl;
      //gr_weight->GetPoint(i,w,w);
      //r = TMath::Sqrt( (x - cx_reco) * (x - cx_reco) + (y - cy_reco) * (y - cy_reco) );
      //h1_ring_radius->Fill(r,w);
      //h1_delta_ring_radius->Fill(r - r_reco, w);
    }
  }
  //
  //
  //
  //
  //
  //
  for(Int_t i = 0;i<gr_data->GetN();i++){
    gr_data->GetPoint(i,x,y);
    v_shift.SetX(x - cx_reco);
    v_shift.SetY(y - cy_reco);
    if( (v_shift_zero.Mod() > 0.0) && (v_shift.Mod() > 0.0)){
      d_angle = TMath::ACos(v_shift_zero*v_shift/v_shift_zero.Mod()/v_shift.Mod());
      d_angle = TMath::Abs(d_angle)*180.0/TMath::Pi();
    }
    cout<<"d_angle = "<<d_angle<<endl;
    if(d_angle<15.0){
      gr_weight->GetPoint(i,w,w);
      //if(w>10){
      r = TMath::Sqrt( (x - cx_reco) * (x - cx_reco) + (y - cy_reco) * (y - cy_reco) );
      h1_ring_radius->Fill(r,w);
      h1_delta_ring_radius->Fill(r - r_reco, w);
      gr_test->SetPoint(gr_test->GetN(),x,y);
      //}
    }
  }
  gr_test->SetPoint(gr_test->GetN(),2000,2000);
  gr_test->SetPoint(gr_test->GetN(),-2000,-2000);
}
