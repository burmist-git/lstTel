//C, C++
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include <time.h>

#include "TMinuit.h"

using namespace std;

const Double_t _lst_Rmax_m = 12.12;                      // m
const Double_t _lst_effective_focal_length_m = 29.30565; // m

void get_gr_from_h2(TGraph *gr, TH2D *h2, Double_t phot_cut);

void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R);
void get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R);
void get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg);
void gen_Graph3D(TGraph *gr,TGraph2DErrors *gr2D);
Double_t plane_3D(Double_t *x, Double_t *par);
TGraph *_gr_to_fit = new TGraph();
TH1D *_h1_to_fit;
double equation_of_circle(double x, double y, double *par);
void fcn(int &npar, double *gin, double &f, double *par, int iflag);

Double_t get_angles_from_cam_pos( Double_t val_m, Double_t cam_m_min, Double_t cam_m_max);
Double_t get_canonical_phi_from_Vphi(Double_t phiv);

void fit_ring_with_Minuit(Double_t x0in, Double_t y0in, Double_t Rin,
			  Double_t &x0out, Double_t &y0out, Double_t &Rout,
			  Double_t &x0outerr, Double_t &y0outerr, Double_t &Routerr);

void fit_phi_with_Minuit( int if_local_ring_in, int if_mu_inclination_in, double ksi_deg_in,
			  double mu_opt_axis_angla_deg_in, double theta_c_deg_in, double R_mirror_m_in,
			  double Ampl_in, double rho_m_in, double phi0_deg_in,
			  double &Ampl_out, double &rho_m_out, double &phi0_deg_out,
			  double &Ampl_outerr, double &rho_m_outerr, double &phi0_deg_outerr);

//void test_D_rho_phi0(TH1D *h1_phi_test, int if_local_ring_in, double theta_c_deg_in, double R_mirror_m_in, double Ampl_in, double rho_m_in, double phi0_deg_in);
void test_D_rho_phi0(TH1D *h1_phi_test, int if_local_ring_in, int if_mu_inclination_in, double ksi_deg_in, double mu_opt_axis_angla_deg_in, double theta_c_deg_in,
		     double R_mirror_m_in, double Ampl_in, double rho_m_in, double phi0_deg_in);

double D_rho_phi0( double phi_deg, double *par);

Int_t fit_muon_ring(Int_t evID = 1999){
  //
  TRandom3 *rnd = new TRandom3(12312312);
  //
  TString fileN01;
  fileN01 = "./histSingle.root";
  //
  TFile *f01 = new TFile(fileN01.Data());
  //
  
  //
  TString h2_cam_ev_str = "camimage/h2_cam_";
  h2_cam_ev_str += evID;
  h2_cam_ev_str += "_ev";
  //
  TString h1_phi_ev_str = "camimage/h1_phi_";
  h1_phi_ev_str += evID;
  h1_phi_ev_str += "_ev";
  //
  TString h1_true_ev_str = "camimage/h1_true_";
  h1_true_ev_str += evID;
  h1_true_ev_str += "_ev";
  //
  TH2D *h2 = (TH2D*)f01->Get(h2_cam_ev_str.Data());
  _h1_to_fit = (TH1D*)f01->Get(h1_phi_ev_str.Data());
  TH1D *h1_true_info = (TH1D*)f01->Get(h1_true_ev_str.Data());  
  //
  //TH2D *h2 = (TH2D*)f01->Get("camimage/h2_cam_2_ev");
  //_h1_to_fit = (TH1D*)f01->Get("camimage/h1_phi_2_ev");
  //TH1D *h1_true_info = (TH1D*)f01->Get("camimage/h1_true_2_ev");  
  //
  Double_t trk_mirror_impact_point_X = h1_true_info->GetBinContent(1);
  Double_t trk_mirror_impact_point_Y = h1_true_info->GetBinContent(2);
  Double_t trk_mirror_impact_point_Z = h1_true_info->GetBinContent(3);
  Double_t trk_theta_deg = h1_true_info->GetBinContent(4);
  Double_t trk_phi_deg = h1_true_info->GetBinContent(5);
  //
  cout<<"trk_mirror_impact_point_X "<<trk_mirror_impact_point_X<<endl
      <<"trk_mirror_impact_point_Y "<<trk_mirror_impact_point_Y<<endl
      <<"trk_mirror_impact_point_Z "<<trk_mirror_impact_point_Z<<endl
      <<"trk_theta_deg             "<<trk_theta_deg<<endl
      <<"trk_phi_deg               "<<trk_phi_deg<<endl;

  //
  TGraph *gr = new TGraph();
  get_gr_from_h2( gr, h2, 40);
  //
  TH1D *h1_phi_test = new TH1D("h1_phi_test","h1_phi_test", 40,-200.0,200.0);
  
  //
  Double_t x0app, y0app, Rapp;
  Double_t x0app_average = 0.0;
  Double_t y0app_average = 0.0;
  Double_t Rapp_average = 0.0;
  vector<Double_t> x0app_v;
  vector<Double_t> y0app_v;
  vector<Double_t> Rapp_v;
  for(Int_t i = 0;i<50;i++){
    get_approximate_ring_parameters( gr, rnd, x0app, y0app, Rapp);
    x0app_v.push_back(x0app);
    y0app_v.push_back(y0app);
    Rapp_v.push_back(Rapp);
    x0app_average += x0app;
    y0app_average += y0app;
    Rapp_average += Rapp;
    //cout<<"x0app "<<x0app<<endl
    //    <<"y0app "<<y0app<<endl
    //    <<"Rapp  "<<Rapp<<endl;
  }
  //
  x0app_average /= x0app_v.size();
  y0app_average /= x0app_v.size();
  Rapp_average /= x0app_v.size();
  //
  Double_t x0out, y0out, Rout;
  Double_t x0outerr, y0outerr, Routerr;
  fit_ring_with_Minuit( x0app_average, y0app_average, Rapp_average,
			x0out, y0out, Rout,
			x0outerr, y0outerr, Routerr);
  //
  TGraph *gr_app_average_r0 = new TGraph();
  gr_app_average_r0->SetNameTitle("gr_app_average_r0","gr_app_average_r0");  
  gr_app_average_r0->SetPoint( 0, x0app_average, y0app_average);
  gr_app_average_r0->SetMarkerStyle(43);
  gr_app_average_r0->SetMarkerColor(kMagenta+3);
  gr_app_average_r0->SetMarkerSize(3.0);
  //
  TGraph *gr_app_average_ring = new TGraph();
  gr_app_average_ring->SetNameTitle("gr_app_average_ring","gr_app_average_ring");  
  gr_app_average_ring->SetMarkerStyle(7);
  gr_app_average_ring->SetMarkerColor(kMagenta+3);
  gr_app_average_ring->SetMarkerSize(3.0);
  gen_ring(gr_app_average_ring, 360, x0app_average, y0app_average, Rapp_average);
  //
  TGraph *gr_fit_ring = new TGraph();
  gr_fit_ring->SetNameTitle("gr_fit_ring","gr_fit_ring");
  gr_fit_ring->SetMarkerStyle(7);
  gr_fit_ring->SetMarkerColor(kRed);
  gr_fit_ring->SetMarkerSize(2.0);
  gen_ring(gr_fit_ring, 360, x0out, y0out, Rout);
  //
  TGraph *gr_app_r0 = new TGraph();
  gr_app_r0->SetNameTitle("gr_app_r0","gr_app_r0");
  for(unsigned int ii = 0;ii<x0app_v.size();ii++)
    gr_app_r0->SetPoint( ii, x0app_v.at(ii), y0app_v.at(ii));
  //
  gr->SetMarkerStyle(20);
  gr_app_r0->SetMarkerStyle(43);
  gr_app_r0->SetMarkerColor(kBlue);
  gr_app_r0->SetMarkerSize(2.0);
  //  
  TGraph *gr_R_app_average = new TGraph();
  gr_R_app_average->SetNameTitle("gr_R_app_average","gr_R_app_average");
  TH1D *h1_theta_app_average_deg = new TH1D("h1_theta_app_average_deg","h1_theta_app_average_deg", 36, 0.0, 360.0); 
  get_Rvs_theta_and_theta_dist( gr, x0app_average, y0app_average, Rapp_average, gr_R_app_average, h1_theta_app_average_deg);
  //
  TCanvas *c1 = new TCanvas("c1","c1",10,10,1800,600);
  c1->Divide(3,1);
  c1->cd(1);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  //gStyle->SetOptStat(kFALSE);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(gr);
  //mg->Add(gr_frame);
  mg->Add(gr_app_r0);  
  mg->Add(gr_app_average_r0);
  mg->Add(gr_app_average_ring);
  mg->Add(gr_fit_ring);
  mg->Draw("AP");
  //
  c1->cd(2);
  gr_R_app_average->SetMarkerStyle(20);
  gr_R_app_average->Draw("AP");  
  //
  c1->cd(3);
  h1_theta_app_average_deg->SetLineColor(kBlack);
  h1_theta_app_average_deg->SetLineWidth(2.0);
  h1_theta_app_average_deg->Draw();
  //
  //
  //
  //
  //
  TGraph2DErrors *gr2D = new TGraph2DErrors();
  TGraph2DErrors *gr2D_average = new TGraph2DErrors();
  gen_Graph3D( gr,gr2D);
  gen_Graph3D( gr_app_average_ring,gr2D_average);
  //
  //
  Double_t boxL = 5;
  TGraph2D *gr2D_box = new TGraph2D();
  gr2D_box->SetTitle("gr2D_box");
  gr2D_box->SetName("gr2D_box");
  //
  gr2D_box->SetPoint(gr2D_box->GetN(),-boxL/2.0,-boxL/2.0,-boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),-boxL/2.0, boxL/2.0,-boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),-boxL/2.0,-boxL/2.0, boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),-boxL/2.0, boxL/2.0, boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),boxL/2.0,-boxL/2.0,-boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),boxL/2.0, boxL/2.0,-boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),boxL/2.0,-boxL/2.0, boxL/2.0);
  gr2D_box->SetPoint(gr2D_box->GetN(),boxL/2.0, boxL/2.0, boxL/2.0);
  //
  //

  TCanvas *c2 = new TCanvas("c2","c2",10,10,1400,600);
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gROOT->ForceStyle();
  gStyle->SetStatColor(kWhite);
  gStyle->SetOptStat(kFALSE);
  
  //c2->Divide(1,2);

  int if_local_ring_in = 1;
  int if_mu_inclination_in = 0;
  double ksi_deg_in = 0.0;
  double mu_opt_axis_angla_deg_in = 0.0;
  //double theta_c_deg_in = TMath::ATan(Rout/_lst_effective_focal_length_m)*180.0/TMath::Pi();
  double theta_c_deg_in = Rout;
  double R_mirror_m_in = _lst_Rmax_m;
  double Ampl_in = _h1_to_fit->Integral()/36.0;
  //double rho_m_in = TMath::Sqrt(trk_mirror_impact_point_X*trk_mirror_impact_point_X + trk_mirror_impact_point_Y*trk_mirror_impact_point_Y)*0.0;
  //double rho_m_in = TMath::Sqrt(trk_mirror_impact_point_X*trk_mirror_impact_point_X + trk_mirror_impact_point_Y*trk_mirror_impact_point_Y);
  double rho_m_true = TMath::Sqrt(trk_mirror_impact_point_X*trk_mirror_impact_point_X + trk_mirror_impact_point_Y*trk_mirror_impact_point_Y);
  double rho_m_in = 10.0;
  //double rho_m_in = rho_m_true;
  double phi0_deg_in = _h1_to_fit->GetBinCenter(_h1_to_fit->GetMaximumBin());
  
  double Ampl_out, rho_m_out, phi0_deg_out;
  double Ampl_outerr, rho_m_outerr, phi0_deg_outerr;
  
  //test_D_rho_phi0(h1_phi_test, if_local_ring_in, ksi_deg, theta_c_deg_in, R_mirror_m_in, Ampl_in, rho_m_in, phi0_deg_in);

  //c2->cd(1);
  _h1_to_fit->SetLineWidth(2);
  h1_phi_test->SetLineWidth(2);
  _h1_to_fit->SetLineColor(kBlack);
  h1_phi_test->SetLineColor(kRed);
  _h1_to_fit->Draw();
  h1_phi_test->Draw("same");
  
  //fit_phi_with_Minuit( if_local_ring_in, theta_c_deg_in, R_mirror_m_in,
  // 		       Ampl_in, rho_m_in, phi0_deg_in,
  //		       Ampl_out, rho_m_out, phi0_deg_out,
  //		       Ampl_outerr, rho_m_outerr, phi0_deg_outerr);  

  fit_phi_with_Minuit( if_local_ring_in, if_mu_inclination_in, ksi_deg_in,
  		       mu_opt_axis_angla_deg_in, theta_c_deg_in, R_mirror_m_in, Ampl_in, rho_m_in, phi0_deg_in,
  		       Ampl_out, rho_m_out, phi0_deg_out,
  		       Ampl_outerr, rho_m_outerr, phi0_deg_outerr);  

  //
  //
  TVector2 impact_par_v(trk_mirror_impact_point_X,trk_mirror_impact_point_Y);
  TVector2 impact_par_rot_v = impact_par_v.Rotate(-TMath::Pi());
  //cout<<"rho_m_true   "<<rho_m_true<<endl
  //   <<"Phi0_true    "<<impact_par_rot_v.Phi()*180.0/TMath::Pi()<<endl;
  cout<<"rho_m_true   "<<rho_m_true<<endl
      <<"Phi0_true    "<<get_canonical_phi_from_Vphi(impact_par_rot_v.Phi())*180.0/TMath::Pi()<<endl;

  //
  //
  cout<<"info: ID  "<<evID<<endl
      <<"info: rho  "<<rho_m_true<<endl
      <<"info: delta rho  "<<rho_m_true - rho_m_out<<endl
      <<"info: delta Phi0 "<<get_canonical_phi_from_Vphi(impact_par_rot_v.Phi())*180.0/TMath::Pi() - phi0_deg_out<<endl;
  //
  //
  
  test_D_rho_phi0(h1_phi_test, if_local_ring_in, if_mu_inclination_in, ksi_deg_in,
  		  mu_opt_axis_angla_deg_in, theta_c_deg_in, R_mirror_m_in, Ampl_out, rho_m_out, phi0_deg_out);
  //test_D_rho_phi0(h1_phi_test, if_local_ring_in, if_mu_inclination_in, ksi_deg_in,
  //		  mu_opt_axis_angla_deg_in, theta_c_deg_in, R_mirror_m_in, Ampl_in, rho_m_in, phi0_deg_in);


  //TVector2 test_vec(1,-1);
  //cout<<"test_vec.Phi() "<<test_vec.Phi()*180.0/TMath::Pi()<<endl;
  
  /*
  //
  c2->cd(2);
  gr2D_average->SetTitle("");
  gr2D_box->Draw("P");
  gr2D_average->Draw("sameP");
  //
  c2->cd(3);
  gr2D_average->SetTitle("");
  gr2D_box->Draw("P");
  gr2D_average->Draw("sameP");
  gr2D->Draw("sameP");
  */
  //
  //
  //Double_t x_fit_min = -1.2;
  //Double_t x_fit_max = -1.2;
  //Double_t y_fit_min = -1.2;
  //Double_t y_fit_max = -1.2;
  //Int_t ndim = 2;
  //
  //TF2 TF2(const char* name, Double_t(*)(Double_t*,Double_t*) fcn, Double_t xmin = 0, Double_t xmax = 1, Double_t ymin = 0, Double_t ymax = 1, Int_t npar = 0, Int_t ndim = 2)
  //TF2 *plane_fit_f = new TF2("plane_fit_f", plane_3D, x_fit_min, x_fit_max, y_fit_min, y_fit_max, npar, ndim);
  //plane_fit_f->SetParameter(0,x0app_average);
  //plane_fit_f->SetParameter(1,y0app_average);
  //plane_fit_f->SetParameter(2,(x0app_average*x0app_average + y0app_average*y0app_average - Rapp_average*Rapp_average));
  //plane_fit_f->FixParameter(0,x0app_average);
  //plane_fit_f->FixParameter(1,y0app_average);
  //
  //gr2D_average->Fit(plane_fit_f,"");
  //gr2D->Fit(plane_fit_f);
  //
  //
  return 0;
}

Double_t get_angles_from_cam_pos(Double_t val_m, Double_t cam_m_min, Double_t cam_m_max){
  Double_t cam_angle_max = TMath::ATan(cam_m_max/_lst_effective_focal_length_m)*180.0/TMath::Pi();
  Double_t cam_angle_min = -cam_angle_max;
  Double_t k = (cam_angle_max - cam_angle_min)/(cam_m_max - cam_m_min);
  return k*val_m;
}

void get_gr_from_h2(TGraph *gr, TH2D *h2, Double_t phot_cut){
  Double_t val;
  Double_t x;
  Double_t y;
  Double_t x_deg;
  Double_t y_deg;
  Double_t cam_m_min = h2->GetXaxis()->GetBinLowEdge(1);
  Double_t cam_m_max = h2->GetXaxis()->GetBinUpEdge(h2->GetNbinsX());    
  for(Int_t i = 1;i<=h2->GetNbinsX();i++){
    for(Int_t j = 1;j<=h2->GetNbinsY();j++){
      x = h2->GetXaxis()->GetBinCenter(i);
      y = h2->GetYaxis()->GetBinCenter(j);
      val = h2->GetBinContent(i,j);
      if(val>phot_cut){
	x_deg = get_angles_from_cam_pos(x,cam_m_min,cam_m_max);
	y_deg =	get_angles_from_cam_pos(y,cam_m_min,cam_m_max);
	gr->SetPoint(gr->GetN(),x_deg,y_deg);
	_gr_to_fit->SetPoint(_gr_to_fit->GetN(),x_deg,y_deg);
      }
    }
  }
}
//
//GetBinWidth
//

void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R){
  Double_t phi = 0.0;
  TVector2 rc(x0,y0);
  for(Int_t i = 0;i<np;i++){
    TVector2 p;
    p.SetMagPhi(R,2*TMath::Pi()/(np-1)*i);
    TVector2 pt = rc + p;
    gr->SetPoint( i, pt.X(), pt.Y());
  }
}

void get_approximate_ring_parameters(TGraph *gr, TRandom3 *rnd, Double_t &x0, Double_t &y0, Double_t &R){
  Double_t p1x0, p1y0;
  Double_t p2x0, p2y0; 
  Double_t maxDist = 0.0;
  Int_t p1_id = (Int_t)rnd->Uniform(0,gr->GetN());
  Int_t p2_id = (Int_t)rnd->Uniform(0,gr->GetN());
  gr->GetPoint(p1_id, p1x0, p1y0);
  TVector2 p1(p1x0, p1y0);
  gr->GetPoint(p2_id, p2x0, p2y0);
  TVector2 p2(p2x0, p2y0);
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i, p2x0, p2y0);
    p2.Set(p2x0, p2y0);
    TVector2 dp = p1 - p2;
    if(maxDist<dp.Mod()){
      maxDist = dp.Mod();
      p2_id = i;
    }
  }
  //
  gr->GetPoint(p2_id, p2x0, p2y0);
  p2.Set(p2x0, p2y0);
  //  
  TVector2 pm = (p1 + p2)/2.0;
  x0 = pm.X();
  y0 = pm.Y();
  R = maxDist/2.0;
}

void get_Rvs_theta_and_theta_dist(TGraph *gr, Double_t x0, Double_t y0, Double_t R, TGraph *gr_R, TH1D *h1_theta_deg){
  Double_t x, y;
  TVector2 pc(x0, y0);
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i, x, y);
    TVector2 p(x,y);
    TVector2 dp = p - pc;
    gr_R->SetPoint(i,dp.Phi()*180.0/TMath::Pi(),(dp.Mod() - R)/R);
    h1_theta_deg->Fill(dp.Phi()*180.0/TMath::Pi());
  }
}

void gen_Graph3D(TGraph *gr,TGraph2DErrors *gr2D){
  Double_t x, y;
  Double_t xt, yt,zt;
  for(Int_t i = 0;i<gr->GetN();i++){
    gr->GetPoint(i,x,y);
    xt = -2.0*x;
    yt = -2.0*y;
    zt = x*x + y*y;
    gr2D->SetPoint(i,xt,yt,zt);
    gr2D->SetPointError(i,100,100,100);
  }
}

Double_t plane_3D(Double_t *x, Double_t *par){
  //
  //Double_t X0 = par[0];
  //Double_t Y0 = par[1];
  //Double_t Z0 = 1.0;
  //Double_t D  = par[2];
  //
  //Double_t x = x[0];
  //Double_t y = x[1];
  //
  // X0x + Y0y + 1.0*z + D = 0
  return -par[0]*x[0] - par[1]*x[2] - par[2];
}

void fcnD(int &npar, double *gin, double &f, double *par, int iflag){
   double chisq = 0;
   double phi_deg;
   double npotonst;
   double npotonst_err;
   double delta;
   for (int i = 1; i<=_h1_to_fit->GetNbinsX(); i++){
     npotonst = _h1_to_fit->GetBinContent(i);
     npotonst_err = _h1_to_fit->GetBinError(i);
     phi_deg = _h1_to_fit->GetBinCenter(i);
     //delta = (npotonst - D_rho_phi0( phi_deg, par))/npotonst_err;
     if(npotonst_err>0)
       delta = (npotonst - D_rho_phi0( phi_deg, par))/npotonst_err;
     else
       delta = (npotonst - D_rho_phi0( phi_deg, par));
     chisq += delta*delta;
   }
   f = chisq;
}

void fcn(int &npar, double *gin, double &f, double *par, int iflag){
   double chisq = 0;
   double x, y;
   double delta;
   for (int i = 0; i<_gr_to_fit->GetN(); i++){
     _gr_to_fit->GetPoint( i, x, y);
     delta = equation_of_circle(x, y, par);
     if(delta>0.0)
       delta = delta*0.5;
     chisq += TMath::Abs(delta);
   }
   f = chisq;
}

double equation_of_circle(double x, double y, double *par){
  return par[2]*par[2] - (x-par[0])*(x-par[0]) - (y-par[1])*(y-par[1]);
}

double D_rho_phi0( double phi_deg, double *par){
  if(phi_deg<-180.0)
    return 0.0;
  if(phi_deg>180.0)
    return 0.0;
  int    if_local_ring = (int)par[0];
  int    if_mu_inclination = (int)par[1];
  double ksi_deg = par[2];
  double ksi = ksi_deg/180.0*TMath::Pi();
  double mu_opt_axis_angle_deg = par[3];
  double mu_opt_axis_angle = mu_opt_axis_angle_deg/180.0*TMath::Pi();
  double theta_c_deg = par[4];
  double theta_c = theta_c_deg/180.0*TMath::Pi();
  double R_mirror_m = par[5];
  double Ampl = par[6];
  double rho_m = par[7];
  double rho_R_m = rho_m/R_mirror_m;
  double phi0_deg = par[8];
  double phi0 = phi0_deg/180.0*TMath::Pi();
  double phi = phi_deg/180.0*TMath::Pi();
  //
  double dNdPhi = 0.0;
  double phi_m_phi0_max = 0.0;
  double D_zero = 0.0;
  double mu_inclination_corr = 0.0;
  double mu_inclination_corr_denominator = 1.0;
  // local_ring
  if(if_local_ring == 1){
    if(if_mu_inclination == 0){
      if(rho_R_m<=1.0){
	dNdPhi = 2.0*Ampl*(TMath::Sin(2.0*theta_c)*R_mirror_m*TMath::Sqrt(1.0 - rho_R_m*rho_R_m/4.0*TMath::Sin(phi-phi0)*TMath::Sin(phi-phi0)) +
			   rho_R_m/2.0*TMath::Cos(phi-phi0));
      }
      else{
	dNdPhi = 0.0;
      }
    }
    else{
      if(rho_R_m<=1.0){
	D_zero = R_mirror_m*TMath::Sqrt(1.0 - rho_R_m*rho_R_m*TMath::Sin(phi-phi0)*TMath::Sin(phi-phi0)) + rho_R_m*TMath::Cos(phi-phi0);
	mu_inclination_corr = (1.0 + TMath::Tan(theta_c)*TMath::Tan(mu_opt_axis_angle)*TMath::Cos(ksi-phi));
	mu_inclination_corr_denominator = TMath::Sqrt(1.0 + TMath::Tan(mu_opt_axis_angle)*TMath::Tan(mu_opt_axis_angle)*TMath::Cos(ksi-phi)*TMath::Cos(ksi-phi));
	dNdPhi = Ampl*TMath::Sin(theta_c)*TMath::Sin(theta_c)*D_zero/TMath::Tan(theta_c)*mu_inclination_corr/mu_inclination_corr_denominator;
      }
      else{
	dNdPhi = 0.0;
      }
    }
  }
  else{
    if(rho_R_m>1.0){
      phi_m_phi0_max = TMath::ASin(1.0/rho_R_m);
      if(TMath::Abs(phi-phi0)<=phi_m_phi0_max){
	dNdPhi = Ampl*TMath::Sin(2.0*theta_c)*2.0*R_mirror_m*TMath::Sqrt(1.0 - rho_R_m*rho_R_m*TMath::Sin(phi-phi0)*TMath::Sin(phi-phi0));
      }
      else{
	dNdPhi = 0.0;
      }
    }
    else{
      dNdPhi = 0.0;
    }  
  }
  return dNdPhi;
}

void fit_ring_with_Minuit(Double_t x0in, Double_t y0in, Double_t Rin,
			  Double_t &x0out, Double_t &y0out, Double_t &Rout,
			  Double_t &x0outerr, Double_t &y0outerr, Double_t &Routerr){
  //
  Int_t npar = 3;
  TMinuit *gMinuit = new TMinuit(npar);
  //gMinuit->SetPrintLevel(-1.0);
  gMinuit->SetFCN(fcn); 
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // 
  // Set starting values and step sizes for parameters
  gMinuit->mnparm(0, "x0", x0in, 0.01, 0,0,ierflg);
  gMinuit->mnparm(1, "y0", y0in, 0.01, 0,0,ierflg);
  gMinuit->mnparm(2, "R", Rin, 0.01, 0,0,ierflg);
  //

  // Now ready for minimization step
  arglist[0] = 1000000;
  arglist[1] = 1.0;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  // Print results
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  //
  gMinuit->GetParameter(0, x0out, x0outerr);
  gMinuit->GetParameter(1, y0out,  y0outerr);
  gMinuit->GetParameter(2, Rout, Routerr);
  //
  cout<<"x0out "<<x0out<<endl
      <<"y0out "<<y0out<<endl
      <<"Rout  "<<Rout<<endl;
  //
}

void test_D_rho_phi0(TH1D *h1_phi_test, int if_local_ring_in, int if_mu_inclination_in, double ksi_deg_in, double mu_opt_axis_angla_deg_in, double theta_c_deg_in, double R_mirror_m_in, double Ampl_in, double rho_m_in, double phi0_deg_in){
  //
  double phi_val;
  double par[9];
  //
  par[0] = if_local_ring_in;
  par[1] = if_mu_inclination_in;
  par[2] = ksi_deg_in;
  par[3] = mu_opt_axis_angla_deg_in;
  par[4] = theta_c_deg_in;
  par[5] = R_mirror_m_in;
  par[6] = Ampl_in;
  par[7] = rho_m_in;
  par[8] = phi0_deg_in;
  //
  for(Int_t i = 1;i<=h1_phi_test->GetNbinsX();i++){
    phi_val = h1_phi_test->GetBinCenter(i);
    h1_phi_test->SetBinContent(i, D_rho_phi0( phi_val, par));
  }
}

void fit_phi_with_Minuit( int if_local_ring_in, int if_mu_inclination_in, double ksi_deg_in,
			  double mu_opt_axis_angla_deg_in, double theta_c_deg_in, double R_mirror_m_in,
			  double Ampl_in, double rho_m_in, double phi0_deg_in,
			  double &Ampl_out, double &rho_m_out, double &phi0_deg_out,
			  double &Ampl_outerr, double &rho_m_outerr, double &phi0_deg_outerr){
  //
  //
  Int_t npar = 9;
  TMinuit *gMinuit = new TMinuit(npar);
  gMinuit->SetPrintLevel(1.0);
  gMinuit->SetFCN(fcnD); 
  double arglist[10];
  int ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  // 
  // Set starting values and step sizes for parameters

  gMinuit->mnparm(0, "if_local", if_local_ring_in, 0.0, 0,0,ierflg);
  gMinuit->mnparm(1, "if_mu_inclination", if_mu_inclination_in, 0.0, 0,0,ierflg);
  gMinuit->mnparm(2, "ksi_deg", ksi_deg_in, 0.0, 0,0,ierflg);
  gMinuit->mnparm(3, "mu_opt_axis_angle_deg", mu_opt_axis_angla_deg_in, 0.0, 0,0,ierflg);
  gMinuit->mnparm(4, "theta_c_deg", theta_c_deg_in, 0.0, 0,0,ierflg);
  gMinuit->mnparm(5, "R_mirror", R_mirror_m_in, 0.0, 0,0,ierflg);
  //
  gMinuit->mnparm(6, "Ampl", Ampl_in, 0.01, 0,0,ierflg);
  gMinuit->mnparm(7, "rho_m", rho_m_in, 0.01, 0,0,ierflg);
  gMinuit->mnparm(8, "phi0_deg", phi0_deg_in, 0.01, 0,0,ierflg);
  //

  // Now ready for minimization step
  arglist[0] = 50000;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
  
  // Print results
  double amin,edm,errdef;
  int nvpar,nparx,icstat;
  gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //gMinuit->mnprin(3,amin);
  //
  gMinuit->GetParameter(6, Ampl_out, Ampl_outerr);
  gMinuit->GetParameter(7, rho_m_out, rho_m_outerr);
  gMinuit->GetParameter(8, phi0_deg_out, phi0_deg_outerr);
  //
  cout<<"Ampl_out     "<<Ampl_out<<endl
      <<"rho_m_out    "<<rho_m_out<<endl
      <<"phi0_deg_out "<<phi0_deg_out<<endl;
  //
}

Double_t get_canonical_phi_from_Vphi(Double_t phiv){
  if( phiv >= 0.0 && phiv <= TMath::Pi())
    return phiv;
  return TMath::Pi()-phiv;
}
