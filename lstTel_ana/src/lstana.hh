#ifndef lstana_hh
#define lstana_hh

//My
#include "lstanabase.hh"

//root
#include <TROOT.h>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;
class TRandom3;

class lstana: public lstanabase {
public:
  lstana(TString fileList) : lstanabase(fileList)
  {
  }

  lstana(TString file, Int_t key) : lstanabase(file, key)
  {
  }

  void Loop(TString histOut);
  void Loop01(TString histOut);

  void get_trk_theta(Double_t &trk_theta, Double_t &trk_phi);
  Double_t get_canonical_phi_from_Vphi(Double_t phiv);
  Double_t get_trk_mom_GeV();
  
  const Double_t _lst_Rmax_m = 12.12;                      // m
  const Double_t _lst_effective_focal_length_m = 29.30565; // m
  bool get_trk_mirror_impact_point(Double_t &trk_mirror_impact_point_X, Double_t &trk_mirror_impact_point_Y, Double_t trk_mirror_impact_point_Z);
  
  void simulate_optical_systems(TH2D *h2, TH1D *h1_phi);
  
  void save_to_csv( Double_t trk_mir_impact_X, Double_t trk_mir_impact_Y, Double_t trk_mir_impact_Z,
		    Double_t save_percentage, Int_t event_counter);

  Double_t get_ring_Completeness(TH1D *h1, Double_t threshold = 50.0);

  void TH2D_divide( TH2D *h2_w, TH2D *h2, TH2D *h2_norm);
  void TH1D_divide( TH1D *h1_w, TH1D *h1, TH1D *h1_norm);
  
  void get_weight_data_from_TH2D( TH2D *h2, TGraph *gr_data, TGraph *gr_weight, Double_t treshold);
  void get_optimum_circular_fit_Chaudhuri(TGraph *gr_data, TGraph *gr_weight, Double_t &cx_reco,  Double_t &cy_reco,  Double_t &r_reco);
  void gen_ring(TGraph *gr, Int_t np, Double_t x0, Double_t y0, Double_t R);
  void get_phi_dist_and_clean_ring(TH1D *h1_deg, TH2D *h2_PosY_vs_PosX_clean, Float_t cx_reco, Float_t cy_reco, Float_t r_reco, Float_t d_r_reco);

  void predict_phi_dist(TH1D *h1_phi_deg, TH1D *h1_PosZ, Double_t n_norm);
  void predict_phi_dist_camera(TH1D *h1_phi_deg, TH1D *h1_PosZ, Double_t n_norm);
  void reco_muon_impact_point( Double_t phot_z_hit, Double_t sensitive_dZ, Double_t &impact_x, Double_t &impact_y, Double_t &impact_z);
  Double_t d_zero( Double_t impact_x, Double_t impact_y, Double_t R_mirror, Double_t norm, Double_t phi);
  void get_ring_radius(TH1D *h1_ring_radius, TH1D *h1_delta_ring_radius, TGraph *gr_data, TGraph *gr_weight, Double_t cx_reco, Double_t cy_reco, Double_t r_reco);
  void get_partial_ring_radius(TH1D *h1_ring_radius, TH1D *h1_delta_ring_radius, TGraph *gr_test, TGraph *gr_data, TGraph *gr_weight, Double_t cx_reco, Double_t cy_reco, Double_t r_reco);
  
  TRandom3 *_rnd;
  
};

#endif
