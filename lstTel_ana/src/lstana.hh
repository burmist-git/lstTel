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

  void get_trk_theta(Double_t &trk_theta, Double_t &trk_phi);
  Double_t get_canonical_phi_from_Vphi(Double_t phiv);
  Double_t get_trk_mom_GeV();
  
  const Double_t _lst_Rmax_m = 12.12;                      // m
  const Double_t _lst_effective_focal_length_m = 29.30565; // m
  bool get_trk_mirror_impact_point(Double_t &trk_mirror_impact_point_X, Double_t &trk_mirror_impact_point_Y, Double_t trk_mirror_impact_point_Z);
  
  void simulate_optical_systems(TH2D *h2, TH1D *h1_phi);
  
  void save_to_csv( Double_t trk_mir_impact_X, Double_t trk_mir_impact_Y, Double_t trk_mir_impact_Z,
		    Double_t save_percentage, Int_t event_counter);
  
  TRandom3 *_rnd;
  
};

#endif
