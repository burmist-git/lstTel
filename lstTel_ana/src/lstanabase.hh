#ifndef lstanabase_hh
#define lstanabase_hh

#include <TROOT.h>

class TChain;
class TFile;
class TTree;
class TString;
class TBranch;
class TGraph;
class TH1D;
class TH2D;
class TProfile;

class lstanabase {

public :
  lstanabase(TString fileList);
  lstanabase(TString inFileName, Int_t keyID);
  ~lstanabase();
  Int_t GetEntry(Long64_t entry);
  Long64_t LoadTree(Long64_t entry);
  void Init(TTree *tree);
  void Loop(TString histOut);
  Bool_t Notify();
  void Show(Long64_t entry = -1);
  Int_t Cut(Long64_t entry);

protected :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  //Int_t           evt;
  //Int_t           run;
  //Float_t         pValue;
  //...
  //...
   Int_t           EventID;
   Int_t           BunchXID;
   Int_t           NTotPhot;
   Int_t           Nhits;
   Int_t           primType;
   Double_t        primMomX;
   Double_t        primMomY;
   Double_t        primMomZ;
   Double_t        primPosX;
   Double_t        primPosY;
   Double_t        primPosZ;
   Double_t        primTime;
   Double_t        trigTopL;
   Double_t        trigTopEdep;
   Double_t        trigBotL;
   Double_t        trigBotEdep;
   Int_t           nPhot;
   static const Int_t nMaxPhotons = 2000000;
   Int_t           TrackID[nMaxPhotons];   //[nPhot]
   Int_t           ParentID[nMaxPhotons];   //[nPhot]
   Double_t        Energy[nMaxPhotons];   //[nPhot]
   Double_t        Wavelength[nMaxPhotons];   //[nPhot]
   Double_t        Time[nMaxPhotons];   //[nPhot]
   Double_t        photPathLen[nMaxPhotons];   //[nPhot]
   Int_t           SecID[nMaxPhotons];   //[nPhot]
   Int_t           chID[nMaxPhotons];   //[nPhot]
   Double_t        PosX[nMaxPhotons];   //[nPhot]
   Double_t        PosY[nMaxPhotons];   //[nPhot]
   Double_t        PosZ[nMaxPhotons];   //[nPhot]
   Double_t        MomX[nMaxPhotons];   //[nPhot]
   Double_t        MomY[nMaxPhotons];   //[nPhot]
   Double_t        MomZ[nMaxPhotons];   //[nPhot]
   Double_t        trkMomX[nMaxPhotons];   //[nPhot]
   Double_t        trkMomY[nMaxPhotons];   //[nPhot]
   Double_t        trkMomZ[nMaxPhotons];   //[nPhot]
   Double_t        trkPosX[nMaxPhotons];   //[nPhot]
   Double_t        trkPosY[nMaxPhotons];   //[nPhot]
   Double_t        trkPosZ[nMaxPhotons];   //[nPhot]
   Double_t        trkT[nMaxPhotons];   //[nPhot]
   Double_t        trkLength[nMaxPhotons];   //[nPhot]
  //
  //---------------------------------------------------
  // ADD HERE :
  //Tree name
  //const TString treeName = "arich";
  const TString treeName = "T";
  static const Int_t nChannels = 10;



  //---------------------------------------------------
  
  // List of branches
  //TBranch        *b_evt;
  //TBranch        *b_run;
  //TBranch        *b_pValue;
  //...
  //...
  //
  //---------------------------------------------------
  // ADD HERE :
   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_BunchXID;   //!
   TBranch        *b_NTotPhot;   //!
   TBranch        *b_Nhits;   //!
   TBranch        *b_primType;   //!
   TBranch        *b_primMomX;   //!
   TBranch        *b_primMomY;   //!
   TBranch        *b_primMomZ;   //!
   TBranch        *b_primPosX;   //!
   TBranch        *b_primPosY;   //!
   TBranch        *b_primPosZ;   //!
   TBranch        *b_primTime;   //!
   TBranch        *b_trigTopL;   //!
   TBranch        *b_trigTopEdep;   //!
   TBranch        *b_trigBotL;   //!
   TBranch        *b_trigBotEdep;   //!
   TBranch        *b_nPhot;   //!
   TBranch        *b_TrackID;   //!
   TBranch        *b_ParentID;   //!
   TBranch        *b_Energy;   //!
   TBranch        *b_Wavelength;   //!
   TBranch        *b_Time;   //!
   TBranch        *b_photPathLen;   //!
   TBranch        *b_SecID;   //!
   TBranch        *b_chID;   //!
   TBranch        *b_PosX;   //!
   TBranch        *b_PosY;   //!
   TBranch        *b_PosZ;   //!
   TBranch        *b_MomX;   //!
   TBranch        *b_MomY;   //!
   TBranch        *b_MomZ;   //!
   TBranch        *b_trkMomX;   //!
   TBranch        *b_trkMomY;   //!
   TBranch        *b_trkMomZ;   //!
   TBranch        *b_trkPosX;   //!
   TBranch        *b_trkPosY;   //!
   TBranch        *b_trkPosZ;   //!
   TBranch        *b_trkT;   //!
   TBranch        *b_trkLength;   //!
  //---------------------------------------------------
  void tGraphInit(TGraph *gr[nChannels], TString grName, TString grTitle);
  void h1D1Init(TH1D *h1D1[nChannels],TString h1name, TString h1Title,
		Int_t Nbin, Float_t Vmin, Float_t Vmax);
  void h2D2Init(TH2D *h2D1[nChannels],TString h2name, TString h2Title,
                Int_t Nbin1, Float_t Vmin1, Float_t Vmax1,
                Int_t Nbin2, Float_t Vmin2, Float_t Vmax2);
  void tProfInit(TProfile *tprof[nChannels],TString prname, TString prTitle,
                 Int_t Nbin, Float_t Vmin, Float_t Vmax);
  double getUnixTimeFromTime(double d_year, double d_month, double d_day, double d_hour, double d_min, double d_sec);  
  //
  
};

#endif
