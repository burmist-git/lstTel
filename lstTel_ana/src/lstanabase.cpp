//my
#include "lstanabase.hh"

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

//C, C++
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <iomanip>
#include <time.h>
#include <bits/stdc++.h>

using namespace std;

lstanabase::lstanabase(TString fileList) : fChain(0) 
{
  ifstream indata;
  TString rootFileName;
  //TChain *theChain = new TChain("arich");
  TChain *theChain = new TChain(treeName.Data());
  indata.open(fileList.Data()); 
  assert(indata.is_open());  
  while (indata  >> rootFileName ){
    if(indata.eof()){
      std::cout<<"EOF"<<std::endl;
      break;
    }
    cout<<"        adding "<<rootFileName<<endl;
    theChain->Add(rootFileName.Data(),-1);
  }
  indata.close();
  Init(theChain);
}

lstanabase::lstanabase(TString inFileName, Int_t keyID) : fChain(0) 
{
  if(keyID == 1){
    ifstream indata;
    //TChain *theChain = new TChain("T");
    TChain *theChain = new TChain(treeName.Data());
    cout<<"        adding "<<inFileName<<endl;
    theChain->Add(inFileName.Data(),-1);
    Init(theChain);
  }
  else
    assert(0);
}

void lstanabase::tGraphInit(TGraph *gr[nChannels], TString grName, TString grTitle){
  Int_t i;
  TString grNameh;
  TString grTitleh;
  for(i = 0;i<nChannels;i++){
    grNameh = grName; grNameh += "_"; grNameh += "ch_"; grNameh += i;
    grTitleh = grTitle; grTitleh += " "; grTitleh += "ch "; grTitleh += i;
    gr[i] = new TGraph();
    gr[i]->SetTitle(grTitleh.Data());
    gr[i]->SetName(grNameh.Data());
  }
}

void lstanabase::h1D1Init(TH1D *h1D1[nChannels],TString h1name, TString h1Title,
		      Int_t Nbin, Float_t Vmin, Float_t Vmax){
  Int_t i;
  TString h1nameh;
  TString h1Titleh;
  for(i = 0;i<nChannels;i++){    
    h1nameh = h1name; h1nameh += "_"; h1nameh += "ch_"; h1nameh += i;
    h1Titleh = h1Title; h1nameh += " "; h1Titleh += "ch "; h1Titleh += i;
    h1D1[i] = new TH1D(h1nameh.Data(), h1Titleh.Data(),
                       Nbin, Vmin, Vmax);
  }
}

void lstanabase::h2D2Init(TH2D *h2D1[nChannels],TString h2name, TString h2Title,
		      Int_t Nbin1, Float_t Vmin1, Float_t Vmax1,
		      Int_t Nbin2, Float_t Vmin2, Float_t Vmax2){
  Int_t i;
  TString h2nameh;
  TString h2Titleh;
  for(i = 0;i<nChannels;i++){    
    h2nameh = h2name; h2nameh += "_"; h2nameh += "ch_"; h2nameh += i;
    h2Titleh = h2Title; h2nameh += " "; h2Titleh += "ch "; h2Titleh += i;
    h2D1[i] = new TH2D(h2nameh.Data(), h2Titleh.Data(),
                       Nbin1, Vmin1, Vmax1,
		       Nbin2, Vmin2, Vmax2);
  }  
}

void lstanabase::tProfInit(TProfile *tprof[nChannels],TString prname, TString prTitle,
		       Int_t Nbin, Float_t Vmin, Float_t Vmax){
  Int_t i;
  TString prnameh;
  TString prTitleh;
  for(i = 0;i<nChannels;i++){    
    prnameh = prname; prnameh += "_"; prnameh += "ch_"; prnameh += i;
    prTitleh = prTitle; prnameh += " "; prTitleh += "ch "; prTitleh += i;
    tprof[i] = new TProfile(prnameh.Data(), prTitleh.Data(), Nbin, Vmin, Vmax,"");
  }
}

double lstanabase::getUnixTimeFromTime(double d_year, double d_month, double d_day, double d_hour, double d_min, double d_sec){
  //time_t timer;
  struct tm y2k = {0};
  y2k.tm_year = d_year - 1900; y2k.tm_mon = d_month-1; y2k.tm_mday = d_day;
  y2k.tm_hour = d_hour;   y2k.tm_min = d_min; y2k.tm_sec = d_sec;
  return difftime(mktime(&y2k),0);
}

lstanabase::~lstanabase(){
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

void lstanabase::Loop(TString histOut){
}

Int_t lstanabase::GetEntry(Long64_t entry){
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t lstanabase::LoadTree(Long64_t entry){
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
    Notify();
  }
  return centry;
}

void lstanabase::Init(TTree *tree){
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  //fChain->SetBranchAddress("evt", &evt, &b_evt);
  //fChain->SetBranchAddress("run", &run, &b_run);
  //fChain->SetBranchAddress("pValue", &pValue, &b_pValue);
  //...
  //...
  //
  //---------------------------------------------------
  // ADD HERE :

   fChain->SetBranchAddress("EventID", &EventID, &b_Event);
   fChain->SetBranchAddress("BunchXID", &BunchXID, &b_BunchXID);
   fChain->SetBranchAddress("NTotPhot", &NTotPhot, &b_NTotPhot);
   fChain->SetBranchAddress("Nhits", &Nhits, &b_Nhits);
   fChain->SetBranchAddress("primType", &primType, &b_primType);
   fChain->SetBranchAddress("primMomX", &primMomX, &b_primMomX);
   fChain->SetBranchAddress("primMomY", &primMomY, &b_primMomY);
   fChain->SetBranchAddress("primMomZ", &primMomZ, &b_primMomZ);
   fChain->SetBranchAddress("primPosX", &primPosX, &b_primPosX);
   fChain->SetBranchAddress("primPosY", &primPosY, &b_primPosY);
   fChain->SetBranchAddress("primPosZ", &primPosZ, &b_primPosZ);
   fChain->SetBranchAddress("primTime", &primTime, &b_primTime);
   fChain->SetBranchAddress("trigTopL", &trigTopL, &b_trigTopL);
   fChain->SetBranchAddress("trigTopEdep", &trigTopEdep, &b_trigTopEdep);
   fChain->SetBranchAddress("trigBotL", &trigBotL, &b_trigBotL);
   fChain->SetBranchAddress("trigBotEdep", &trigBotEdep, &b_trigBotEdep);
   fChain->SetBranchAddress("nPhot", &nPhot, &b_nPhot);
   fChain->SetBranchAddress("TrackID", TrackID, &b_TrackID);
   fChain->SetBranchAddress("ParentID", ParentID, &b_ParentID);
   fChain->SetBranchAddress("Energy", Energy, &b_Energy);
   fChain->SetBranchAddress("Wavelength", Wavelength, &b_Wavelength);
   fChain->SetBranchAddress("Time", Time, &b_Time);
   fChain->SetBranchAddress("photPathLen", photPathLen, &b_photPathLen);
   fChain->SetBranchAddress("SecID", SecID, &b_SecID);
   fChain->SetBranchAddress("chID", chID, &b_chID);
   fChain->SetBranchAddress("PosX", PosX, &b_PosX);
   fChain->SetBranchAddress("PosY", PosY, &b_PosY);
   fChain->SetBranchAddress("PosZ", PosZ, &b_PosZ);
   fChain->SetBranchAddress("MomX", MomX, &b_MomX);
   fChain->SetBranchAddress("MomY", MomY, &b_MomY);
   fChain->SetBranchAddress("MomZ", MomZ, &b_MomZ);
   fChain->SetBranchAddress("trkMomX", trkMomX, &b_trkMomX);
   fChain->SetBranchAddress("trkMomY", trkMomY, &b_trkMomY);
   fChain->SetBranchAddress("trkMomZ", trkMomZ, &b_trkMomZ);
   fChain->SetBranchAddress("trkPosX", trkPosX, &b_trkPosX);
   fChain->SetBranchAddress("trkPosY", trkPosY, &b_trkPosY);
   fChain->SetBranchAddress("trkPosZ", trkPosZ, &b_trkPosZ);
   fChain->SetBranchAddress("trkT", trkT, &b_trkT);
   fChain->SetBranchAddress("trkLength", trkLength, &b_trkLength);
   
  //---------------------------------------------------
  Notify();
}

Bool_t lstanabase::Notify(){
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.
  return kTRUE;
}

void lstanabase::Show(Long64_t entry){
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t lstanabase::Cut(Long64_t entry){
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}
