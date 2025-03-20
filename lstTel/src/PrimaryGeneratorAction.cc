//my
#include "PrimaryGeneratorAction.hh"
#include "VolumeStructures.hh"

//G4
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"

//root
#include "TMath.h"
#include "TVector3.h"

using namespace CLHEP;

PrimaryGeneratorAction::PrimaryGeneratorAction() :
  _particleGun(0),
  _particleName("pi+"),
  _particleMomentum(3.0*GeV),
  _PhiAngle(0.0*deg),
  _ThetaAngle(0.0*deg),
  _singlePhoton(false)
{
  _particleGun = new G4ParticleGun(1);  
  _BunchXID = 0;

  //backGen = new backgroundGen();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete _particleGun;
  //delete backGen;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle;
  G4double xInit, yInit, zInit;
  G4double dX, dY, dZ;
  G4double Ekin, m_mass;
  G4ThreeVector dir;
  //
  Int_t generation_type;
  G4double particleMomentum_min = _particleMomentum;
  G4double particleMomentum_max = 200.0*GeV;
  G4double R_impact_max = 17.0*m;
  G4double theta_max = 2.0*deg;
  G4double power_m_two_dist;
  G4double theta_cos_flat;
  G4double phi_flat;
  G4double xIntersection;
  G4double yIntersection;
  G4double zIntersection;

  //
  //G4cout<<"_particleMomentum = "<<_particleMomentum<<G4endl
  //	<<"m                 = "<<m<<G4endl
  //	<<"_ThetaAngle       = "<<_ThetaAngle*180.0/TMath::Pi()<<G4endl
  //	<<"_PhiAngle         = "<<_PhiAngle*180.0/TMath::Pi()<<G4endl;
  //_ThetaAngle = _ThetaAngle + (-1 + 2*G4UniformRand())*2*TMath::Pi()/360;//mearing of one degree
  //_PhiAngle   = _PhiAngle   + (-1 + 2*G4UniformRand())*2*TMath::Pi()/360;//mearing of one degree
  //
  //------------------
  //generation_type = 0; // baseline
  generation_type = 1;   // Emin (_particleMomentum_min) - Emax (_particleMomentum_max) baseline 1/E^2, (0 - R_max), (0 - theta_max)
  //------------------
  //
  if(generation_type == 0){
    xInit = 0.0*cm;
    yInit = 0.0*cm;
    zInit = 730.0*m;
    //
    _BunchXID++;
    particle = particleTable->FindParticle(_particleName);
    m_mass = particle->GetPDGMass();
    Ekin = (TMath::Sqrt(_particleMomentum*_particleMomentum + m_mass*m_mass) - m_mass);
    //
    dX =  TMath::Sin(_ThetaAngle)*TMath::Cos(_PhiAngle);
    dY =  TMath::Sin(_ThetaAngle)*TMath::Sin(_PhiAngle);
    dZ =  TMath::Cos(_ThetaAngle);
    //
    dir.set(dX, dY, dZ);
    _particleGun->SetParticleDefinition(particle);
    _particleGun->SetParticleMomentumDirection(dir);
    _particleGun->SetParticleEnergy(Ekin);  
    _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
    _particleGun->GeneratePrimaryVertex(anEvent);
  }
  else if(generation_type == 1){
    xInit = 0.0*cm;
    yInit = 0.0*cm;
    zInit = 730.0*m;
    zIntersection = -749.8*m;
    //
    gen_XY_Intersection( R_impact_max, xIntersection, yIntersection);
    _BunchXID++;
    particle = particleTable->FindParticle(_particleName);
    m_mass = particle->GetPDGMass();
    //
    power_m_two_dist = genEdist( particleMomentum_min, particleMomentum_max);
    theta_cos_flat = genThetadist((_ThetaAngle - theta_max), _ThetaAngle);
    phi_flat = twopi*G4UniformRand();
    gen_XY_generation( xIntersection, yIntersection, zIntersection, theta_cos_flat, phi_flat, zInit, xInit, yInit);
    //    
    Ekin = (TMath::Sqrt(power_m_two_dist*power_m_two_dist + m_mass*m_mass) - m_mass);
    dX =  TMath::Sin(theta_cos_flat)*TMath::Cos(phi_flat);
    dY =  TMath::Sin(theta_cos_flat)*TMath::Sin(phi_flat);
    dZ =  TMath::Cos(theta_cos_flat);
    //
    //
    dir.set(dX, dY, dZ);
    _particleGun->SetParticleDefinition(particle);
    _particleGun->SetParticleMomentumDirection(dir);
    _particleGun->SetParticleEnergy(Ekin);  
    _particleGun->SetParticlePosition(G4ThreeVector(xInit, yInit, zInit));
    _particleGun->GeneratePrimaryVertex(anEvent);
  }  
}

G4int PrimaryGeneratorAction::GenFlatInt(G4int iMin,G4int iMax){
  G4int val;
  val = (G4int)floor((iMax - iMin + 1)*G4UniformRand() + iMin);
  return val;
}

void PrimaryGeneratorAction::generateThetaAndPhi(){
  _PhiAngle = G4UniformRand()*2*TMath::Pi();
  _ThetaAngle = TMath::Pi() - genCos2dist();
}

G4double PrimaryGeneratorAction::genCos2dist(){
  G4double theta = -999.0;//deg 
  G4double x = -999.0;
  G4double y = -999.0;
  while(theta==-999.0){
    x = G4UniformRand()*(70.0*TMath::Pi()/180.0); //rad
    y = G4UniformRand();
    if(TMath::Power(TMath::Cos(x),1.85)>y){
      theta = x;
    }
  }  
  return theta;
}

// 1.0/mom^2
G4double PrimaryGeneratorAction::genEdist( G4double particleMomentum_min, G4double particleMomentum_max){
  return 1.0/((1.0/particleMomentum_min - 1.0/particleMomentum_max)*G4UniformRand() + 1.0/particleMomentum_max);
}

// cos flat
G4double PrimaryGeneratorAction::genThetadist( G4double Theta_min, G4double Theta_max){
  return TMath::ACos((TMath::Cos(Theta_min) - TMath::Cos(Theta_max))*G4UniformRand() + TMath::Cos(Theta_max));
}

void PrimaryGeneratorAction::gen_XY_Intersection(G4double R, G4double &x, G4double &y){
  G4double rr = -999.0;
  while(rr == -999.0){
    x = 2*R*G4UniformRand() - R;
    y = 2*R*G4UniformRand() - R;
    rr = x*x + y*y;
    if(rr>=R*R)
      rr = -999.0;
  }
}

bool PrimaryGeneratorAction::gen_XY_generation( G4double xIntersection, G4double yIntersection, G4double zIntersection,
						G4double theta_cos_flat, G4double phi_flat, G4double zInit,
						G4double &xInit, G4double &yInit){
  //
  TVector3 trk_r0(xIntersection, yIntersection, zIntersection);
  TVector3 trk_v;
  trk_v.SetMagThetaPhi(1.0, theta_cos_flat, phi_flat);  
  TVector3 plane_r0(0.0, 0.0, zInit);
  TVector3 plane_n(0.0, 0.0, 1.0);
  //
  G4double div = trk_v.Dot(plane_n); 
  if(div == 0){
    if(plane_r0.Dot(plane_n) == trk_r0.Dot(plane_n)){ 
      xInit = xIntersection;
      yInit = yIntersection;
    }
    else{
      xInit = -999.0;
      yInit = -999.0;
    }
    return false;
  }
  G4double t = (plane_r0.Dot(plane_n) - trk_r0.Dot(plane_n))/div;
  TVector3 trk_gen = trk_v*t + trk_r0;
  xInit = trk_gen.x();
  yInit = trk_gen.y();
  //
  //G4cout<<"xIntersection  "<<xIntersection<<G4endl
  //	<<"yIntersection  "<<yIntersection<<G4endl
  //	<<"zIntersection  "<<zIntersection<<G4endl;
  //
  //G4cout<<"theta_cos_flat "<<theta_cos_flat<<G4endl
  //	<<"phi_flat       "<<phi_flat<<G4endl
  //	<<"theta_cos_flat "<<theta_cos_flat/deg<<G4endl
  //	<<"phi_flat       "<<phi_flat/deg<<G4endl;
  //
  //G4cout<<"t              "<<t<<G4endl
  //	<<"xInit          "<<xInit<<G4endl
  //	<<"yInit          "<<yInit<<G4endl
  //	<<"trk_gen.z()    "<<trk_gen.z()<<G4endl
  //	<<"zInit          "<<zInit<<G4endl;
  //
  return true;
}
