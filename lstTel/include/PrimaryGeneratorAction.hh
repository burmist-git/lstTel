#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4String.hh"

class G4ParticleGun;
class G4Event;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction();
  
public:
  void GeneratePrimaries(G4Event* anEvent);
  
public:
  void SetParticleName(G4String particleName) {_particleName = particleName;}

  G4String GetParticleName() {return _particleName;}
  G4double GetParticleMomentum() {return _particleMomentum;}
  G4double GetPhiAngle() { return _PhiAngle; }
  G4double GetThetaAngle() { return _ThetaAngle; }

  void SetParticleMomentum(G4double momentum) {_particleMomentum = momentum;}
  void SetPhiAngle(G4double val) {_PhiAngle = val;}
  void SetThetaAngle(G4double val) {_ThetaAngle = val;}

  void SetDisp(G4double valx, G4double valy) { _x_disp_cm = valx; _y_disp_cm = valy;}

  
  void SetSinglePhoton(G4bool singlePhoton) { _singlePhoton = singlePhoton;}
  G4bool SinglePhotonGenerator() { return _singlePhoton;}
  G4int GetBunchXID(){return _BunchXID;};
private:
  G4ParticleGun* _particleGun;
  G4String _particleName;
  G4double _particleMomentum;
  G4double _PhiAngle;
  G4double _ThetaAngle;
  G4double _x_disp_cm;
  G4double _y_disp_cm;

  G4bool _singlePhoton;

  G4int _BunchXID;

  G4int GenFlatInt(G4int iMin,G4int iMax);
  void generateThetaAndPhi();
  G4double genCos2dist();

  G4double genEdist( G4double particleMomentum_min, G4double particleMomentum_max);
  G4double genThetadist( G4double Theta_min, G4double Theta_max);
  void gen_XY_Intersection(G4double R, G4double &x, G4double &y);
  bool gen_XY_generation( G4double xIntersection, G4double yIntersection, G4double zIntersection,
			  G4double theta_cos_flat, G4double phi_flat, G4double zInit,
			  G4double &xInit, G4double &yInit);
  
};

#endif
