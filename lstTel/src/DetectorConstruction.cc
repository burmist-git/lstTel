//My
#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"

//G4
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Para.hh"
#include "G4Paraboloid.hh"
#include "G4EllipticalTube.hh"
#include "G4ExtrudedSolid.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4VisAttributes.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Color.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"
#include "globals.hh"
//magnetic field
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4ChordFinder.hh"
#include "G4UserLimits.hh"
//GDML
//#include <G4GDMLParser.hh>

//root 
#include "TMath.h"

DetectorConstruction::DetectorConstruction()
{
  //magField = new MagneticField();
  worldVisAtt = new G4VisAttributes();
  quartzVisAtt = new G4VisAttributes();
  aerogelVisAtt = new G4VisAttributes();
  sensitiveVisAtt = new G4VisAttributes();
  pmtboxVisAtt = new G4VisAttributes();
  absVisAtt = new G4VisAttributes();
  // Define Materials to be used
  DefineMaterials();
}

DetectorConstruction::~DetectorConstruction()
{
  //delete magField;
  delete worldVisAtt;
  delete quartzVisAtt;
  delete sensitiveVisAtt;
  delete pmtboxVisAtt;
  delete absVisAtt;
  delete stepLimit;
}

void DetectorConstruction::DefineMaterials()
{
  G4String symbol;
  G4double a, z, density;
  G4int ncomponents, natoms;
  G4double fractionmass;

  // Define elements
  G4Element* H = 
    new G4Element("Hydrogen", symbol = "H", z = 1., a = 1.01*g/mole);
  G4Element* C = 
    new G4Element("Carbon",   symbol = "C", z = 6., a = 12.01*g/mole);
  G4Element* N = 
    new G4Element("Nitrogen", symbol = "N", z = 7., a = 14.01*g/mole);
  G4Element* O =
    new G4Element("Oxygen",   symbol = "O", z = 8., a = 16.00*g/mole);
  G4Element* Si = 
    new G4Element("Silicon",  symbol = "Si", z = 14., a = 28.09*g/mole);
  G4Element* Al = 
    new G4Element("Aluminum", symbol = "Al", z = 13., a = 26.98*g/mole);
  G4Element* F = 
    new G4Element("Fluorine", symbol = "F", z = 9., a = 19.0*g/mole);

  // Quartz Material (SiO2_cladd)
  SiO2_cladd = new G4Material("quartzCladd", density = 2.200*g/cm3, ncomponents = 2);
  SiO2_cladd->AddElement(Si, natoms = 1);
  SiO2_cladd->AddElement(O , natoms = 2);

  // Quartz Material (SiO2_coat)
  SiO2_coat = new G4Material("quartzCoat", density = 2.200*g/cm3, ncomponents = 2);
  SiO2_coat->AddElement(Si, natoms = 1);
  SiO2_coat->AddElement(O , natoms = 2);

  // Quartz Material (SiO2)
  SiO2 = new G4Material("quartz", density = 2.200*g/cm3, ncomponents = 2);
  SiO2->AddElement(Si, natoms = 1);
  SiO2->AddElement(O , natoms = 2);

  C12H26 = new G4Material("mineraloil", density = 0.870*g/cm3, ncomponents = 2);
  C12H26->AddElement(C, natoms = 12);
  C12H26->AddElement(H, natoms = 26);
  
  // Aerogel Material (SiO2)
  Aerogel = new G4Material("Aerogel", density = 2000*g/m3, ncomponents = 2);
  Aerogel->AddElement(Si, natoms = 1);
  Aerogel->AddElement(O , natoms = 2);

  // C4F10
  C4F10 = new G4Material("fluorocarbon", density = 1.8/1000*g/cm3, ncomponents = 2);
  C4F10->AddElement(C, natoms = 4);
  C4F10->AddElement(F, natoms = 10);

  // Air
  //Air = new G4Material("Air", density = 1.290*mg/cm3, ncomponents = 2);
  Air = new G4Material("Air", density = 1.000*mg/cm3, ncomponents = 2);
  Air->AddElement(N, fractionmass = 0.7);
  Air->AddElement(O, fractionmass = 0.3);
  // Air Zero
  AirZero = new G4Material("AirZero", density = 0.010*mg/cm3, ncomponents = 2);
  AirZero->AddElement(N, fractionmass = 0.7);
  AirZero->AddElement(O, fractionmass = 0.3);

  // Aluminum
  //Aluminum = new G4Material("Aluminum", density = 2.7*g/cm3, ncomponents = 1);
  Aluminum = new G4Material("Aluminum", density = 2700*g/cm3, ncomponents = 1);
  Aluminum->AddElement(Al, fractionmass = 1.0);

  AluminumZero = new G4Material("AluminumZero", density = 0.0000027*g/cm3, ncomponents = 1);
  AluminumZero->AddElement(Al, fractionmass = 1.0);

  // Aluminum of the mirror
  AluminumMirr = new G4Material("mirrAluminum", density = 2.7*g/cm3, ncomponents = 1);
  AluminumMirr->AddElement(Al, fractionmass = 1.0);

  /*
  // Assign Materials
  world.material = Air;
  secA.material = SiO2;
  secB.material = SiO2;
  secC.material = SiO2;
  secWin.material = SiO2;
  fiberCorr.material = SiO2;
  fiberClad.material = SiO2_cladd;
  //fiberCoat.material = SiO2_coat;
  fiberCoat.material = Aluminum;
  sensitive.material = Aluminum;
  abs1.material = Aluminum;
  */

  //
  // Generate and Add Material Properties Table
  //						
  const G4int num = 36;
  G4double WaveLength[num];
  G4double Absorption[num];      // Default value for absorption
  G4double AirAbsorption[num];
  G4double AirRefractiveIndex[num];
  G4double PhotonEnergy[num];
  G4double AbsorptionDummy[num];
  G4double RefractiveIndexDummy[num];

  // Absorption of quartz per 1m
  G4double QuartzAbsorption[num] =
    {0.999572036,0.999544661,0.999515062,0.999483019,0.999448285,
     0.999410586,0.999369611,0.999325013,0.999276402,0.999223336,
     0.999165317,0.999101778,0.999032079,0.998955488,0.998871172,
     0.998778177,0.99867541 ,0.998561611,0.998435332,0.998294892,
     0.998138345,0.997963425,0.997767484,0.997547418,0.99729958 ,
     0.99701966 ,0.99670255 ,0.996342167,0.995931242,0.995461041,
     0.994921022,0.994298396,0.993577567,0.992739402,0.991760297,
     0.990610945};

  G4double C4F10RefractiveIndex[num];
  G4double AerogelRefractiveIndex[num];

  for (int i=0; i<num; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    // Aerogel and C4F10
    Absorption[i] = 100*m;
    // Air
    AirAbsorption[i] = 4000.0*m;
    //AirAbsorption[i] = 4.0*mm;
    AirRefractiveIndex[i] = 1.0 + 0.22*1.0e-3;
    AbsorptionDummy[i] = 0.0*mm;
    //RefractiveIndexDummy[i] = AirRefractiveIndex[i];
    RefractiveIndexDummy[i] = 0.0;
    // C4F10
    C4F10RefractiveIndex[i] = 1.0014;
    // Aerogel
    AerogelRefractiveIndex[i] = 1.03;
    PhotonEnergy[num - (i+1)] = twopi*hbarc/WaveLength[i];
    // Absorption is given per length and G4 needs mean free path
    // length, calculate it here
    // mean free path length - taken as probablility equal 1/e
    // that the photon will be absorbed
    QuartzAbsorption[i] = (-1)/log(QuartzAbsorption[i])*100*cm;
    //QuartzAbsorption[i] = 4.0*cm;
    //EpotekAbsorption[i] = (-1)/log(EpotekAbsorption[i])*
    //epotekBarJoint.thickness;
  }

  G4double QuartzRefractiveIndex[num] =
    {1.456535,1.456812,1.4571  ,1.457399,1.457712,1.458038,
     1.458378,1.458735,1.459108,1.4595  ,1.459911,1.460344,
     1.460799,1.46128 ,1.461789,1.462326,1.462897,1.463502,
     1.464146,1.464833,1.465566,1.46635 ,1.46719 ,1.468094,
     1.469066,1.470116,1.471252,1.472485,1.473826,1.475289,
     1.476891,1.478651,1.480592,1.482739,1.485127,1.487793};

  G4double CladdingRefractiveIndex[num];

  for(int i=0; i<num; i++){
    CladdingRefractiveIndex[i] = TMath::Sqrt(QuartzRefractiveIndex[i]*QuartzRefractiveIndex[i]-0.22*0.22); 
  }

  // Assign absorption and refraction to materials

  // Quartz
  G4MaterialPropertiesTable* QuartzMPT = new G4MaterialPropertiesTable();
  QuartzMPT->AddProperty("RINDEX", PhotonEnergy, QuartzRefractiveIndex, num);
  QuartzMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);
  
  // C4F10
  G4MaterialPropertiesTable* C4F10MPT = new G4MaterialPropertiesTable();
  C4F10MPT->AddProperty("RINDEX", PhotonEnergy, C4F10RefractiveIndex, num);
  //C4F10MPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption, num);
  C4F10MPT->AddProperty("ABSLENGTH", PhotonEnergy, AbsorptionDummy, num);
  
  // Aerogel
  G4MaterialPropertiesTable* AerogelMPT = new G4MaterialPropertiesTable();
  AerogelMPT->AddProperty("RINDEX", PhotonEnergy, AerogelRefractiveIndex, num);
  AerogelMPT->AddProperty("ABSLENGTH", PhotonEnergy, Absorption, num);
  
  // Cladding (of the fiber) only for the fiber aplication
  G4MaterialPropertiesTable* CladdingMPT = new G4MaterialPropertiesTable();
  CladdingMPT->AddProperty("RINDEX", PhotonEnergy, CladdingRefractiveIndex, num);
  CladdingMPT->AddProperty("ABSLENGTH", PhotonEnergy, QuartzAbsorption, num);

  // Air
  G4MaterialPropertiesTable* AirMPT = new G4MaterialPropertiesTable();
  AirMPT->AddProperty("RINDEX", PhotonEnergy, AirRefractiveIndex, num);
  AirMPT->AddProperty("ABSLENGTH", PhotonEnergy, AirAbsorption, num);

  // Dummy
  G4MaterialPropertiesTable* DummyMPT = new G4MaterialPropertiesTable();
  DummyMPT->AddProperty("RINDEX", PhotonEnergy, RefractiveIndexDummy, num);
  DummyMPT->AddProperty("ABSLENGTH", PhotonEnergy, AbsorptionDummy, num);

  // Assign this material to the bars
  SiO2->SetMaterialPropertiesTable(QuartzMPT);
  C12H26->SetMaterialPropertiesTable(QuartzMPT);  
  SiO2_cladd->SetMaterialPropertiesTable(CladdingMPT);
  C4F10->SetMaterialPropertiesTable(C4F10MPT);
  Aerogel->SetMaterialPropertiesTable(AerogelMPT);
  Air->SetMaterialPropertiesTable(AirMPT);
  //AluminumZero->SetMaterialPropertiesTable(DummyMPT);
  //Aluminum->SetMaterialPropertiesTable(C4F10MPT);
  
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  bool overlapsChecking = false;
  bool only_senstaive = true;
  bool build_camera_frame = false;
  
  //
  //
  G4double world_sizeX =  500.0*m;
  G4double world_sizeY =  500.0*m;
  G4double world_sizeZ = 1500.0*m;
  //
  //G4double world_AirZero_sizeX = world_sizeX - 1.0*m;
  //G4double world_AirZero_sizeY = world_sizeY - 1.0*m;
  //G4double world_AirZero_sizeZ = world_sizeZ/2.0 - 5.0*m;
  //
  G4double world_AirZero_X0 = 0.0;
  G4double world_AirZero_Y0 = 0.0;
  G4double world_AirZero_Z0 = -world_sizeZ/4.0;
  //
  //

  G4double tel_Z_shift = -world_sizeZ/2.0 + 10.0*m;
  //G4double tel_Z_shift = 0.0;
  
  //
  //G4double parabolicMirror_f  = 29.30565*m;
  G4double parabolicMirror_f  = 28.00*m;
  G4double parabolicMirror_R1 = 0.0*m;
  G4double parabolicMirror_R2 = 12.393*m;
  G4double parabolicMirror_H  = parabolicMirror_R2*parabolicMirror_R2/4.0/parabolicMirror_f;
  G4double parabolicMirror_X0 = 0.0;
  G4double parabolicMirror_Y0 = 0.0;
  //G4double parabolicMirror_Z0 = 0.0*cm + tel_Z_shift;
  G4double parabolicMirror_Z0 = 0.0*cm;
  //
  G4double parabolicMirrorSub_R1 = parabolicMirror_R1;
  G4double parabolicMirrorSub_H  = parabolicMirror_H + 20.0*cm;
  G4double parabolicMirrorSub_R2 = 2.0*TMath::Sqrt(parabolicMirrorSub_H*parabolicMirror_f);  
  //
  G4double parabolicMirrorFrame_minimum_wall_thickness = 5.0*mm;
  G4double parabolicMirrorFrame_sizeX = parabolicMirror_R2*2 + parabolicMirrorFrame_minimum_wall_thickness*2.0;
  G4double parabolicMirrorFrame_sizeY = parabolicMirrorFrame_sizeX;
  G4double parabolicMirrorFrame_sizeZ = parabolicMirror_H + parabolicMirrorFrame_minimum_wall_thickness;
  //  
  G4double sensitive_sizeX = 3.1400*m;
  G4double sensitive_sizeY = 3.1400*m;
  G4double sensitive_sizeZ = 1.0*mm;
  G4double sensitive_X0 = 0.0*cm;
  G4double sensitive_Y0 = 0.0*cm;
  G4double sensitive_Z0 = 1.0/(1.0/parabolicMirror_f - 1.0/(12000*m)) - parabolicMirror_H/2.0 + tel_Z_shift;
  //  
  G4double camear_box_sizeX = 3.1400*m;
  G4double camear_box_sizeY = 3.1400*m;
  G4double camear_box_sizeZ = 1.7024*m;
  G4double camear_box_X0 = sensitive_X0;
  G4double camear_box_Y0 = sensitive_Y0;
  G4double camear_box_Z0 = sensitive_Z0 + sensitive_sizeZ/2.0 + 5.0*mm + camear_box_sizeZ/2.0;
  //
  G4double sensitive_full_sizeX = world_sizeX - 20*cm;
  G4double sensitive_full_sizeY = world_sizeY - 20*cm;
  G4double sensitive_full_sizeZ = 1.0*mm;
  G4double sensitive_full_X0 = 0.0*cm;
  G4double sensitive_full_Y0 = 0.0*cm;
  G4double sensitive_full_Z0 = -world_sizeZ/2.0 + sensitive_full_sizeZ/2.0 + 20*cm;

  
  //
  G4double mirror_shield_solid_sizeZ = 10.0*cm;
  G4double mirror_shield_X0 = 0.0;
  G4double mirror_shield_Y0 = 0.0;
  G4double mirror_shield_Z0 = parabolicMirror_Z0 + parabolicMirror_H/2.0 + mirror_shield_solid_sizeZ/2.0 + 10*cm;
  //
  G4double parabolicMirror_env_sizeX = parabolicMirror_R2*2.0*TMath::Sqrt(2.0) + 30*cm;
  G4double parabolicMirror_env_sizeY = parabolicMirror_R2*2.0*TMath::Sqrt(2.0) + 30*cm;
  G4double parabolicMirror_env_sizeZ = parabolicMirror_H + 1.0*m;
  G4double parabolicMirror_env_X0 = 0.0;
  G4double parabolicMirror_env_Y0 = 0.0;
  G4double parabolicMirror_env_Z0 = 0.0 + tel_Z_shift;
  //
  
  /*
  //  
  G4double sensitiveShield_sizeX = mineralOil_sizeX-1.0*mm;
  G4double sensitiveShield_sizeY = mineralOil_sizeY-1.0*mm;
  G4double sensitiveShield_sizeZ = 2*mm;
  G4double sensitiveShield_X0 = 0.0*cm;
  G4double sensitiveShield_Y0 = 0.0*cm;
  G4double sensitiveShield_Z0 = sensitive_Z0 + sensitiveShield_sizeZ/2.0 + sensitive_sizeZ/2.0;
  */

  //
  // Define World Volume
  //
  G4VSolid *world_solid = new G4Box("World",world_sizeX/2.0,world_sizeY/2.0,world_sizeZ/2.0);
  G4LogicalVolume *world_logical = new G4LogicalVolume(world_solid,Air,"World");
  G4VPhysicalVolume *world_physical = new G4PVPlacement(0,G4ThreeVector(),world_logical,"World",0,false,0);

  //
  //G4VSolid *world_AirZero_solid = new G4Box("World_AirZero",world_AirZero_sizeX/2.0,world_AirZero_sizeY/2.0,world_AirZero_sizeZ/2.0);
  //G4LogicalVolume *world_AirZero_logical = new G4LogicalVolume(world_AirZero_solid,AirZero,"World_AirZero");
  //new G4PVPlacement(0,G4ThreeVector(world_AirZero_X0,world_AirZero_Y0,world_AirZero_Z0),world_AirZero_logical,"World_AirZero",world_logical,false,0);

  //
  // sensitive
  //  
  G4RotationMatrix Ra_sensitive;
  G4ThreeVector Ta_sensitive;
  G4Transform3D Tr_sensitive;
  G4VSolid *sensitive_solid;
  G4LogicalVolume *sensitive_logical;
  if (!only_senstaive){
    sensitive_solid = new G4Box("Sensitive", sensitive_sizeX/2.0, sensitive_sizeY/2.0, sensitive_sizeZ/2.0);
    //sensitive_logical = new G4LogicalVolume(sensitive_solid, AirZero, "Sensitive");
    //sensitive_logical = new G4LogicalVolume(sensitive_solid, C12H26, "Sensitive");
    sensitive_logical = new G4LogicalVolume(sensitive_solid, C4F10, "Sensitive");
    Ta_sensitive.setX(sensitive_X0);
    Ta_sensitive.setY(sensitive_Y0);
    Ta_sensitive.setZ(sensitive_Z0);
    Tr_sensitive = G4Transform3D(Ra_sensitive, Ta_sensitive);
    new G4PVPlacement(Tr_sensitive,       //Transformation
		      sensitive_logical,  //its logical volume				 
		      "Sensitive",        //its name
		      world_logical,      //its mother  volume
		      false,	        //no boolean operation
		      0,                  //copy number
		      overlapsChecking);  //overlapsChecking
  }
  else{
    sensitive_solid = new G4Box("Sensitive", sensitive_full_sizeX/2.0, sensitive_full_sizeY/2.0, sensitive_full_sizeZ/2.0);
    sensitive_logical = new G4LogicalVolume(sensitive_solid, C4F10, "Sensitive");
    Ta_sensitive.setX(sensitive_full_X0);
    Ta_sensitive.setY(sensitive_full_Y0);
    Ta_sensitive.setZ(sensitive_full_Z0);
    Tr_sensitive = G4Transform3D(Ra_sensitive, Ta_sensitive);
    new G4PVPlacement(Tr_sensitive,       //Transformation
		      sensitive_logical,  //its logical volume				 
		      "Sensitive",        //its name
		      world_logical,      //its mother  volume
		      false,	        //no boolean operation
		      0,                  //copy number
		      overlapsChecking);  //overlapsChecking
  }
  
  //
  // camear_box
  //
  G4RotationMatrix Ra_camear_box;
  G4ThreeVector Ta_camear_box;
  G4Transform3D Tr_camear_box;
  G4VSolid *camear_box_solid = new G4Box("camear_box_solid", camear_box_sizeX/2.0, camear_box_sizeY/2.0, camear_box_sizeZ/2.0);
  G4LogicalVolume *camear_box_logical = new G4LogicalVolume(camear_box_solid, AluminumZero, "camear_box_logical");
  Ta_camear_box.setX(camear_box_X0);
  Ta_camear_box.setY(camear_box_Y0);
  Ta_camear_box.setZ(camear_box_Z0);
  Tr_camear_box = G4Transform3D(Ra_camear_box, Ta_camear_box);
  if(build_camera_frame)
    new G4PVPlacement(Tr_camear_box,      //Transformation
		      camear_box_logical, //its logical volume
		      "camear_box",       //its name
		      world_logical,      //its mother volume
		      false,	        //no boolean operation
		      0,                  //copy number
		      overlapsChecking);  //overlapsChecking


  //
  // Parabolic mirror envelop
  //
  G4RotationMatrix Ra_parabolicMirror_env;
  G4ThreeVector Ta_parabolicMirror_env;
  G4Transform3D Tr_parabolicMirror_env;
  G4VSolid *parabolicMirror_env_solid = new G4Box("parabolicMirror_env", parabolicMirror_env_sizeX/2.0, parabolicMirror_env_sizeY/2.0, parabolicMirror_env_sizeZ/2.0);
  G4LogicalVolume *parabolicMirror_env_logical = new G4LogicalVolume(parabolicMirror_env_solid,Air,"parabolicMirror_env");
  Ta_parabolicMirror_env.setX(parabolicMirror_env_X0);
  Ta_parabolicMirror_env.setY(parabolicMirror_env_Y0);
  Ta_parabolicMirror_env.setZ(parabolicMirror_env_Z0);
  Tr_parabolicMirror_env = G4Transform3D(Ra_parabolicMirror_env, Ta_parabolicMirror_env);
  if (!only_senstaive)
    new G4PVPlacement(Tr_parabolicMirror_env,      //Transformation
		      parabolicMirror_env_logical, //its logical volume
		      "parabolicMirror_env",       //its name
		      world_logical,               //its mother volume
		      false,	                 //no boolean operation
		      0,                           //copy number
		      overlapsChecking);           //overlapsChecking

  
  
  //
  // Parabolic mirror
  //
  G4Box* parabolicMirrorFrame_solid = new G4Box("parabolicMirrorFrame", parabolicMirrorFrame_sizeX/2.0, parabolicMirrorFrame_sizeY/2.0, parabolicMirrorFrame_sizeZ/2.0);
  G4Paraboloid* parabolicMirrorSub_solid = new G4Paraboloid("parabolicMirrorSub", parabolicMirrorSub_H/2.0, parabolicMirrorSub_R1, parabolicMirrorSub_R2);
  G4RotationMatrix* rotMatrix = new G4RotationMatrix();
  G4ThreeVector transVector(0.0,
			    0.0,
			    parabolicMirrorFrame_minimum_wall_thickness-(parabolicMirrorFrame_sizeZ - parabolicMirrorSub_H)/2.0);
  G4SubtractionSolid *parabolicMirror_solid = new G4SubtractionSolid("parabolicMirror_solid",parabolicMirrorFrame_solid,parabolicMirrorSub_solid,rotMatrix,transVector);
  G4LogicalVolume *parabolicMirror_logical = new G4LogicalVolume(parabolicMirror_solid, Aluminum,"parabolicMirror");
  G4RotationMatrix Ra_parabolicMirror;
  G4ThreeVector Ta_parabolicMirror;
  G4Transform3D Tr_parabolicMirror;
  Ta_parabolicMirror.setX(parabolicMirror_X0);
  Ta_parabolicMirror.setY(parabolicMirror_Y0);
  Ta_parabolicMirror.setZ(parabolicMirror_Z0);
  Tr_parabolicMirror = G4Transform3D(Ra_parabolicMirror, Ta_parabolicMirror);
  if (!only_senstaive)
    new G4PVPlacement(Tr_parabolicMirror,      //Transformation
		      parabolicMirror_logical, //its logical volume				 
		      "parabolicMirror",       //its name
		      parabolicMirror_env_logical, //its mother  volume (old world_logical)
		      false,	             //no boolean operation
		      0,                       //copy number
		      overlapsChecking);       //overlapsChecking
  

  G4VSolid *mirror_shield_solid = new G4Tubs("mirror_shield_solid", //name
					     parabolicMirror_R2,                  //pRMin
					     parabolicMirror_R2*TMath::Sqrt(2.0), //pRMax
					     mirror_shield_solid_sizeZ,           //pDz
					     0,                                   //pSPhi
					     twopi);                              //pDPhi
  G4LogicalVolume *mirror_shield_logical = new G4LogicalVolume(mirror_shield_solid, AluminumZero,"mirror_shield_logical");
  G4RotationMatrix Ra_mirror_shield;
  G4ThreeVector Ta_mirror_shield;
  G4Transform3D Tr_mirror_shield;
  Ta_mirror_shield.setX(mirror_shield_X0);
  Ta_mirror_shield.setY(mirror_shield_Y0);
  Ta_mirror_shield.setZ(mirror_shield_Z0);
  Tr_mirror_shield = G4Transform3D(Ra_mirror_shield, Ta_mirror_shield);
  if (!only_senstaive)
    new G4PVPlacement(Tr_mirror_shield,        //Transformation
		      mirror_shield_logical,   //its logical volume				 
		      "mirror_shield_logical", //its name
		      parabolicMirror_env_logical, //its mother volume (old world_logical)
		      false,	             //no boolean operation
		      0,                       //copy number
		      overlapsChecking);       //overlapsChecking

  
  //
  // Set Visualization Attributes
  //
  G4Color blue        = G4Color(0., 0., 1.);
  G4Color green       = G4Color(0., 1., 0.);
  G4Color red         = G4Color(1., 0., 0.);
  G4Color white       = G4Color(1., 1., 1.);
  G4Color DircColor   = G4Color(0.0, 0.0, 1.0, 0.2);
  G4Color SensColor   = G4Color(0.0, 1.0, 1.0, 0.1);

  worldVisAtt->SetColor(white);
  worldVisAtt->SetVisibility(true);
  //quartzVisAtt->SetColor(DircColor);
  //quartzVisAtt->SetVisibility(true);  
  sensitiveVisAtt->SetColor(SensColor);
  sensitiveVisAtt->SetVisibility(true);

  sensitive_logical->SetVisAttributes(sensitiveVisAtt);
  world_logical->SetVisAttributes(worldVisAtt);

  //sphericalMirror_body_logical->SetVisAttributes(quartzVisAtt);
  //flatMirror_logical->SetVisAttributes(quartzVisAtt);
  //al_body_logical->SetVisAttributes(absVisAtt);
  //c4f10_body_logical->SetVisAttributes(absVisAtt);
  //aerogel_body_logical->SetVisAttributes(quartzVisAtt);


  //
  // Define Optical Borders
  //
  // Define mirror surface
  const G4int num2 = 36;
  //G4double EfficiencyMirrors[num2];
  G4double WaveLength[num2];
  G4double PhotonEnergy[num2];
  //G4double MirrorReflectivity[num2];
  for (G4int i=0; i<num2; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    PhotonEnergy[num2 - (i+1)] = twopi*hbarc/WaveLength[i];
    //EfficiencyMirrors[i] = 0.0;
    //MirrorReflectivity[i]=0.85;
  }
  //
  G4double MirrorReflectivity[num2]=
    {0.87,0.88,0.885,0.89,0.895,0.9,0.905,0.91,0.915,0.92,0.923,0.9245,
     0.926,0.928,0.93,0.935,0.936,0.937,0.938,0.94,0.94,0.939,0.9382,
     0.938,0.937,0.937,0.936,0.935,0.934,0.932,0.93,0.928,0.926,0.924,
     0.922,0.92};
  //
  G4MaterialPropertiesTable* MirrorMPT = new G4MaterialPropertiesTable();
  MirrorMPT->AddProperty("REFLECTIVITY", PhotonEnergy, MirrorReflectivity, num2);
  //MirrorMPT->AddProperty("EFFICIENCY" , PhotonEnergy, EfficiencyMirrors,  num2);

  G4OpticalSurface* OpMirrorSurface = new G4OpticalSurface("MirrorSurface");
  OpMirrorSurface->SetType(dielectric_metal);
  OpMirrorSurface->SetFinish(polished);
  OpMirrorSurface->SetModel(glisur);
  OpMirrorSurface->SetMaterialPropertiesTable(MirrorMPT);

  new G4LogicalSkinSurface("MirrorSurfT",
  			   parabolicMirror_logical, OpMirrorSurface);
  ///////////////////////////////////////


  
  /*

  // Surface for killing photons at borders
  const G4int num1 = 2;
  G4double Ephoton[num1] = {1.5*eV, 5.8*eV};

  G4OpticalSurface* OpVolumeKillSurface =
    new G4OpticalSurface("VolumeKillSurface");
  OpVolumeKillSurface->SetType(dielectric_metal);
  OpVolumeKillSurface->SetFinish(polished);
  OpVolumeKillSurface->SetModel(glisur);


  G4double ReflectivityKill[num1] = {0., 0.};
  G4double EfficiencyKill[num1] = {1., 1.};
  G4MaterialPropertiesTable* VolumeKill = new G4MaterialPropertiesTable();
  VolumeKill->AddProperty("REFLECTIVITY", Ephoton, ReflectivityKill, num1);
  VolumeKill->AddProperty("EFFICIENCY",   Ephoton, EfficiencyKill,   num1);
  OpVolumeKillSurface->SetMaterialPropertiesTable(VolumeKill);
  new G4LogicalSkinSurface("SensitiveSurface", 
  			   sensitive_logical, OpVolumeKillSurface);
  



  // Define mirror surface
  const G4int num2 = 36;
  G4double EfficiencyMirrors[num2];
  G4double WaveLength[num2];
  G4double PhotonEnergy[num2];
  //G4double MirrorReflectivity[num2];
  for (G4int i=0; i<num2; i++) {
    WaveLength[i] = (300 + i*10)*nanometer;
    PhotonEnergy[num2 - (i+1)] = twopi*hbarc/WaveLength[i];
    EfficiencyMirrors[i] = 0.0;
    //MirrorReflectivity[i]=0.85;
  }
  //
  //G4double MirrorReflectivity[num2]=
  //  {0.87,0.88,0.885,0.89,0.895,0.9,0.905,0.91,0.915,0.92,0.923,0.9245,
  //   0.926,0.928,0.93,0.935,0.936,0.937,0.938,0.94,0.94,0.939,0.9382,
  //   0.938,0.937,0.937,0.936,0.935,0.934,0.932,0.93,0.928,0.926,0.924,
  //   0.922,0.92};
  //
  G4MaterialPropertiesTable* MirrorMPT = new G4MaterialPropertiesTable();
  //MirrorMPT->AddProperty("REFECTIVITY", PhotonEnergy, MirrorReflectivity, num2);
  MirrorMPT->AddProperty("EFFICIENCY" , PhotonEnergy, EfficiencyMirrors,  num2);

  G4OpticalSurface* OpMirrorSurface = new G4OpticalSurface("MirrorSurface");
  OpMirrorSurface->SetType(dielectric_metal);
  OpMirrorSurface->SetFinish(polished);
  OpMirrorSurface->SetModel(glisur);
  new G4LogicalSkinSurface("MirrorSurfT",
			   sphericalMirror_body_logical, OpMirrorSurface);
  new G4LogicalSkinSurface("MirrorSurfT",
			   flatMirror_logical, OpMirrorSurface);

  OpMirrorSurface->SetMaterialPropertiesTable(MirrorMPT);
  ///////////////////////////////////////

  */
  
  // 
  // Sensitive detector definition
  //
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  SensitiveDetector* aSD = new SensitiveDetector("fTOF");
  SDman->AddNewDetector(aSD);
  sensitive_logical->SetSensitiveDetector(aSD);
  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  G4double maxStep   = 10*mm;
  G4double maxLength = 10*mm;
  G4double maxTime   = 1.0*ns; 
  G4double minEkin   = 100.0*MeV;
  G4double mionRang  = 100.0*cm;
  stepLimit = new G4UserLimits(maxStep,maxLength,maxTime,minEkin,mionRang);
  parabolicMirror_logical->SetUserLimits(stepLimit);
  sensitive_logical->SetUserLimits(stepLimit);
  camear_box_logical->SetUserLimits(stepLimit);
  mirror_shield_logical->SetUserLimits(stepLimit);

  G4double maxStepFine   = 0.1*mm;
  G4double maxLengthFine = 2.0*m;
  G4double maxTimeFine   = 20.0*ns; 
  G4double minEkinFine   = 1.0*MeV;
  G4double mionRangFine  = 1.0*mm;
  stepLimitFine = new G4UserLimits(maxStepFine,maxLengthFine,maxTimeFine,minEkinFine,mionRangFine);
  parabolicMirror_env_logical->SetUserLimits(stepLimitFine);



  /*
  secA.logical->SetUserLimits(stepLimit);
  secB.logical->SetUserLimits(stepLimit);
  secC.logical->SetUserLimits(stepLimit);
  secWin.logical->SetUserLimits(stepLimit);
  */

  //G4GDMLParser parser;
  //parser.Write("CpFM.gdml", world.physical);

  return world_physical;
}
