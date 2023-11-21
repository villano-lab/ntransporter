#ifndef MINIMALUSERINITIALIZATION_H
#define MINIMALUSERINITIALIZATION_H

// MinimalUserInitialization.hh - interface file for Geant4-required classes
//
// includes four classes:
//
//     NeutronPhysicsList      - G4VModularPhysicsList with NeutronHP data for 
//                               neutron interactions elastic, inelastic, 
//                               capture, and fission
//
//     MinimalDetector         - G4VUserDetectorConstruction that creates 
//                               worldbox of G4_Galactic
//
//     MinimalInitialization   - G4VUserActionInitialization that uses 
//                               MinimalPrimaryGenerator
//
//     MinimalPrimaryGenerator - G4VUserPrimaryGeneratorAction that implements 
//                               empty GeneratePrimaries() method
//


#include "G4VModularPhysicsList.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"


class NeutronPhysicsList : public G4VModularPhysicsList {
public:
    NeutronPhysicsList();
    ~NeutronPhysicsList() {;}
    void ConstructProcess();
    void ConstructHadronics();
    void ConstructParticle();
};


class MinimalDetector : public G4VUserDetectorConstruction {
public:
    MinimalDetector();
    ~MinimalDetector() {;}
    G4VPhysicalVolume* Construct();
};



class MinimalInitialization : public G4VUserActionInitialization {
public:
    MinimalInitialization();
    ~MinimalInitialization() {;}
    void Build() const;
};


class MinimalPrimaryGenerator : public G4VUserPrimaryGeneratorAction {
public:
    MinimalPrimaryGenerator();
    ~MinimalPrimaryGenerator() {;}
    void GeneratePrimaries(G4Event *);
};


#endif