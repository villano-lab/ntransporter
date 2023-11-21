// MinimalUserInitialization.cc - source code for MinimalUserInitialization.hh

#include "MinimalUserInitialization.hh"

#include "G4Neutron.hh"
#include "G4NistManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

// G4HadronicProcess objects
#include "G4HadronElasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronFissionProcess.hh"

// G4VCrossSectionDataSet objects
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPFissionData.hh"
#include "G4ParticleHPInelasticData.hh"

// G4HadronicInteraction objects
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPFission.hh"
#include "G4ParticleHPCapture.hh"


NeutronPhysicsList::NeutronPhysicsList() : G4VModularPhysicsList() {;}


void NeutronPhysicsList::ConstructHadronics() {
    
    G4PhysicsListHelper *ph = G4PhysicsListHelper::GetPhysicsListHelper();

    G4ParticleDefinition *theNeutron = G4Neutron::Definition();
    
    ph->RegisterProcess(new G4HadronElasticProcess(), theNeutron);
    ph->RegisterProcess(new G4HadronInelasticProcess("NeutronInelastic", theNeutron), theNeutron);
    ph->RegisterProcess(new G4HadronCaptureProcess(), theNeutron);
    ph->RegisterProcess(new G4HadronFissionProcess(), theNeutron);

    G4ProcessVector *procs = theNeutron->GetProcessManager()->GetProcessList();

    static_cast<G4HadronicProcess*>((*procs)[1])->AddDataSet(new G4ParticleHPElasticData);
    static_cast<G4HadronicProcess*>((*procs)[2])->AddDataSet(new G4ParticleHPInelasticData);
    static_cast<G4HadronicProcess*>((*procs)[3])->AddDataSet(new G4ParticleHPCaptureData);
    static_cast<G4HadronicProcess*>((*procs)[4])->AddDataSet(new G4ParticleHPFissionData);

    static_cast<G4HadronicProcess*>((*procs)[1])->RegisterMe(new G4ParticleHPElastic);
    static_cast<G4HadronicProcess*>((*procs)[2])->RegisterMe(new G4ParticleHPInelastic);
    static_cast<G4HadronicProcess*>((*procs)[3])->RegisterMe(new G4ParticleHPCapture);
    static_cast<G4HadronicProcess*>((*procs)[4])->RegisterMe(new G4ParticleHPFission);
}




void NeutronPhysicsList::ConstructProcess() {
    AddTransportation();
    ConstructHadronics();
}


void NeutronPhysicsList::ConstructParticle() {
    G4ParticleDefinition *theNeutron = G4Neutron::NeutronDefinition();
    G4ParticleTable::GetParticleTable()->SetGenericIon(theNeutron);
}


MinimalDetector::MinimalDetector() : G4VUserDetectorConstruction() {;}

G4VPhysicalVolume* MinimalDetector::Construct() {
    G4Material *mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
    G4double l = 10.*m;
    G4Box *worldBox = new G4Box("WorldBox", l,l,l);
    return new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), new G4LogicalVolume(worldBox, mat, "World"), "World", nullptr, false, 0, false);
}



MinimalInitialization::MinimalInitialization() : G4VUserActionInitialization() {;}

void MinimalInitialization::Build() const {
    SetUserAction(new MinimalPrimaryGenerator());
}


MinimalPrimaryGenerator::MinimalPrimaryGenerator() : G4VUserPrimaryGeneratorAction() {;}



void MinimalPrimaryGenerator::GeneratePrimaries(G4Event *anEvent) {
    return;
}


