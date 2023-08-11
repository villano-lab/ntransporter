// hypothetically, script to grab cross section data from Geant
// first attempt at script - grab for neutrons incident on silicon

#include <iostream>
#include <vector>

#include "G4MTRunManager.hh"
#include "G4UImanager.hh"

#include "CDMSGeometryManager.hh"

#include "SuperSim_Main.hh"

#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "CDMSMaterialTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"

#include "G4Neutron.hh"
#include "G4ParticleDefinition.hh"
#include "Shielding.hh"

#include "G4ProcessManager.hh"
#include "G4PDefManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcess.hh"


int main(int argc, char **argv) {
  
  // needed to initialize correct world geometry and default region
  //G4RunManager* rm = G4RunManager::GetRunManager();
  //G4Region *worldRegion = new G4Region("DefaultRegionForTheWorld");
  //worldRegion->SetProductionCuts(new G4ProductionCuts());
  //G4RegionStore::Register(worldRegion);
  
  //G4MTRunManager *runManager = new G4MTRunManager();

  //runManager->SetUserInitialization(new 
  //runManager->SetUserInitialization(new Shielding());

  //runManager->InitializePhysics();
  //G4UImanager *UI = G4UImanager::GetUIpointer();

  //CDMSGeometryManager *geomManager = CDMSGeometryManager::Instance();
  SuperSim_Main *sMain = new SuperSim_Main();

  sMain->SetVerboseLevel(0);
  sMain->runManager.SetVerboseLevel(0);
  sMain->Configure("full"); // full Shielding physics list, etc
  //CDMSRunManager runManager = sMain->runManager;

  //G4UImanager *UI = sMain->UImanager;

  sMain->UImanager->ApplyCommand("/CDMS/lab NoLab");
  sMain->UImanager->ApplyCommand("/CDMS/detector ZIP");

  sMain->runManager.Initialize();

  //if (UI) {
	//UI->ApplyCommand("/run/initialize");
  //} else {
	//std::cout << "nope on the UI" << std::endl;
  //}
  

  G4Neutron *theNeutron = G4Neutron::Definition();
  //theNeutron->SetParticleDefinitionID(1);
  const CDMSMaterialTable *theTable = CDMSMaterialTable::GetInstance();
  
  G4Material *rock = theTable->GetMaterial("Norite");
  G4int nElm = rock->GetNumberOfElements();
  const G4ElementVector *elmVector = rock->GetElementVector();
  const G4double *fracVector = rock->GetFractionVector();

  //G4NeutronHPElasticData *elastic;
  //G4NeutronHPCaptureData *capture;
  //G4NeutronHPInelasticData *inelastic;
  
  // "Shielding" list should be chosen by default in SuperSim_Main
  //Shielding *shield = new Shielding(verbose); // create Shielding physics list object  
  //shield->ConstructProcess(); // construct processes

  //typedef std::vector<G4PhysicsConstructor*> G4PhysConstVectorData;
  //G4PhysConstVectorData* physicsVector = 
  //                  shield->GetSubInstanceManager().offset[0].physicsVector;

  
  // get ProcessManager for the neutron
  G4ProcessManager *theMan = theNeutron->GetProcessManager();
  //G4ProcessManager *theMan = theNeutron->GetMasterProcessManager();

  //const G4PDefManager pdef = G4ParticleDefinition::GetSubInstanceManager();
  //G4PDefData dat = pdef.offset()[1];

  std::cout << "n ID = " << theNeutron->GetInstanceID() << std::endl;
  //std::cout << "pdef offset = " << dat << std::endl;
  //G4ParticleDefinition::GetSubInstanceManager().CreateSubInstance();
  //G4ProcessManager *theMan = (G4ParticleDefinition::GetSubInstanceManager().offset()[theNeutron->GetInstanceID()]).theProcessManager;

  //std::cout << (theMan ? "Y\n" : "N\n") << std::endl;
  
  // get element from material
  G4int nProc = theMan->GetProcessListLength();
  G4ProcessVector *processes = theMan->GetProcessList();
  const G4Element *elm = (*elmVector)[6]; // silicon
  
  // variables to store cross sections, neutron energy
  G4double x, Eneut = 2.*keV;


  std::cout << "Cross sections for " << elm->GetName() << " with " << Eneut 
            << " MeV neutron" << std::endl;

  G4DynamicParticle *dynamicNeutron = new G4DynamicParticle(theNeutron, 
                                              G4ThreeVector(0.,0.,1.), Eneut);
  
  std::cout << nProc << " processes" << std::endl;
  for (G4int i = 0; i < nProc; ++i) {
    G4HadronicProcess *thisProc = dynamic_cast<G4HadronicProcess*>((*processes)[i]);
    //std::cout << "i = " << i << std::endl;
    if (thisProc) {
      x = thisProc->GetElementCrossSection(dynamicNeutron, elm, rock);
      std::cout << "Cross section for " << thisProc->GetProcessName() << " : " << x 
                << std::endl;
    }
  }

  std::cout << "Done" << std::endl;
	
  return 0;
}
