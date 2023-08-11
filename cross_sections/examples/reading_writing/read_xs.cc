// script to grab cross section data from Geant
// first attempt at script - print process cross sections for specific energy 

#include <iostream>
#include <vector>

#include "SuperSim_Main.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

#include "CDMSMaterialTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ElementVector.hh"

#include "G4Neutron.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4HadronicProcess.hh"


int main(int argc, char **argv) {
  
  // manager object to configure everything
  SuperSim_Main *sMain = new SuperSim_Main();

  // set verbosity
  sMain->SetVerboseLevel(0);
  sMain->runManager.SetVerboseLevel(0);

  sMain->Configure("full"); // full Shielding physics list, etc
  
  // configure geometry and detector (required for initialization)
  sMain->UImanager->ApplyCommand("/CDMS/lab NoLab");
  sMain->UImanager->ApplyCommand("/CDMS/detector ZIP");
  
  // initialize: among other things, build physics processes and attach to 
  // particle process managers
  sMain->runManager.Initialize();

  // neutron singleton
  G4Neutron *theNeutron = G4Neutron::Definition();

  // pull material table from SuperSim
  const CDMSMaterialTable *theTable = CDMSMaterialTable::GetInstance();
  
  // pull material data
  G4Material *rock = theTable->GetMaterial("Norite");
  G4int nElm = rock->GetNumberOfElements();
  const G4ElementVector *elmVector = rock->GetElementVector();
  const G4double *fracVector = rock->GetFractionVector();
  const G4Element *elm = (*elmVector)[6]; // silicon

  
  // get ProcessManager for the neutron
  G4ProcessManager *theMan = theNeutron->GetProcessManager();
  
  // print neutron particle instance ID (should be >0)
  if (0) std::cout << "n ID = " << theNeutron->GetInstanceID() << std::endl;
  
  // get vector of neutron processes
  G4int nProc = theMan->GetProcessListLength();
  G4ProcessVector *processes = theMan->GetProcessList();
  
  // variables to store cross sections, neutron energy
  G4double x, Eneut = 2.*keV;


  std::cout << "Cross sections for " << elm->GetName() << " with " << Eneut 
            << " MeV neutron" << std::endl;
  
  // dynamic particle: set energy, momentum
  G4DynamicParticle *dynamicNeutron = new G4DynamicParticle(theNeutron, 
                                              G4ThreeVector(0.,0.,1.), Eneut);
  
  std::cout << nProc << " processes" << std::endl;
  for (G4int i = 0; i < nProc; ++i) {
    // if process is a hadronic process, print cross section
    G4HadronicProcess *thisProc = dynamic_cast<G4HadronicProcess*>((*processes)[i]);
    if (thisProc) {
      // calculate and print
      x = thisProc->GetElementCrossSection(dynamicNeutron, elm, rock);
      std::cout << "Cross section for " << thisProc->GetProcessName() 
                << " : " << x << std::endl;
    }
  }

  std::cout << "Done" << std::endl;
	
  return 0;
}
