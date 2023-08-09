// hypothetically, script to grab cross section data from Geant
// first attempt at script - grab for neutrons incident on silicon

#include <iostream>
#include <vector>

#include "G4ThreeVector"
#include "G4SystemOfUnits.hh"

#include "CDMSMaterialTable.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4ElementVector"

#include "G4Neutron.hh"
#include "Shielding.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"
#include "G4HadronicProcess.hh"


int main(int argc, char **argv) {

  G4int verbose = 0;

  G4Neutron *theNeutron = G4Neutron::Definition();
  const CDMSMaterialTable *theTable = CDMSMaterialTable::GetInstance();
  
  G4Material *rock = theTable->GetMaterial("Norite");
  G4int nElm = rock->GetNumberOfElements();
  G4ElementVector *elmVector = rock->GetElementVector();
  G4double *fracVector = rock->GetFractionVector();

  //G4NeutronHPElasticData *elastic;
  //G4NeutronHPCaptureData *capture;
  //G4NeutronHPInelasticData *inelastic;

  Shielding *shield = new Shielding(verbose); // create Shielding physics list object  
  shield->ConstructProcess(); // construct processes

  //typedef std::vector<G4PhysicsConstructor*> G4PhysConstVectorData;
  //G4PhysConstVectorData* physicsVector = 
  //                  shield->GetSubInstanceManager().offset[0].physicsVector;


  // get ProcessManager for the neutron
  G4ProcessManager *theMan = theNeutron->GetProcessManager();

  // get element from material
  G4int nProc = theMan->GetProcessListLength();
  G4ProcessVector *processes = theMan->GetProcessList();
  G4Element *elm = elmVector[6]; // silicon
  
  // variables to store cross sections, neutron energy
  G4double x, Eneut = 2.*keV;


  std::cout << "Cross sections for " << elm->GetName() << " with " << Eneut 
            << " MeV neutron" << std:endl;

  G4DynamicParticle *dynamicNeutron = new G4DynamicParticle(theNeutron, 
                                              G4ThreeVector(0.,0.,1.), Eneut);

  for (size_t i; i < nProc; ++i) {
    G4HadronicProcess *thisProc = dynamic_cast<G4HadronicProcess> processes[i];
    if (thisProc) {
      x = thisProc->GetElementCrossSection(dynamicNeutron, elm, rock);
      std::cout << "Cross section for " << thisProc->GetName() << " : " << x 
                << std::endl;
    }
  }

  std::cout << "Done" << std::endl;
	
  return 0;
}
