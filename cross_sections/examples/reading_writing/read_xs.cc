// hypothetically, script to grab cross section data from Geant
// first attempt at script - grab for neutrons incident on silicon

//#include "G4RunManagerFactory.hh"
//#include "G4UImanager.hh"

//#include "DetectorConstruction.hh"
//#include "G4PhysListFactory.hh"
//#include "G4VModularPhysicsList.hh"
//#include "PrimaryGeneratorAction.hh"
//#include "ActionInitialization.hh"
//#include "G4ParticleTable.hh"

//#include "G4QHadronInelasticDataSet.hh"
//#include "G4NeutronHPBGGNucleonInelasticXS.hh"
//#include "G4NeutronHPJENDLHEInelasticData.hh"
//#include "G4HadronNucleonXsc.hh"i
//#include "G4BGGNucleonInelasticXS.hh"

//#include "G4Neutron.hh"
#include "CDMSMaterialTable.hh"

#include "G4Neutron.hh"

#include "G4Element.hh"
#include "G4Material.hh"

#include <iostream>


int main(int argc, char **argv) {

  G4Neutron *theNeutron = G4Neutron::Definition();
  
  G4Material *rock = CDMSMaterialTable::GetMaterial("Norite");


  std::cout << "Did it" << std::endl;
	
  return 0;
}
