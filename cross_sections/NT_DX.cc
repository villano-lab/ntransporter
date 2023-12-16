// estimate differential cross sections for given materials
//
// usage: ./NT_DX output_file_base_path path_to_ntransporter_base   
//                material1 material2 ... [-n ngroups1=100 ngroups2 ...]


#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4NistManager.hh"

#include "MinimalUserInitialization.hh"
#include "NTUtilities.hh"

#include "G4HadronElasticProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4HadProjectile.hh"
#include "G4HadFinalState.hh"
#include "G4Track.hh"

#include "G4SystemOfUnits.hh"

#include "G4ParticleHPManager.hh"


typedef std::vector<G4double> doubles;
typedef std::vector<G4int> ints;
typedef std::vector<G4String> strings;


void writeRows(G4String fileName, G4double **data, G4int G) {
    std::ofstream os(fileName);
    os << std::setprecision(17);
    for (int g1 = 0; g1 < G+2; ++g1) {
        os << g1 << ' ';
        for (int g2 = g1; g2 < G+2; ++g2) {
            if (data[g1][g2] > 0) {
                os << data[g1][g2] << ' ';
            } else {
                break;
            }
        }
        os << G4endl;
    }
    os.close();
}


void writeAllRows(G4String fileName, G4double **data, G4int G) {
    std::ofstream os(fileName);
    os << std::setprecision(17);
    for (int g1 = 0; g1 < G+2; ++g1) {
        for (int g2 = 0; g2 < G+2; ++g2) {
            os << data[g1][g2] << ' ';
        }
        os << G4endl;
    }
    os.close();
}


void writeColumns(G4String fileName, G4double **data, G4int G) {
    std::ofstream os(fileName);
    os << std::setprecision(17);
    for (int g2 = 0; g2 < G+2; ++g2) {
        os << g2;
        for (int g1 = g2; g1 > 0; --g1) {
            if (data[g1][g2] > 0) {
                os << ' ' << data[g1][g2];
            } else {
                break;
            }
        }
        os << G4endl;
    }
    os.close();
}


G4double reciprocalDistribution(G4double a, G4double b) {
    // draw from 1/x distribution between a and b
    return a*std::pow(b/a, G4UniformRand());
}

const G4String ERRORSTRING = "\n=========================================\n" 
                             "                 ERROR"  
                             "\n=========================================\n";

const G4String WARNINGSTRING = "\n=========================================\n" 
                             "                 WARNING"  
                             "\n=========================================\n";


int main(int argc, char **argv) {

    CLHEP::HepRandom::setTheSeed(time(NULL));

#ifdef G4MULTITHREADED  
    G4MTRunManager * runManager = new G4MTRunManager(); 
    G4cout << runManager->GetNumberOfThreads() << " threads" << G4endl;
#else
    G4RunManager * runManager = new G4RunManager();
    G4cout << "sequential mode" << G4endl;
#endif

    G4String output_file_base, NT_path_base;
    std::vector<G4String> material_names;
    ints Gs;

    // parse command line args
    if (argc < 4)  {
        G4cout << ERRORSTRING 
            << "Error in NT_DX: Not enough arguments." 
            << "\n    Usage: ./NT_DX output_file_base_path path_to_ntransporter_base "
            "material1 material2 ... [-n ngroups1=100 ngroups2 ...]" << G4endl;
        return 1;
    } else {
        try {

            output_file_base = argv[1];
            NT_path_base = argv[2];

            size_t i = 3;
            while (i < argc && argv[i][0] != '-' && argv[i][1] != 'n') {
                material_names.push_back(argv[i]);
                ++i;
            }
            for (size_t j = i+1; j < argc; ++j) {
                Gs.push_back(std::stoi(argv[j]));
            }
            if (Gs.size() == 0) {
                Gs.push_back(100);
            }
        } catch (std::invalid_argument) {
            G4cout << ERRORSTRING 
                << "Error in NT_DX: Invalid arguments. "
                << "ngroups must be numeric."
                << "\n    Usage: ./NT_DX output_file_base_path path_to_ntransporter_base "
            "material1 material2 ... [-n ngroups1=100 ngroups2 ...]" << G4endl;
            return 2;
        }
    }

    G4cout << "Pulling data from nTransporter at " << NT_path_base << G4endl;

    G4cout << "Gs: ";
    for (int G : Gs) {
        G4cout << G << " ";
    }
    G4cout << G4endl;

    G4cout << "Materials: ";
    for (G4String material_name : material_names) {
        G4cout << material_name << " ";
    }
    G4cout << G4endl;

    runManager->SetUserInitialization(new NeutronPhysicsList);
    runManager->SetUserInitialization(new MinimalDetector);
    runManager->SetUserInitialization(new MinimalInitialization);



    G4ThreeVector origin = G4ThreeVector(0.,0.,0.);


    // build materials
    G4Material *material;
    G4NistManager *nist = G4NistManager::Instance();
    for (G4String material_name : material_names) {
        material = nist->FindOrBuildMaterial(material_name);
        nist->BuildMaterialWithNewDensity(material_name+"_0", material_name, material->GetDensity(), 0., material->GetPressure());
    }
    G4cout << "Materials built" << G4endl;

    // bounds of fast energies
    const G4double Emin = 0.1*eV, Emax = 20.*MeV; 

    G4cout << "fast energy range " << Emin << " to " << Emax << " MeV" << G4endl;

    
    const G4ParticleDefinition *theNeutron = G4Neutron::Definition();

    // initialize physics list
    runManager->Initialize();
    G4cout << "Initialized" << G4endl;

    G4ProcessManager *theMan = theNeutron->GetProcessManager();
    G4ProcessVector *procList = theMan->GetProcessList();
    G4HadronElasticProcess *elasticProc = static_cast<G4HadronElasticProcess*>((*procList)[1]);
    //std::vector<G4HadronicInteraction*> interactionList = elasticProc->GetHadronicInteractionList();
    G4HadronicInteraction* elasticInteraction = elasticProc->GetHadronicInteractionList()[0];

   
    G4DynamicParticle *dynamicNeutron = new G4DynamicParticle(theNeutron, G4ThreeVector(0.,0.,1.), 0.);

    G4StepPoint *thePoint = new G4StepPoint();
    thePoint->SetPosition(origin);
    G4Step *theStep = new G4Step();
    theStep->SetPreStepPoint(thePoint);
    G4Track *neutronTrack = new G4Track(dynamicNeutron, 0., origin);
    neutronTrack->SetStep(theStep);
    G4HadProjectile *projectile = new G4HadProjectile(*neutronTrack);

    G4Nucleus *materialNucleus = new G4Nucleus;

    G4CrossSectionDataStore *elasticDataStore = elasticProc->GetCrossSectionDataStore();

    //G4Material *material;
    G4HadFinalState *neutronFS;
    G4double initialEnergy, finalEnergy;

    const G4int Gmax = *max_element(Gs.begin(), Gs.end());
    G4double Eg[Gmax+2], alpha;

    G4int **counts;
    G4double **probs, **sigs, **uncers;

    const G4Element *theElement;
    
    G4int gf, totalCounts, gmax, gdummy, numthrown;

    // variables for total cross section data
    G4String xs_tots_filename;
    G4double *xs_tots, Edummy, xdummy;

    //G4Element *theElement;

    // allocate space
    counts = (G4int**)malloc((Gmax+2)*sizeof(G4int*)); // MC counts
    probs = (G4double**)malloc((Gmax+2)*sizeof(G4double*)); // MC probabilities
    sigs = (G4double**)malloc((Gmax+2)*sizeof(G4double*)); // differential cross sections
    uncers = (G4double**)malloc((Gmax+2)*sizeof(G4double*)); // relative uncertainties
    xs_tots = (G4double*)malloc((Gmax+2)*sizeof(G4double)); // total scattering cross section

    for (int g = 0; g < Gmax + 2; ++g) {
        counts[g] = (G4int*)malloc((Gmax+2)*sizeof(G4int));
        probs[g] = (G4double*)malloc((Gmax+2)*sizeof(G4double));
        sigs[g] = (G4double*)malloc((Gmax+2)*sizeof(G4double));
        uncers[g] = (G4double*)malloc((Gmax+2)*sizeof(G4double));
    }

    // materials loop
    for (G4String material_name : material_names) {
        
        //material = nist->FindOrBuildMaterial(material_name);
        material = nist->FindOrBuildMaterial(material_name + "_0"); // zero-temp version
        
        thePoint->SetMaterial(material);

        G4cout << "Set material to " << material_name << G4endl;
        //G4cout << "Material temp = " << material->GetTemperature()/CLHEP::kelvin << " kelvin" << G4endl; // 293.15

        // loop over G
        for (G4int G : Gs) {

            // reset counts
            for (int g1 = 0; g1 < G+2; ++g1) {
                for (int g2 = 0; g2 < G+2; ++g2) {
                    counts[g1][g2] = 0;
                    probs[g1][g2] = 0;
                    sigs[g1][g2] = 0;
                    uncers[g1][g2] = 0;
                }
            }

            G4cout << "Reset arrays" << G4endl;

            // calculate Eg
            alpha = std::pow(Emin/Emax, 1./G);
            //G4cout << "alpha = " << alpha << G4endl;
            Eg[0] = Emax;

            //G4cout << "Calculated group edges" << G4endl;
            //G4cout << " Emax = " << Emax << G4endl;
            //G4cout << "Eg[0] = " << Eg[0] << G4endl;

            for (G4int g = 1; g < G+1; ++g) {
                Eg[g] = alpha*Eg[g-1];
                //G4cout << Eg[g] << G4endl;
            }
            Eg[G+1] = 1e-6*Eg[G];

            G4cout << "numthrown: ";

            


            // initial energy loop
            for (int g = 1; g < G+2; g++) {

                gmax = g;

                numthrown = 0;


                
                // throw until desired statistics reached
                do {

                    ++numthrown;

                    initialEnergy = reciprocalDistribution(Eg[g], Eg[g-1]);
                    //G4cout << "initial energy = " << initialEnergy << G4endl;
                    dynamicNeutron->SetKineticEnergy(initialEnergy);
                    
                    
                    
                    
                    //neutronTrack = new G4Track(dynamicNeutron, 0., origin);
                    neutronTrack->SetKineticEnergy(initialEnergy);
                    
                    
                    
                    
                    
                    //neutronTrack->SetStep(theStep);
                    projectile->Initialise(*neutronTrack);

                    //G4cout << "Set energy" << G4endl;
                    //G4cout << Emin << " < " << Eg[g] << " < Projectile energy = " << projectile->GetKineticEnergy() << " < " << Eg[g-1] << " < " << Emax << G4endl;

                    // select which nucleus gets hit
                    //elasticDataStore->ComputeCrossSection(dynamicNeutron, material); // need to compute cross sections in material before sampling Z/A 
                    //theElement = elasticDataStore->SampleZandA(dynamicNeutron, material, *materialNucleus);

                    //G4cout << "Sampled Z and A" << G4endl;
                    //G4cout << "Nucleus A = " << materialNucleus->GetA_asInt() << G4endl;


                    neutronFS = elasticInteraction->ApplyYourself(*projectile, *materialNucleus);

                    //G4cout << "Applied self" << G4endl;
                    finalEnergy = neutronFS->GetEnergyChange();

                    //G4cout << "Set final energy" << G4endl;
                    //G4cout << "Final energy = " << finalEnergy << G4endl;

                    /*if (finalEnergy > initialEnergy) {
                        G4cout << WARNINGSTRING << "upscattering\ng = " << g << ", gf = " << gf << G4endl;
                        //G4cout << "initialEnergy = " << initialEnergy << G4endl;
                        //G4cout << "finalEnergy = " << finalEnergy << G4endl;
                        //return 4;
                    }//*/

                    if (finalEnergy < Emin) { // thermal energy
                        gf = G+1;
                    } else {
                        
                        gf = static_cast<G4int>(std::ceil(std::log(finalEnergy/Eg[0])/std::log(alpha)));

                        if (gf == 0) {
                            G4cout << ERRORSTRING << "ZERO" << G4endl;
                            return 4;
                        }


                        // double check gf
                        if ((Eg[gf] > finalEnergy) || (finalEnergy >= Eg[gf-1])) {
                            G4cout << ERRORSTRING 
                                <<"Group calculation failed: finalEnergy = " 
                                + std::to_string(finalEnergy) << ", gf = " << gf << G4endl;
                            return 3;
                        }
                    }

                    if (gf > gmax) {
                        gmax = gf;
                    }

                    //G4cout << "Calculated gf" << G4endl;

                    // increment counter
                    ++counts[g][gf];

                    //G4cout << g << "/" << gf << " : counts: " << counts[g][g] << G4endl; 

                    //delete theElement;
                    //delete neutronTrack;

                } while (counts[g][gmax] < 1000 && counts[g][g] < 100000);
                // when desired error reached, break

                G4cout << ' ' << numthrown << ' ';

            } // end loop over initial group

            G4cout << G4endl;

            // calculate probs from counts
            for (int g1 = 1; g1 < G+2; ++g1) {
                totalCounts = 0;
                for (int g2 = 1; g2 < G+2; ++g2) {
                    totalCounts += counts[g1][g2];
                }
                for (int g2 = 0; g2 < G+2; ++g2) {
                    probs[g1][g2] = static_cast<G4double>(counts[g1][g2])/totalCounts;
                }
            }//*/


            // read in total cross section data
            xs_tots_filename = NT_path_base + "/cross_sections/data/V1/data_"
                                + material_name + "_"
                                + std::to_string(G)
                                + "_20_xs.dat";

            std::ifstream xs_stream(xs_tots_filename);

            G4cout << "Reading cross section data from " << xs_tots_filename << G4endl;

            for (int g = 0; g < G+2; ++g) {
                if (!(xs_stream >> gdummy >> Edummy >> xs_tots[g] >> xdummy)) {
                    G4cout << ERRORSTRING << "Error in NT_DX: reading variables for g = " << g << " failed in " << xs_tots_filename << G4endl;
                    return 101;
                }
                if (gdummy != g) {
                    G4cout << ERRORSTRING << "Error in NT_DX: line indices in " << xs_tots_filename << " don't match, g = " << g << G4endl;
                    return 102;
                }
                if (Eg[g] > 0 && g < G+1 && (std::abs(Edummy - Eg[g])/Eg[g] > 1e-12)) {
                    G4cout << ERRORSTRING << "Error in NT_DX: group boundary mismatch with " << xs_tots_filename << ", g = " << g << G4endl;
                    return 103;
                }
            }

            xs_stream.close();

            for (int g1 = 1; g1 < G+2; ++g1) {
                for (int g2 = 1; g2 < G+2; ++g2) {
                    sigs[g1][g2] = xs_tots[g1]*probs[g1][g2];
                }
            }


            // calculate stats uncertainties
            for (int g1 = 1; g1 < G+2; ++g1) {
                for (int g2 = 1; g2 < G+2; ++g2) {
                    if (counts[g1][g2] > 0) {
                        uncers[g1][g2] = 1./std::sqrt(counts[g1][g2]);
                    }
                }
            }

            G4cout << "Calculated counts" << G4endl;
            
            // write this material/group number to file
            G4String fileName = output_file_base + "_" 
                                + material_name + "_" 
                                + std::to_string(G);

            G4cout << "Writing differential cross section data to " << fileName << "_dx.dat" << G4endl;

            writeColumns(fileName + "_dx.dat", sigs, G);
            writeAllRows(fileName + "_uncers.dat", uncers, G);
            writeAllRows(fileName + "_alldx.dat", sigs, G);
            //writeAllProbs("test_comp.txt", counts, G);
            //writeAllProbs(output_file_base+"_"+material_name+"_"+std::to_string(G)+"_alldx.dat", counts, G);


        } // end loop over G

    } // end loop over materials
    

    delete runManager;

    G4cout << "Done" << G4endl;

    return 0;
}
