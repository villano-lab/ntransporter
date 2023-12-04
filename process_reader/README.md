# `ntransporter/process_reader/`

Directory for test code/code to get more info about hadronic processes and other information in `SuperSim` instance


## Executables

This directory generates the executable `PROCINFO` if the cmake option `BUILD_PROCINFO` is set.

### `PROCINFO`

`PROCINFO` prints process information and descriptions from `SuperSim`. It is run without any arguments and prints its output to shell. This output should be the same every time (notwithstanding very major changes to the physics lists used by SuperSim). 


Going through the full output, first the program prints some status reports as it boots up a SuperSim instance: 

```
Starting SuperSim_Main
Configuring physics list
Configuring environment geometry
Initializing
```

The program then calls up the `G4Neutron` and fetches its process manager, which is responsible for keeping track of and choosing what processes influence a particular neutron at every simulation step.

```
Fetching neutron singleton
Fetching neutron process manager
```


From the process manager, the program accesses the list of all processes implemented to affect the neutron. It prints the number of processes:

```
9 processes
```

Then it goes through each process, printing the index of the process in the list of processes, the name, and type, and then calls the `DumpInfo` and `ProcessDescription` methods of the process. This results in the following:

```
Process 0 :  name = CoupledTransportation, type = Transportation
This process has not yet been described
```
```
Process 1 :  name = Decay, type = Decay
Decay: Decay of particles. 
kinematics of daughters are dertermined by DecayChannels  or by PreAssignedDecayProducts
```
```
Process 2 :  name = hadElastic, type = Hadronic
G4HadronElasticProcess handles the elastic scattering of 
hadrons by invoking the following hadronic model(s) and 
hadronic cross section(s).
```
```
Process 3 :  name = neutronInelastic, type = Hadronic
G4NeutronInelasticProcess handles the inelastic scattering of
neutrons from nuclei by invoking one or more hadronic models and
one or more hadronic cross section sets.
```
```
Process 4 :  name = nCapture, type = Hadronic
G4HadronCaptureProcess handles the capture of neutrons by nuclei
following by gamma/electron de-excitation cascade. One or more
hadronic models and hadronic cross section sets may be invoked.
```
```
Process 5 :  name = nFission, type = Hadronic
G4HadronFissionProcess handles neutron-induced fission of nuclei
by invoking one or more hadronic models and one or more hadronic
cross sections.
```
```
Process 6 :  name = nKiller, type = General
This process has not yet been described
```
```
Process 7 :  name = UserSpecialCut, type = General
This process has not yet been described
```
```
Process 8 :  name = Scorers, type = ------
This process has not yet been described
```

The first process is a `Transportation` process, which all particles have in `Geant4`. This covers unimpeded transportation through materials and across material boundaries. Process 2 is the natural decay process, again present in almost all particles. For neutrons the decay process is negligible, as the half-life of the neutron is on the order of 15 minutes (compare to ~ns time-of-flight between elastic scatters). The `nKiller`, `UserSpecialCut`, and `Scorers` processes are special user-controlled processes to kill or track neutrons and are not relevant to us. The four remaining processes represent real physical interactions that neutrons undergo during transport, and are the ones of interest to us:

- `hadElastic` : elastic scattering
- `neutronInelastic` : inelastic scattering
- `nCapture` : thermal capture, $(n,\gamma)$
- `nFission` : neutron-induced fission


### Full Output





The full output is:


```
Starting SuperSim_Main
Configuring physics list
Configuring environment geometry
Initializing
Fetching neutron singleton
Fetching neutron process manager
9 processes
Process 0 :  name = CoupledTransportation, type = Transportation
This process has not yet been described
Process 1 :  name = Decay, type = Decay
Decay: Decay of particles. 
kinematics of daughters are dertermined by DecayChannels  or by PreAssignedDecayProducts
Process 2 :  name = hadElastic, type = Hadronic
G4HadronElasticProcess handles the elastic scattering of 
hadrons by invoking the following hadronic model(s) and 
hadronic cross section(s).
Process 3 :  name = neutronInelastic, type = Hadronic
G4NeutronInelasticProcess handles the inelastic scattering of
neutrons from nuclei by invoking one or more hadronic models and
one or more hadronic cross section sets.
Process 4 :  name = nCapture, type = Hadronic
G4HadronCaptureProcess handles the capture of neutrons by nuclei
following by gamma/electron de-excitation cascade. One or more
hadronic models and hadronic cross section sets may be invoked.
Process 5 :  name = nFission, type = Hadronic
G4HadronFissionProcess handles neutron-induced fission of nuclei
by invoking one or more hadronic models and one or more hadronic
cross sections.
Process 6 :  name = nKiller, type = General
This process has not yet been described
Process 7 :  name = UserSpecialCut, type = General
This process has not yet been described
Process 8 :  name = Scorers, type = ------
This process has not yet been described
```