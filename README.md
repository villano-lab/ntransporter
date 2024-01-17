# `ntransporter v2.0`
Neutron transport code for background analysis.


## Overview

This code calculates the neutron flux in an infinite slab of material using the multigroup diffusion approximation with zeroth-order evaluated group constants (see below for description of this approximation).

## Directory Structure

This repository is divided into several directories by functionality. The overall repository structure is:

```
|- ntransporter/
|---- boundary_conditions/
|---- cross_sections/
|---- include/
|---- process_reader/
|---- sources/
|---- src/
```

More details on each top-level subdirectory can be found in individual `README`'s.


## Build Instructions

### Dependencies

The following packages must be installed/configured before building `ntransporter`:

- **`Geant4`** : required by `SuperSim`, also used for differential cross section calculations. Versions of `SuperSim` used by `ntransporter` (at time of writing, `<= V11.00.01`) use [`geant4-v10.6.3`](https://geant4.web.cern.ch/download/10.6.3.html). Installation instructions can be found [here](https://geant4-userdoc.web.cern.ch/UsersGuides/InstallationGuide/html/) for the newest version, but should still work for 10.6.3.


> [!NOTE]
> `Geant4` must be built in multithreaded mode in order to interface correctly with `ntransporter`. One should set the cmake flag `GEANT4_BUILD_MULTITHREADED=ON` during compilation of `Geant4`

- **`SuperSim`** : SuperCDMS's main simulations software, used for configuring CDMS-specific physics lists and calculating total cross sections. Data under the `supersim` directory is used for source calculations. The entire source code must be available to this program, and the environment should be configured as when building/compiling the software yourself (see [this Confluence page](https://confluence.slac.stanford.edu/display/CDMS/Running+Simulations+at+Specific+Sites)). The code is available on the [Gitlab](https://gitlab.com/supercdms/Simulations/supersim). 

    Linking to `SuperSim` is managed by the `SuperSim_include` and `SuperSim_exclude` files. An example of each is included in the top level of this repository (see [Top-level files](#top-level-files) section below). The `SuperSim_include` file gives a list of directories (under the `$CDMS_SUPERSIM` environment variable) to link and compile files in. The `SuperSim_exclude` file gives a list of file names (not including directory paths) to exclude from explicit linking. By default, the build system will parse all files matching `*.cc` by Linux [`find`](https://man7.org/linux/man-pages/man1/find.1.html) command unless the full filename matches "`/X`" for `X` any line given in the `SuperSim_exclude` file. 

    
    At the time of compilation, the environment must be configured just as one would in trying to build `SuperSim`. For example, if one wishes to run `ntransporter` on Cedar, follow the instructions on [this Confluence page](https://confluence.slac.stanford.edu/display/CDMS/Running+Simulations+on+ComputeCanada#RunningSimulationsonComputeCanada-CompilingSuperSimonCedar), and run the `setup_supersim.sh` script before compiling `ntransporter`. 


Depending on the system, `SuperSim` may require some additional packages to be correctly installed/configured. See also the [external package guide](https://confluence.slac.stanford.edu/display/CDMS/SuperSim+Reference+Manual#SuperSimReferenceManual-ExternalPackages). These include:


- **`ROOT`** : CERN's data processing language. Find installation instructions [here](https://root.cern/install/). Current versions of `ntransporter` have been tested on `ROOT` version `6.24/04`, but should work for any version `>= 6.02`. One should also note the C++ standard used to compile `ROOT` for setting the `CXX_STD` `cmake` variable

- `G4CMP` : Mike Kelsey's solid-state physics extension to `Geant4`. Available on [Github](https://github.com/kelseymh/G4CMP). Required by certain components of `SuperSim`, though not used directly by any `ntransporter` functionality. If one wishes to avoid this dependency, they may investigate modifying the `SuperSim_include` and `SuperSim_exclude`

- `CVODE` : LLNL's differential equation solver package. Can be downloaded as a component of the `SUNDIALS` package (find `v5.1.0` [here](https://github.com/LLNL/sundials/releases/tag/v5.1.0)), and find installation instructions [here](https://sundials.readthedocs.io/en/latest/Install_link.html). Note that `cmake` assumes the `CVODE` object files are built with the suffix `.so` unless the compiler matches "`AppleClang`", in which case it looks for `.dylib` files. The base `CMakeLists.txt` file must be modified directly if this is not the case.

- `uuid-dev` : utility package for generating 128-bit University Unique IDs (UUIDs) for random number generation. Can be downloaded with `sudo apt-get install uuid-dev`



### Environment variables

A few environment variables (denoted by `$` before the variable name) must be set before compiling `ntransporter` (in addition to the ones either set or required by the `SuperSim` configuration process, e.g., `$CDMS_INSTALL`, `$ROOTSYS`, or `$G4WORKDIR`). These are:

- `$G4CMPINSTALL` : path to the system's `G4CMP` installation (note: it does not appear to be necessary to source the `$G4CMPINSTALL/g4cmp_env.sh` script), e.g., `/home/ajbiffl3/G4CMP`. This one should also be set *before* sourcing the `SuperSim` `g4setup` script.

- `$CVODE_HOME` : path to installation directory of `CVODE`, e.g., (if following the installation guide) `/home/ajbiffl3/SUNDIALS/sundials-5.1.0/instdir`

- `$LIBUUID_OBJ` : (only used if `cmake` option `FORCE_LIBUUID_LINK` is set) path to `uuid-dev` object file, e.g., `/usr/lib/libuuid.so`

### `cmake` variables

`cmake` requires some variables giving paths to important files or other programs be set correctly in order to compile:

- `SUPERSIM_INCLUDE` : the full path of the `SuperSim_include` file. If one wishes to use `SuperSim_include.txt` included at the top level of this repository, the default value ("`SuperSim_include.txt`") will be fine.

- `SUPERSIM_EXCLUDE` : the full path of the `SuperSim_exclude` file. Once again, the default value ("`SuperSim_exclude.txt`") will use the `SuperSim_exclude` file included in this repository.

- `CXX_STD` : integer (default `14`) specifying the C++ standard to compile the program with. It is recommended to compile against the same version used by `ROOT`. The default value of `14` is typically used to compile `SuperSim`. Conda installations of `ROOT` use `17`.

- `FORCE_LIBUUID_LINK` boolean option (default `OFF`) to force linking of the `libuuid` object file specified by the `$LIBUUID_OBJ` environment variable (if set, `$LIBUUID_OBJ` must be specified in the current environment to point to the `libuuid.so` object file on the current system, e.g., `/usr/lib/libuuid.so`). 

- `CVODE_LIB_FOLDER` : name of folder under environment variable `$CVODE_HOME` containing object files for `CVODE` (probably either `lib`, the default value, or `lib64`). 


Some other optional cmake options can be used to personalize the build:

- `REGEN_CDMSVERSION` : boolean option (default `OFF`) to force re-generation of the `CDMSVersion.hh` header file. This option will always reset to `OFF` after running, so must be explicitly set every time one wants to regenerate the version header file. This should be done any time the `SuperSim` version is changed.

- `IGNORE_WARNINGS` : boolean option (default `OFF`) to ignore compiler warnings (passes `-w` flag to compiler if set)

- `SUBDIRS_VERBOSE` : boolean option (default `OFF`) to print additional information about linking executables in first-level subdirectories

- `EXCLUDE_VERBOSE` : boolean option (default `ON`) to print files excluded from compilation/linking by the `SuperSim_exclude` file

- `BUILD_PROCINFO` : boolean option (default `OFF`) to build the `PROCINFO` executable from the `process_reader` subdirectory (see that `README` for more info)

These options can be set when calling `cmake` or by modifying the `CMakeCache.txt` file in the build directory directly. All but `REGEN_CDMSVERSION` will be saved in `CMakeCache.txt` for subsequent calls of `cmake`.

### Building the program

The `ntransporter` program uses `cmake` utility to compile and build the executables in the program. Note that several steps in the cmake script require using the linux shell.

To build the program, make a build directory. For these examples, it's assumed the build directory is located inside the top-level `ntransporter` directory:

```
$ cd ntransporter
$ mkdir build
$ cd build
```


Once the environment variables are correctly set up, the dependencies configured, and the necessary `cmake` variables set, `cmake` can be run, e.g.,

```
cmake -DFORCE_LIBUUID_LINK=ON -DCVODE_LIB_FOLDER="lib64" -DCXX_STD=17 ../
```



If this is successful, `make` all targets:

```
make all
```

This will build all executables in `ntransporter` into the build directory.

Alternatively, one can build each executable individually, e.g.,

```
make NT_Src
```

## Executables Quick Guide

An overview of the executables built by `ntransporter` is given here. In general, the executables built by `ntransporter` begin with "`NT_`".

- `NT_XS` : calculate group cross sections for one material/group number pair

    Usage: `./NT_XS output_file_base_path material [ngroups=100] [points_per_group=10]`

- `NT_XX` : calculate group cross sections for multiple material/group number pairs (note the results change depending on the order of pairs - use is not recommended)

    Usage: `./NT_XX output_file_base_path material1 material2 ... [-n ngroups1=100 ngroups2 ...]`

- `NT_DX` : calculate differential cross sections for multiple material/group number pairs

    Usage: `./NT_DX output_file_base_path path_to_ntransporter_base material1 material2 ... [-n ngroups1=100 ngroups2 ...]`

- `NT_Src` : calculate group sources 

    Usage: `./NT_Src output_file_base_path material path_to_supersim [ngroups=100]`

- `NT_BC` : calculate boundary conditions (infinite slab flux)

    Usage: `./NT_BC output_file_base_path path_to_ntransporter_base material [ngroups=100]`

- `PROCINFO` : print hadronic process info for the neutron in the "Shielding" physics list in `SuperSim`

    Usage: `./PROCINFO`

### Note About Outputs

The repository currently includes many output files from the executables in `data/V1` and `data/V2` directories under the subdirectory of the corresponding executable, e.g., `ntransporter/sources/data/V1` contains files generated with the `NT_Src` executable.

## Top-level files

The top-level `ntransporter/` directory contains the following files:

- `.gitignore` : standard `.gitignore` file

- `CMakeLists.txt` : top-level CMakeLists file for the project. Responsible for compiling all executables in the repository

- `LICENSE` :  Standard MIT license

- `SuperSim_include.txt` : `SuperSim_include` file, list of directories within SuperSim needed for project

- `SuperSim_exclude.txt` : `SuperSim_exclude` file, list of filenames to exclude from linking if found in any of the `SuperSim_include` directories

- `SuperSim_Main.hh` : overloaded `SuperSim_Main` header file from `SuperSim` - needed to redefine members of the `SuperSim_Main` class as public so they can be accessed by `ntransporter`

## Version 2 Physics - Zero-Temperature Multigroup Diffusion Approximation with Zeroth-Order Evaluated Group Constants


The full neutron transport equation is: 

$\frac{1}{v}\frac{\partial \psi}{\partial t} + \boldsymbol{\hat{\Omega}} \cdot \nabla \psi + \Sigma_t \psi = s + \int_{4\pi} d\boldsymbol{\Omega'} \int_0^\infty dE' \Sigma_s(E'\rightarrow E, \boldsymbol{\hat{\Omega}'}\rightarrow\boldsymbol{\hat{\Omega}}) \psi (E',\boldsymbol{\hat{\Omega}'})$

where

$\boldsymbol{r}=$ position

$E =$ energy

$\boldsymbol{\hat{\Omega}} =$ direction of neutron travel

$t =$ time

$v=$ neutron velocity

$\psi = \psi(\boldsymbol{r}, E, \boldsymbol{\hat{\Omega}},t) =$ angular neutron flux (neutrons per unit time per unit area per unit energy per unit solid angle) 

$\Sigma_t = \Sigma_t(\boldsymbol{r}, E) =$ total neutron cross section (scattering plus absorption), units of inverse length

$s=s(\boldsymbol{r},E,\boldsymbol{\hat{\Omega}})=$ neutron source (neutrons per unit time per unit volume per unit energy per unit solid angle) 

$\Sigma_s(E'\rightarrow E,\boldsymbol{\hat{\Omega}'}\rightarrow\boldsymbol{\hat{\Omega}})$ $=\Sigma_s(E,E',\boldsymbol{\hat{\Omega}'}\cdot\boldsymbol{\hat{\Omega}})=$ scattering cross section from primed to unprimed state (differential over final direction and final energy)

To move to the **neutron diffusion equation**, we first define the scalar flux $\phi$:

$\phi(E) \equiv \int_{4\pi} \psi(E,\boldsymbol{\hat{\Omega}}) d\boldsymbol{\Omega} $

and the neutron current density $\boldsymbol{J}$:

$\boldsymbol{J}(E) \equiv \int_{4\pi} \boldsymbol{\hat{\Omega}} \psi(E,\boldsymbol{\hat{\Omega}}) d\boldsymbol{\Omega}$

(both with units neutrons per unit time per unit area per unit energy) and make two physical assumptions. 


First, assume the angular flux $\psi$ can be modelled sufficiently with linear anisotropy:

$\psi(\boldsymbol{\hat{\Omega}}) \approx \frac{1}{4\pi} \phi + \frac{3}{4\pi} \boldsymbol{J}\cdot \boldsymbol{\hat{\Omega}}$

Second, assume the neutron current density $\boldsymbol{J}$ is proportional to the gradient of the scalar flux ("Fick's Law"):

$\boldsymbol{J} = [3(\Sigma_t-\bar{\mu}_0\Sigma_s)]^{-1} \nabla \phi \equiv -D(E)\nabla\phi(E)$


Additional physical assumptions: isotropic scattering and isotropic sources:


$\Sigma_s(E'\rightarrow E, \boldsymbol{\hat{\Omega}'}\rightarrow\boldsymbol{\hat{\Omega}}) = \frac{1}{4\pi}\Sigma_s(E'\rightarrow E)$

$s(\boldsymbol{\hat{\Omega}},E) = \frac{1}{4\pi}S(E)$

where $S(E)=S(\boldsymbol{r}, E, t)$ is the scalar source.


Now if we integrate over solid angle $\boldsymbol{\Omega}$, the transport equation becomes the diffusion equation:

$\frac{1}{v}\frac{\partial \phi}{\partial t} - \nabla  D(E) \nabla \phi(E) + \Sigma_t(E) \phi(E) = S(E) + \int_0^\infty dE' \Sigma_s(E'\rightarrow E) \phi (E')$





### The Multigroup Approximation

We now define a series of $G$ energy intervals called "groups" delineated by the energies $E_g$, where $0\leq g\leq G$, such that $E_{g+1} < E_g \forall g$. Group $g$ refers to the range of energies $E_{g} \leq E \leq E_{g-1}$. $E_0$ corresponds to the maximum energy attainable by neutrons in the system (for radiogenic neutrons, $E_0\sim$ 20 MeV), and $E_G$ corresponds to the minimum energy of interest, typically thermal energies. Introduce the notation for integrating over group $g$: 

$\int_g dE \cdots \equiv \int_{E_{g}}^{E_{g-1}} dE \cdots$


Integrating the diffusion equation over group $g$ then yields the **multigroup diffusion equation**:


$\frac{1}{v_g}\frac{\partial \phi_g}{\partial t} - \nabla \cdot D_g \nabla \phi_g + \Sigma_{tg} \phi_g = S_g + \sum\limits_{g'=1}^{\infty} \Sigma_{sg'g} \phi_{g'}$


where we've defined the group flux $\phi_g$:
  
$\phi_g \equiv \int_g dE \phi(E)$

and the group source $S_g$:

$S_g \equiv \int_g dE S(E)$

and where we've defined the following "group constants:"

- The total cross section:

$\hspace{5ex}\Sigma_{tg} \equiv \frac{1}{\phi_g} \int_g dE \Sigma_t(E) \phi(E)$

- The diffusion coefficient on the $j$ component of the flux gradient:

$\hspace{5ex}D_{gj} \equiv \frac{\int_g dE D(E) \nabla_j \phi(E) }{\int_g dE \nabla_j \phi(E)}$


- The inverse neutron speed:

$\hspace{5ex}\frac{1}{v_g} \equiv \frac{1}{\phi_g}\int_g dE \frac{1}{v}\phi(E)$

- The scattering cross section from group $g'$ to group $g$:

$\hspace{5ex}\Sigma_{sg'g} \equiv \frac{1}{\phi_{g'}} \int_g dE \int_{g'} dE' \Sigma_s(E'\rightarrow E) \phi(E')$



Also note the spatial derivative term should be expanded as the following:

$\nabla \cdot D_g \nabla \phi_g = \sum\limits_j \nabla_j D_{gj} \nabla_j \phi_g$


### Note on Group Structure

We consider a group structure with fast groups between 0.1 eV and 20 MeV and a single thermal group below 0.1 eV. The actual numerical value of the lower bound of the thermal group is set to 1e-7 eV


Given these specifications, a grouping (set of group boundaries) is specified by the number of fast groups. In the code this is what `G` refers to, while $G$ in the derivations here usually refers to the *total* number of groups, `G+1`.

With these in place, the bounds of the fast groups are $E_G$ = 0.1 eV and $E_0$ = 20 MeV, and any particular group boundary can be calculated as $E_g=E_0\beta^g$, where $\beta=(E_G/E_0)^{1/G}$

### The Slowing-Down Equation


From the neutron diffusion equation, we consider the steady-state, infinite medium case, and so can drop the time and spatial derivatives, resulting in the "slowing down equation:"

$\Sigma_t(E) \phi(E) = S(E) + \int_0^{\infty} dE' \Sigma_s(E'\rightarrow E) \phi (E')$

where now everything is only a function of energy.


Once again integrating over group $g$, we arrive at the **multigroup slowing down equation**:


$\Sigma_{tg} \phi_g = S_g + \sum\limits_{g'=1}^{\infty} \Sigma_{sg'g} \phi_{g'}$




### Zero-Temperature Approximation

In the zero-temperature approximation, it is impossible for a neutron to increase in energy (decrease in group number) after a collision, so the multigroup slowing down equation becomes:


$(\Sigma_{tg} -\Sigma_{sgg})\phi_g = S_g + \sum\limits_{g'=1}^{g-1} \Sigma_{sg'g} \phi_{g'}$



### Evaluating Group Constants

The evaluation of most group constants requires knowledge of the flux $\phi(E)$. We use a zeroth-order approximation for fast groups of $\phi(E)\propto 1/E$, the derivation of which relies on the following approximations:

- Assume s-wave scattering: $\Sigma_s(E'\rightarrow E)=\frac{\Sigma_s(E')}{(1-\alpha)E'}$ when $E\lt E'\lt E/\alpha$

- Ignore absorption and inelastic scattering

- Consider energies well below source energies

These above assumptions result in an asymptotic solution $\phi(E)\propto \frac{1}{\Sigma_s(E) E}$ (note also this is the form used to derive the approximations for the self-scattering fraction $q$ above). We then neglect the variation in $\Sigma_s(E)$ over the group: $\phi(E)\propto \frac{1}{E}$.

For the thermal group, we assume a Maxwell-Boltzmann distribution at room temperature (0.0257 eV).

More details on evaluating the differential cross sections $\Sigma_{sg'g}$ can be found in the `cross_sections` documentation.
