# `ntransporter v1.0`
Neutron transport code for background analysis.


## Overview

This code calculates the neutron flux in an infinite slab of material using the multigroup diffusion approximation with zeroth-order evaluated group constants (see below for derivation/description of this approximation).

## Directory Structure

This repository is divided into several directories by function. The overall repository structure is:

```
|- ntransporter/
|---- boundary_conditions/
|------- data/
|---- cross_sections/
|------- data/
|------- examples/
|---- process_reader/
|---- include/
|------- 
|---- sources/
|------- 
|---- src/
|------- 
```

More details on each top-level subdirectory can be found in individual `README`'s.


## Build Instructions

The `ntransporter` program uses `cmake` utility to compile and build the executables in the program. Note that several steps in the build process require using the linux shell.

To build the program, make a build directory. For these examples, it's assumed the build directory is located inside the top-level `ntransporter` directory:

```
$ cd ntransporter
$ mkdir build
$ cd build
```

At the time of compilation, the environment must be configured just as one would in trying to build `SuperSim` itself. For example, if one wishes to run `ntransporter` on Cedar, follow the instructions on [this Confluence page](https://confluence.slac.stanford.edu/display/CDMS/Running+Simulations+on+ComputeCanada#RunningSimulationsonComputeCanada-CompilingSuperSimonCedar), and run the `setup_supersim.sh` script before compiling `ntransporter`. This should set all required environment variables. One can check if this is successful by checking the value of `CDMS_SUPERSIM`:

```
$ echo $CDMS_SUPERSIM
```

which should display the full path of the top-level `supersim` directory of `SuperSim`, e.g., `/home/ajbiffl3/supersim`.

In addition, `ntransporter` requires a few variables giving paths to important files or other programs. These are:

- `SUPERSIM_INCLUDE` : the full path of the `SuperSim_include.txt` file. If one wishes to use the `SuperSim_include.txt` file included at the top level of this directory, the default value of "`SuperSim_include.txt`" will be fine.



Some other optional cmake options can be used to personalize the build:

- `REGEN_CDMSVERSION` : boolean option (default `OFF`) to force re-generation of the `CDMSVersion.hh` header file. This option will always reset to `OFF` after running, so must be explicitly set every time one wants to regenerate the version header file. This should be done any time the `SuperSim` version is changed.

- `SUBDIRS_VERBOSE` : boolean option (default `OFF`) to print additional information about linking executables in top-level subdirectories

- `BUILD_PROCINFO` : boolean option (default `OFF`) to build the `PROCINFO` executable from the `process_reader` subdirectory (see that `README` for more info)

These options can be set when calling `cmake`, e.g., from the build directory:

```
cmake -DREGEN_CDMSVERSION=ON -DSUPERSIM_INCLUDE="../SuperSim_include.txt" ../
```

And then `make` all targets:

```
make all
```


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

- `NT_Src` : calculate group sources 

Usage: `./NT_Src output_file_base_path material path_to_supersim [ngroups=100]`

- `NT_BC` : calculate boundary conditions (infinite slab flux)

Usage: `./NT_BC output_file_base_path path_to_ntransporter_base material [ngroups=100]`

- `PROCINFO` : print hadronic process info for the neutron in the "Shielding" physics list in `SuperSim` and calculate total cross sections for neutron of different energies

Usage: `./PROCINFO`

### Note About Outputs

The repository currently includes many output files from the executables in `data/V1` directories under the subdirectory of the corresponding executable, e.g., `ntransporter/sources/data/V1` contains files generated with the `NT_Src` executable.

## Top-level files

The top-level `ntransporter/` directory contains the following files:

- `.gitignore` : standard `.gitignore` file

- `CMakeLists.txt` : top-level CMakeLists file for the project. Responsible for compiling all executables in the repository

- `LICENSE` :  Standard MIT license

- `NTUtilities.hh` : Header file containing functions and constants used in multiple of the top-level directories, ex. thermal energy, the kernel of the Maxwell-Boltzmann energy distribution, etc.

- `SuperSim_include.txt` : list of directories within SuperCDMS [`SuperSim`](https://gitlab.com/supercdms/Simulations/supersim) needed for project

- `SuperSim_Main.hh` : overloaded `SuperSim_Main` header file from `SuperSim` - needed to redefine members of the `SuperSim_Main` class as public so they can be accessed by `ntransporter`

## Version 1 Physics - Direct-Coupled Multigroup Diffusion Approximation with Zeroth-Order Evaluated Group Constants


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

$s=s(\boldsymbol{r}, E, \boldsymbol{\hat{\Omega}}) =$ neutron source (neutrons per unit time per unit volume per unit energy per unit solid angle) 

$\Sigma_s(E'\rightarrow E, \boldsymbol{\hat{\Omega}'}\rightarrow\boldsymbol{\hat{\Omega}})=\Sigma_s(E,E',\boldsymbol{\hat{\Omega}'}\cdot\boldsymbol{\hat{\Omega}}) = $ scattering cross section from primed to unprimed state (differential over final direction and final energy)

To move to the **neutron diffusion equation**, we first define the scalar flux $\phi$:

$\phi(E) \equiv \int_{4\pi} \psi(E,\boldsymbol{\hat{\Omega}}) d\boldsymbol{\Omega} $

and the neutron current density $\boldsymbol{J}$:

$\boldsymbol{J}(E) \equiv \int_{4\pi} \boldsymbol{\hat{\Omega}} \psi(E,\boldsymbol{\hat{\Omega}}) d\boldsymbol{\Omega}$

(both with units neutrons per unit time per unit area per unit energy) and make two physical assumptions. 


First, assume the angular flux $\psi$ can be modelled sufficiently with linear anisotropy:

$\psi \approx \frac{1}{4\pi} \phi + \frac{3}{4\pi} \boldsymbol{J}\cdot \boldsymbol{\hat{\Omega}}$

Second, assume the neutron current density $\boldsymbol{J}$ is proportional to the gradient of the scalar flux ("Fick's Law"):

$\boldsymbol{J} = [3(\Sigma_t-\bar{\mu}_0\Sigma_s)]^{-1} \nabla \phi \equiv -D(E)\nabla\phi(E)$


Additional physical assumptions: isotropic scattering and isotropic sources:


$\Sigma_s(E'\rightarrow E, \boldsymbol{\hat{\Omega}'}\rightarrow\boldsymbol{\hat{\Omega}}) = \frac{1}{4\pi}\Sigma_s(E'\rightarrow E)$

$s(\boldsymbol{\hat{\Omega}},E) = \frac{1}{4\pi}S(E)$

where $S(E)=S(\boldsymbol{r}, E, t)$ is the scalar source.


Now if we integrate over solid angle $\boldsymbol{\Omega}$, the transport equation becomes the diffusion equation:

$\frac{1}{v}\frac{\partial \phi}{\partial t} - \nabla  D(E) \nabla \phi(E) + \Sigma_t(E) \phi(E) = S(E) + \int_0^\infty dE' \Sigma_s(E'\rightarrow E) \phi (E')$





### The Multigroup Approximation

We now define a series of $G$ energy intervals called ``groups" delineated by the energies $E_g$, where $0\leq g\leq G$, such that $E_{g+1} < E_g \forall g$. Group $g$ refers to the range of energies $E_{g} \leq E \leq E_{g-1}$. $E_0$ corresponds to the maximum energy attainable by neutrons in the system (for radiogenic neutrons, $E_0\sim 20\,\text{MeV}$), and $E_G$ corresponds to the minimum energy of interest, typically thermal energies. Introduce the notation for integrating over group $g$: 

$\int_g dE \cdots \equiv \int_{E_{g}}^{E_{g-1}} dE \cdots$


Integrating the diffusion equation over group $g$ then yields the **multigroup diffusion equation**:


$\frac{1}{v_g}\frac{\partial \phi_g}{\partial t} - \nabla \cdot D_g \nabla \phi_g + \Sigma_{tg} \phi_g = S_g + \sum_{g'=1}^\infty \Sigma_{sg'g} \phi_{g'}$


where we've defined the group flux $\phi_g$:
  
$\phi_g \equiv \int_g dE \phi(E)$

and the group source $S_g$:

$S_g \equiv \int_g dE S(E)$

and where we've defined the following "group constants:"

- The total cross section:

$\hspace{5ex}\Sigma_{tg} \equiv \frac{1}{\phi_g} \int_g dE \Sigma_t(E) \phi(E)$

- The diffusion coefficient on the $j$ component of the flux gradient:

$\hspace{5ex}D_{gj} \equiv \frac{\int_g dE D(E) \nabla_j \phi(E) }{\int_g dE \nabla_j \phi(E)}$


- The neutron speed:

$\hspace{5ex}\frac{1}{v_g} \equiv \frac{1}{\phi_g}\int_g dE \frac{1}{v}\phi(E)$

- The scattering cross section from group $g'$ to group $g$:

$\hspace{5ex}\Sigma_{sg'g} \equiv \frac{1}{\phi_{g'}} \int_g dE \int_{g'} dE' \Sigma_s(E'\rightarrow E) \phi(E')$



Also note the spatial derivative term should be expanded as the following:

$\nabla \cdot D_g \nabla \phi_g = \sum_j \nabla_j D_{gj} \nabla_j \phi_g$


### Note on Group Structure

We consider a group structure with fast groups between 0.1 eV and 20 MeV and a single thermal group below 0.1 eV. The actual numerical value of the lower bound of the thermal group is set as machine epsilon times 0.1 eV (note that going this low is unnecessary and causes a few minor numerical problems, and later version change this to 1e-7 eV).


Given these specifications, a grouping (set of group boundaries) is specified by the number of fast groups. In the code this is what `G` in the code refers to, while $G$ in the derivations here refers to the *total* number of groups, `G+1`.

With these in place, the bounds of the fast groups are $E_G$ = 0.1 eV and $E_0$ = 20 MeV, and any particular group boundary can be calculated as $E_g=E_0\beta^g$, where $\beta=(E_G/E_0)^{1/G}$

### The Slowing-Down Equation


From the neutron diffusion equation, we consider the steady-state, infinite medium case, and so can drop the time and spatial derivatives, resulting in the "slowing down equation:"

$\Sigma_t(E) \phi(E) = S(E) + \int_0^\infty dE' \Sigma_s(E'\rightarrow E) \phi (E')$

where now everything is only a function of energy.


Once again integrating over group $g$, we arrive at the **multigroup slowing down equation**:


$\Sigma_{tg} \phi_g = S_g + \sum_{g'=1}^\infty \Sigma_{sg'g} \phi_{g'}$




### Direct Coupling

We split up the scattering term into four separate terms:

$\sum_{g'=1}^\infty \Sigma_{sg'g} \phi_{g'} = \sum_{g'=1}^{g-2} \Sigma_{sg'g} \phi_{g'} +\Sigma_{s,g-1,g}\phi_{g-1} + \Sigma_{sgg}\phi_g + \sum_{g'=g+1}^\infty \Sigma_{sg'g}\phi_{g'}$

It is common to consider neutron energies well above thermal energies, where upscattering can be neglected, and the last term vanishes. 

In the "direct coupling" approximation, it is also common to only consider scattering events that increase a neutron's group number by one, corresponding to the second term in the sum, and the first term vanishes. We also neglect self scattering for now, where a scattering event will leave a neutron's group number unchanged, corresponding to the third term. 

With these approximations, the multigroup slowing down equation becomes:

$\Sigma_{tg} \phi_g = S_g +  \Sigma_{s,g-1,g} \phi_{g-1}$

We also note that we can approximate the scattering cross section as:

$\Sigma_{s,g-1,g} = \frac{1}{\phi_{g-1}} \int_g dE \int_{g-1} dE' \Sigma_s(E'\rightarrow E) \phi(E') $

$\hspace{5ex}= \frac{1}{\phi_{g-1}} \int_{g-1} dE' \phi(E') \int_g dE \Sigma_s(E'\rightarrow E)$

$\hspace{5ex}\approx \frac{1}{\phi_{g-1}} \int_{g-1} dE' \phi(E') \int_0^\infty dE \Sigma_s(E'\rightarrow E)$

$\hspace{5ex}= \frac{1}{\phi_{g-1}} \int_{g-1} dE'\phi(E')\Sigma_s(E') $

$\hspace{5ex}\equiv \Sigma_{s,g-1}$
    
which is just the total scattering cross section integrated over group $g-1$, the analog of $\Sigma_{t,g-1}$ for only scattering events. Going from the second to the third line we used the same approximation as the direct coupling approximation, written slightly differently mathematically. We write the total scattering cross section as a sum over the groups of the differential scattering cross section:


$\Sigma_s(E') = \int_0^\infty dE \Sigma_s(E'\rightarrow E) = \sum_{g=1}^\infty \int_g dE \Sigma_s(E'\rightarrow E) $

$\hspace{5ex}= \sum_{g=1}^{g'-1} \int_g dE \Sigma_s(E'\rightarrow E )
    + \int_{g'} dE \Sigma_s(E'\rightarrow E) $

$\hspace{5ex}+ \int_{g+1} dE \Sigma_s(E'\rightarrow E )
    + \sum_{g=g'+2}^\infty \int_g dE \Sigma_s(E'\rightarrow E) $

The first term represents upscattering, the second is the self-scattering term, the third term is the direct-coupling term, and the fourth term represents larger energy loss. We neglect all but the third term, leaving:

$\Sigma_s(E') \approx \int_{g'+1} dE \Sigma_s(E'\rightarrow E )$

or,

$\Sigma_{s,g-1,g} \approx \Sigma_{s,g-1}$

This is the version of the slowing-down equation implemented in Version 1.0 of `ntransporter`. 


### Including Self Scattering


Now including the self-scattering term in the in-scattering term in the multigroup slowing down equation, we get:

$\Sigma_{tg} \phi_g = S_g +  \Sigma_{s,g-1,g} \phi_{g-1} + \Sigma_{sgg} \phi_g$


Before, we had approximated $\Sigma_{s,g-1,g}$ as $\Sigma_{s,g-1}$, but now we keep the self-scattering term:

$\Sigma_{s,g-1,g} = \frac{1}{\phi_{g-1}} \int_g dE \int_{g-1} dE' \Sigma_s(E'\rightarrow E) \phi(E') $

$\hspace{5ex}= \frac{1}{\phi_{g-1}} \int_{g-1} dE' \phi(E') \int_g dE \Sigma_s(E'\rightarrow E)  $

$\hspace{5ex}\approx \frac{1}{\phi_{g-1}} \int_{g-1} dE' \phi(E') \left[\int_0^\infty dE \Sigma_s(E'\rightarrow E) -  \int_{g-1} dE \Sigma_s(E'\rightarrow E)\right]$

$\hspace{5ex} = \Sigma_{s,g-1} - \Sigma_{s,g-1,g-1}$


The slowing down equation then becomes:

$(\Sigma_{tg} - \Sigma_{sgg})\phi_g  = (\Sigma_{s,g-1} - \Sigma_{s,g-1,g-1})\phi_{g-1} + S_g$


Some additional calculations were done outside of `ntransporter` including the self-scattering term, and assuming the ratio of $\Sigma_{sgg}$ and $\Sigma_{s,g,g+1}$ is a constant at all energies. The resulting ratio is parameterized with the value $q$, $\Sigma_{sgg} = q\Sigma_{sg}$. Then we calculate the following:

$\frac{\Sigma_{sgg}}{\Sigma_{s,g,g+1}} = \frac{q}{1-q} = \frac{A}{B}$

with the values of $A$ and $B$ given by:

|  | $\alpha<\beta^2$           | $\beta^2 < \alpha < \beta$ | $\beta < \alpha$ 
|---| ------------- | ------------- |-------|
|$A$| $\beta-\ln \beta -1 $  | $\beta-\ln \beta -1 $  | $(1-\alpha)\left(\ln \frac{\alpha}{\beta} - 1\right) - \ln \alpha$  |
|$B$| $(1-\beta)^2$  | $1 + \alpha-2\beta-\alpha\ln \frac{\alpha}{\beta^2}$  | $1-\alpha+\alpha\ln\alpha$ |


with $\alpha=\left(\frac{A-1}{A+1}\right)^2$ (for mass number $A$) and $\beta=E_g/E_{g-1}$.

Note that calculations including the self-scattering effects are much more accurate for most heavy materials than the approximation without used in `ntransporter` version 1.


### Evaluating Group Constants

The evaluation of most group constants requires knowledge of the flux $\phi(E)$. We use a zeroth-order approximation for fast groups of $\phi(E)\propto 1/E$, the derivation of which relies on the following approximations:

- Assume s-wave scattering: $\Sigma_s(E'\rightarrow E)=\frac{\Sigma_s(E')}{(1-\alpha)E'}$ when $E<E'<E/\alpha$

- Ignore absorption and inelastic scattering

- Consider energies well below source energies

These above assumptions result in an asymptotic solution $\phi(E)\propto \frac{1}{\Sigma_s(E) E}$ (note also this is the form used to derive the approximations for the self-scattering fraction $q$ above). We then neglect the variation in $\Sigma_s(E)$ over the group: $\phi(E)\propto \frac{1}{E}$.

For the thermal group, we assume a Maxwell-Boltzmann distribution at room temperature (0.0257 eV).
