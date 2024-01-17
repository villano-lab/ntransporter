# `ntransporter/boundary_conditions/`

Directory for calculating infinite-slab neutron flux. In the full neutron transport code, this data will serve as boundary conditions to the system, hence the name of the directory.


## Executables

This directory builds one executable called `NT_BC` for calculating infinite slab flux using direct-coupling approximation without self-scattering. 

### `NT_BC`

`NT_BC` is called with the following signature:

```
./NT_BC output_file_base_path path_to_ntransporter_base material [ngroups=100]
```




**Outputs:**

`NT_BC` outputs one data file with the form

```
<output_file_base>_<material>_<ngroups>_BC_V2.dat
```

This file contains calculated infinite-slab fluxes. Each row looks like:

$g\hspace{5ex}E_g\hspace{5ex}\phi_g$

That is, the group number (including group zero, for which the flux is zero), the group lower bound, and the group flux, separated by spaces.


**Arguments:**

`output_file_base_path` : base of file path of output data files, including base of filename

`path_to_ntransporter_base` : path to top-level `ntransporter` directory, e.g., `/home/ajbiffl3/ntransporter` (useful for out of source builds or when running from a different directory)

`material` : name of material to calculate boundary conditions in

`ngroups` : (optional, default = 100) number of fast groups

> [!Note]
> In order to run, the corresponding files must exist: `<path_to_ntransporter_base>/cross_sections/data/V1/data_<material>_<ngroups>_20_xs.dat`, `<path_to_ntransporter_base>/cross_sections/data/V2/data_<material>_<ngroups>_dx.dat` and `<path_to_ntransporter_base>/sources/data/V1/data_<material>_<ngroups>_Sg.dat`. These paths to the needed cross section and source data files are hard-coded in and can be changed if desired

## Subdirectories


```
|- ntransporter/
|---- boundary_conditions/
|------- data/
|---------- V1/
|               output data from NT_BC v1.0
|---------- V2/
|               output data from NT_BC v2.0
|---------- SuperSim_estimates/
|               estimates of infinite slab flux from SuperSim
```

## Implementation


The zero-temperature multigroup slowing down equation is

$(\Sigma_{tg}-\Sigma_{sgg})\phi_g =S_g+\sum_{g'=1}^{g-1}\Sigma_{sg'g} \phi_{g'}$

This is very straightforwardly implemented with $\phi_0\equiv0$.


## SuperSim Estimates

The estimates of infinite-slab flux in different materials from SuperSim stored in `ntransporter/boundary_conditions/data/SuperSim_estimates/` were recorded using SuperSim simulations of 10^5 neutrons in a 1-km slab of the specified material. Every track of every neutron was summed up and contributed to the flux in the corresponding energy group with weights given by the track length. Overall normalization is given by the sum of the group source data. These data files also contain an additional $\delta\phi_g$ column giving the absolute uncertainty in the group flux due to statistical errors. 


