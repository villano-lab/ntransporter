# `ntransporter/sources/`

Directory for calculating evaluated group source constants in specified material

## Executables

### `NT_Src`

This subdirectory generates the executable `NT_Src`. Its call signature is:

```
./NT_Src output_file_base_path material path_to_supersim [ngroups=100]
```

**Outputs:**

The executable creates one output data file, with the form

```
<output_file_base_path>_<material>_<ngroups>_Sg.dat
```

The file contains the evaluated group source constants in three columns. Each row looks like:

$g\hspace{5ex}E_g\hspace{5ex}S_g$

The group number (starting with group zero, for which the source term is zero), the group lower bound, and the group source separated by spaces.

**Arguments:**

`output_file_base_path` : base of file path of output data files, including base of filename

`material` : name of material to calculate sources for. In Version 1, accepted values are:
- `Norite` (SNOLAB Norite rock)
- `G4_Si` (pure amorphous silicon)
- `G4_Pb` (solid lead)
- `G4_POLYETHYLENE` (high-density polyethylene)
- `BoratedPoly` (borated polyethylene)
- `G4_Cu` (solid copper)

> [!Note] 
> In this version of `ntransporter`, all materials use the norite source files and weights merely for purposes of comparison with SuperSim

`path_to_supersim` : path to base `supersim` directory of SuperSim, e.g., `/home/ajbiffl3/supersim`

`ngroups` : (optional, default = 100) number of fast groups


## Subdirectories

The directory substructure under `sources` is the following:

```
|- ntransporter/
|---- sources/
|------- data/
|---------- V1/
|               output data from NT_Src v1.0
```

## Source Calculations

The group sources are evaluated with radiogenic neutron spectrum calculations from [`SOURCES4A`](https://www.osti.gov/biblio/15215)[^1], stored in data files within SuperSim (in `supersim/CDMSsources/spectra/neutron/`) for different materials.
[^1]: Note the version of `S4` used was likely a modified version provided by Prof. Vitaly Kudryavstev of the University of Sheffield, but we are still untangling all the details

These data files contain (energy, spectrum) information with the "energy" column in MeV and the "spectrum" column representing the total emission rate (neutrons per second per cubic centimeter) in the corresponding bin (with each bin spanning from the energy value in the same row down to the previous energy value). It is assumed the lower bound of the first bin is zero. To return these values to a differential quantity, each spectrum value is divided by its corresponding bin width.

The trapezoidal integral over the data is then taken over each bin, linearly interpolating between points to estimate the value at the group edges. The value over each group is then the total emission rate in each group (neutrons per second per cubic centimeter) for one ppb of each radioactive contaminant. The values are then multiplied by measured contaminant concentrations in the material.


### `Norite`

In this version of `ntransporter`, all calculations used the `S4` predictions for norite, which comes from the following two files:


- `norite_2013_U_1ppb.dat` : net spectrum (($\alpha$,n) and SF) from 1 ppb uranium contamination in norite
- `norite_2013_Th_1ppb.dat` : net spectrum (($\alpha$,n) and SF) from 1 ppb thorium contamination in norite


Measured contamination levels from a 2005 report are given in Table 6 of [this Confluence page](https://confluence.slac.stanford.edu/pages/viewpage.action?pageId=383932067):

|Isotope|Concentration (ppb)|
|----|---|
|U|1095|
|Th|5715|







