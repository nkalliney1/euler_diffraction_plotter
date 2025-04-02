# Diffraction Plotter

## How to use the code

### Factors under consideration
### Types of diffraction
This code can create diffraction plots for four types of diffraction. The type of diffraction must be specified by the user with a flag.

(a) x-ray diffraction without angular dependence (flag: `-xz`)

(b) x-ray diffraction with angular dependence (cromer form factors) (flag: `-xc`)

(c) nuclear scattering (flag: `-n`)

(d) magnetic nuclear scattering (flag: `-nm`)


### File types
This code can read crystal structures in `.cif` and `.vasp` files. The default is `.cif` files, but `.vasp` files can be used by adding the flag `-v`.

### Extra factors
Diffraction is impacted by a temperature factor. By default this is not included, but can be added to calculations if desired by including the flag `-t`.

Unpolarized x-ray diffraction is also impacted by an angular dependence called the Lorentz-Polarization factor. By default this is not included, but can be added to calculations if desired by including the flag `-lp`. This can also be added to other types of diffraction but is not physically meaningful.

### Wavelength
The wavelength of the scattering x-ray or neutron must be specified by the user.

### Partial occupancy
If the crystal has partial occupancy, add the flag `-po`.

### Command structure
When running the code, the terminal command follows this structure:
```
python3 diffraction.py crystal_name wavelength diffraction_type optional_flags
```
The crystal name should be given without the file extension, and wavelength should be in Angstroms.

Example command:
This command produces the magnetic neutron diffraction plot for a crystal titled BCC with partial occupancy. The wavelength of the incident neutrons is 1.54 A.
```
python3 diffraction.py BCC 1.54 -nm -po
```

## Files required
All crystals must be described using a `.cif` or `.vasp` file.

## Partial occupancy
Crystals that are described with partial occupancy should be structured like so:
In the `.cif` file, each location should be filled with a placeholder element, which does not necessarily have to be one of the elements in the crystal, but must be unique.
The user should create a second file, titled `crystalname-header.csv` that should be structured like so:
```
placeholder_element_name element_1 occupancy_fraction_1 element_2 occupancy_fraction_2 [etc]
```
for as many lines as necessary. Lines can be in any order.

For example, to describe AlCrVTi as a BCC lattice where the center position is occupied by V and Ti at equal fractions, and the vertex occupied by Al and Cr at equal fractions, and using Al as the placeholder for the vertices and Ti as the placeholder for the center:
End of `.cif` file:
```
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_symmetry_multiplicity
_atom_site_Wyckoff_label
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_fract_symmform
Al1 Al  1 a 0.00000 0.00000 0.00000 1.00000 0,0,0 
Ti1 Ti  1 b 0.50000 0.50000 0.50000 1.00000 0,0,0 
```

Header file:
```
Al Al 0.5 Cr 0.5
Ti Ti 0.5 V 0.5
```

## Magnetic neutron diffraction
For magnetic neutron diffraction, all atoms with magnetic moments should have their moments specified in any order n a file named `crystalname-moments.csv`. Moments should be in terms of the Bohr magneton. 
For example, for a crystal where Fe is the only atom with a magnetic moment, `Fe-moments.csv` is simply:
```
Fe 2
```