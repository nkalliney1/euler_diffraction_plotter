# Diffraction Plotter

## How to use the code

## Factors under consideration
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

## Partial occupancy
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