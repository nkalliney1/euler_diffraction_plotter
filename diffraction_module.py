import math
import cmath
import numpy as np
import ase.io.vasp
import ase.io.cif
import csv
from openeye import oechem

def get_cromer_f(form_factors, symbol, wavelength, theta_r):
    #get cromer form factors as a function of sine(theta)/wavelength
    nums = form_factors[symbol]
    form_factor = 0
    s = math.sin(theta_r)/wavelength
    for i in range(4):
        ai = float(nums[2*i])
        bi = float(nums[2*i+1])
        form_factor += ai * math.exp(-bi*(s**2))
    form_factor += float(nums[8])
    return form_factor

def get_magnetic_f(symbol, wavelength, theta_r, j0_coeffs, j2_coeffs):
    #get magnetic form factors as a function of sine(theta)/wavelength
    s = math.sin(theta_r)/wavelength
    j0 = 0
    j2 = 0
    if symbol in j0_coeffs:
        j0_atom = j0_coeffs[symbol]
        for i in range(3):
            j0 += float(j0_atom[i]) * math.exp(-1*float(j0_atom[i+1])*s*s)
        j0 += float(j0_atom[6])

    if symbol in j2_coeffs:
        j2_atom = j2_coeffs[symbol]
        for i in range(3):
            j2 += float(j2_atom[i]) * math.exp(-1*float(j2_atom[i+1])*s*s)
        j2 += float(j2_atom[6])
        j2 = j2*s*s

    f = j0 + (-2-3.82608552)*j2/2
    return f

def get_magnetic_moment(symbol, moments):
    return moments[symbol]*np.array([1,0,0])

def get_moments(name):
    #get all magnetic moments
    moments = {}
    with open("crystals/" + name + "-moments.csv", newline='') as f:
        lines = f.readlines()
        for line in lines:
            moments[line[0:2]] = float(line[3:])
    return moments

def get_reciprocal_vectors(crystal):
    #reciprocal lattice vectors
    reciprocal_vectors = []
    cell = crystal.cell
    volume = np.dot(cell[0], np.cross(cell[1], cell[2]))
    reciprocal_vectors.append(2*math.pi*np.cross(cell[1],cell[2])/volume)
    reciprocal_vectors.append(2*math.pi*np.cross(cell[2],cell[0])/volume)
    reciprocal_vectors.append(2*math.pi*np.cross(cell[0],cell[1])/volume)
    return reciprocal_vectors

def calculate_I_xz(crystal, v, partial_occupancy, occupancies):
    #calculate intensity for scattering w/o angular dependence

    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        if partial_occupancy:
            f = 0
            for atom in occupancies[crystal.symbols[i]]:
                f += atom[1]*oechem.OEGetAtomicNum(atom[0])
        else:
            f = crystal.numbers[i]
        SG += f*e
    return abs(SG)**2

def calculate_I_xc(form_factors, crystal, v, wavelength, OH, partial_occupancy, occupancies):
    #calculate intensity for scatting w cromer form factors
    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))

        if partial_occupancy:
            f = 0
            for atom in occupancies[crystal.symbols[i]]:
                f += atom[1]*float(get_cromer_f(form_factors,atom[0],wavelength,OH))
        else:
            f = float(get_cromer_f(form_factors,crystal.symbols[i],wavelength,OH))

        SG += f*e
    return abs(SG)**2

def calculate_I_n(form_factors, crystal, v, partial_occupancy, occupancies):
    #nuclear neutron scattering
    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))

        if partial_occupancy:
            f = 0
            for atom in occupancies[crystal.symbols[i]]:
                f += atom[1]*float(form_factors[atom[0]])
        else:
            f = float(form_factors[crystal.symbols[i]])

        SG += f*e
    return abs(SG)**2

def calculate_I_nm(form_factors, crystal, v, wavelength, OH, j0_coeffs, j2_coeffs, moments, partial_occupancy, occupancies):
    #nuclear and magnetic neutron scattering

    I = calculate_I_n(form_factors, crystal, v, partial_occupancy, occupancies)
    
    F_m = np.array([0,0,0])
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        if partial_occupancy:
            f = 0
            m = 0
            for atom in occupancies[crystal.symbols[i]]:
                f = atom[1]*float(get_magnetic_f(atom[0],wavelength,OH, j0_coeffs, j2_coeffs))
                m = atom[1]*get_magnetic_moment(atom[0], moments)
        else:
            f = float(get_magnetic_f(crystal.symbols[i],wavelength,OH, j0_coeffs, j2_coeffs))
            m = get_magnetic_moment(crystal.symbols[i], moments)

        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        F_m = np.add(F_m, f*e*m)
    k_hat = v / math.sqrt(np.dot(v, v))
    I += np.linalg.norm(np.cross(k_hat, np.cross(F_m, k_hat)))**2

    return I

def calculate_I(type, form_factors,crystal, v, wavelength, OH, j0_coeffs, j2_coeffs, moments, partial_occupancy, occupancies):
    if type == "xz":
        return calculate_I_xz(crystal, v, partial_occupancy, occupancies)
    elif type == "xc":
        return calculate_I_xc(form_factors, crystal, v, wavelength, OH, partial_occupancy, occupancies)
    elif type == "n":
        return calculate_I_n(form_factors, crystal, v, partial_occupancy, occupancies)
    elif type == "nm":
        return calculate_I_nm(form_factors, crystal, v, wavelength, OH, j0_coeffs, j2_coeffs, moments, partial_occupancy, occupancies)

def get_crystal(name, file_type):
    if file_type == "v":
        crystal = ase.io.vasp.read_vasp("crystals/"+name+".vasp")
    else:
        crystal = ase.io.cif.read_cif("crystals/"+name+".cif")
    return crystal

def get_form_factor_array(crystal, diffraction_type, partial_occupancy, occupancies):
    #get form factors depending on type of diffraction
    #only for nuclear neutron or cromer form factors
    
    symbols = []
    if partial_occupancy:
        for location in occupancies:
            for atom in occupancies[location]:
                symbols.append(atom[0])
    else:
        symbols = crystal.symbols

    form_factors = []
    if diffraction_type == "xc":
        with open("params/cromer_factors.csv", newline='') as f:
            cromer_form_factors = {}
            lines = f.readlines()
            for i in range(0,len(lines),2):
                line = lines[i][:-1]
                if line in symbols and line not in cromer_form_factors:
                    nums = lines[i+1][2:-1].split(" ")
                    cromer_form_factors[line]=[float(x) for x in nums]
            form_factors = cromer_form_factors

    elif diffraction_type == "n" or diffraction_type == "nm":
        with open("params/neutron_diffraction_lengths.csv", 'r') as file:
            all_neutron_form_factors = []
            csvreader = csv.reader(file)
            all_neutron_form_factors = next(csvreader)
            neutron_form_factors = {}

            for i in range(len(all_neutron_form_factors)):
                if (oechem.OEGetAtomicSymbol((i)) in symbols):
                    neutron_form_factors[oechem.OEGetAtomicSymbol((i))] = float(all_neutron_form_factors[i])
            form_factors = neutron_form_factors
    elif diffraction_type == "xz":
        return
    else:
        raise ValueError('Not a valid type of diffraction for this code')
    
    return form_factors

def get_j0_coeffs(crystal, partial_occupancy, occupancies):
    #used to calculate nuclear scattering form factor
    symbols = []
    if partial_occupancy:
        for location in occupancies:
            for atom in occupancies[location]:
                symbols.append(atom[0])
    else:
        symbols = crystal.symbols

    j0 = {}
    with open("params/j0.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            if line[0:2] in symbols and line[0:2] not in j0:
                nums = lines[i][4:-1].split("\t")
                j0[line[0:2]]=[float(x) for x in nums]
    return j0

def get_j2_coeffs(crystal, partial_occupancy, occupancies):
    #used to calculate nuclear scattering form factor
    symbols = []
    if partial_occupancy:
        for location in occupancies:
            for atom in occupancies[location]:
                symbols.append(atom[0])
    else:
        symbols = crystal.symbols

    j2 = {}
    with open("params/j2.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            if line[0:2] in symbols and line[0:2] not in j2:
                nums = lines[i][4:-1].split("\t")
                j2[line[0:2]]=[float(x) for x in nums]
    return j2

def get_occupancies(name):
    occupancies = {}

    with open("crystals/"+name+"-header.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            nums = lines[i][3:-1].split(" ")
            occupancies[line[0:2]] = []
            for j in range(0,len(nums),2):
                occupancies[line[0:2]].append([nums[j], float(nums[j+1])])
    return occupancies