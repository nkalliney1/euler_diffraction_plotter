import math
import cmath
import numpy as np
import ase.io.vasp
import ase.io.cif
import csv

def get_cromer_f(form_factors, number, symbol, wavelength, theta_r):
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

    f = j0 + (2-1.91304276)*j2/2
    return f

def get_magnetic_moment(symbol):
    return np.array([1,0,0])

def get_reciprocal_vectors(crystal):
    reciprocal_vectors = []
    cell = crystal.cell
    volume = np.dot(cell[0], np.cross(cell[1], cell[2]))
    reciprocal_vectors.append(2*math.pi*np.cross(cell[1],cell[2])/volume)
    reciprocal_vectors.append(2*math.pi*np.cross(cell[2],cell[0])/volume)
    reciprocal_vectors.append(2*math.pi*np.cross(cell[0],cell[1])/volume)
    return reciprocal_vectors

def calculate_I_xz(crystal, v):
    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        f = crystal.numbers[i]
        SG += f*e
    return abs(SG)**2

def calculate_I_xc(form_factors, crystal, v, wavelength, OH):
    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        f = float(get_cromer_f(form_factors,crystal.numbers[i],crystal.symbols[i],wavelength,OH))
        SG += f*e
    return abs(SG)**2

def calculate_I_n(form_factors, crystal, v):
    SG = 0
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        f = float(form_factors[crystal.numbers[i]-1])
        SG += f*e
    return abs(SG)**2

def calculate_I_nm(form_factors, crystal, v, wavelength, OH, j0_coeffs, j2_coeffs):
    I = calculate_I_n(form_factors, crystal, v)
    
    F_m = np.array([0,0,0])
    positions = crystal.get_scaled_positions()
    for i in range(len(positions)):
        e = cmath.exp(complex(0,-2*math.pi*np.dot(v, positions[i])))
        f = float(get_magnetic_f(crystal.symbols[i],wavelength,OH, j0_coeffs, j2_coeffs))
        F_m = np.add(F_m, f*e*get_magnetic_moment(crystal.symbols[i]))
    k_hat = v / math.sqrt(np.dot(v, v))
    I += np.linalg.norm(np.cross(k_hat, np.cross(F_m, k_hat)))**2

    return I

def calculate_I(type, form_factors,crystal, v, wavelength, OH, j0_coeffs, j2_coeffs):
    if type == "xz":
        return calculate_I_xz(crystal, v)
    elif type == "xc":
        return calculate_I_xc(form_factors, crystal, v, wavelength, OH)
    elif type == "n":
        return calculate_I_n(form_factors, crystal, v)
    elif type == "nm":
        return calculate_I_nm(form_factors, crystal, v, wavelength, OH, j0_coeffs, j2_coeffs)

def get_crystal(name, file_type):
    if file_type == "v":
        crystal = ase.io.vasp.read_vasp("crystals/"+name+".vasp")
    else:
        crystal = ase.io.cif.read_cif("crystals/"+name+".cif")
    return crystal

def get_form_factor_array(crystal,diffraction_type):
    #get form factors depending on type of diffraction
    form_factors = []
    if diffraction_type == "xc":
        with open("cromer_factors.csv", newline='') as f:
            cromer_form_factors = {}
            lines = f.readlines()
            for i in range(0,len(lines),2):
                line = lines[i][:-1]
                if line in crystal.symbols and line not in cromer_form_factors:
                    cromer_form_factors[line]=lines[i+1][2:-1].split(" ")
            form_factors = cromer_form_factors
    elif diffraction_type == "n" or diffraction_type == "nm":
        with open("neutron_diffraction_lengths.csv", 'r') as file:
            neutron_form_factors = []
            csvreader = csv.reader(file)
            neutron_form_factors = next(csvreader)
            form_factors = neutron_form_factors
    elif diffraction_type == "xz":
        form_factors.append("easter egg!")
    else:
        raise ValueError('Not a valid type of diffraction for this code')
    return form_factors

def get_j0_coeffs(crystal):
    j0 = {}
    with open("j0.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            if line[0:2] in crystal.symbols and line[0:2] not in j0:
                j0[line[0:2]]=lines[i][4:-1].split("\t")
    return j0

def get_j2_coeffs(crystal):
    j2 = {}
    with open("j2.csv", newline='') as f:
        lines = f.readlines()
        for i in range(0,len(lines)):
            line = lines[i]
            if line[0:2] in crystal.symbols and line[0:2] not in j2:
                j2[line[0:2]]=lines[i][4:-1].split("\t")
    return j2