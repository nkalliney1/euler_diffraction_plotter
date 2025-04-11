import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from diffraction_module import get_reciprocal_vectors, calculate_I, get_crystal, get_form_factor_array, get_j0_coeffs, get_j2_coeffs
from diffraction_module import get_moments, get_occupancies

name = sys.argv[1]
wavelength = float(sys.argv[2])

file_type = "c"
if "-v" in sys.argv:
    file_type = "v"
if "-p" in sys.argv:
    file_type = "p"

diffraction_type = ""
if "-xz" in sys.argv:
    diffraction_type = "xz"
elif "-xc" in sys.argv:
    diffraction_type = "xc"
elif "-n" in sys.argv:
    diffraction_type = "n"
elif "-nm" in sys.argv:
    diffraction_type = "nm"
elif "-m" in sys.argv:
    diffraction_type = "m"
else:
    raise ValueError("Need valid diffraction type")

p_on = False
t_on = False
l_on = False
partial_occupancy = False
occupancies = {}
if "-p" in sys.argv:
    lp_on = True
if "-l" in sys.argv:
    l_on = True
if "-lp" in sys.argv:
    l_on = True
    p_on = True
if "-t" in sys.argv:
    t_on = True
if "-po" in sys.argv:
    partial_occupancy = True
    occupancies = get_occupancies(name)


crystal = get_crystal(name, file_type)
form_factor_array = get_form_factor_array(crystal, diffraction_type, partial_occupancy, occupancies)
rvectors = get_reciprocal_vectors(crystal)
j0_coeffs = 0
j2_coeffs = 0
moments = 0
if diffraction_type == "nm" or diffraction_type == "m":
    j0_coeffs = get_j0_coeffs(crystal, partial_occupancy, occupancies)
    j2_coeffs = get_j2_coeffs(crystal, partial_occupancy, occupancies)
    moments = get_moments(name)

to_write = "h,k,l,mag(G),S,I,d,2theta\n"
two_theta = []
I_G = []
labels = []
I_G_g = []
two_theta_g = []

#diffraction calculations
lower = -9
upper = 9
for h in range(lower,upper):
    for k in range(lower,upper):
        for l in range(lower,upper):
            if (h == 0 and l == 0 and k == 0):
                continue
            
            G = h*rvectors[0]+k*rvectors[1]+l*rvectors[2]
            hkl = np.array([h, k, l])

            mag_k = 2*math.pi/wavelength
            mag_G = math.sqrt(np.dot(G, G))

            OH = 0.5*math.sqrt(np.dot(G, G))/mag_k
            #test if 2theta is in range
            if (OH < -1 or OH > 1):
                theta = -2
            else:
                theta = math.degrees(math.asin(OH))

            theta_r = math.radians(theta)
            L = 1
            T = 1
            P = 1
            if l_on:
                L = 1/((math.sin(theta_r)**2)*4*math.cos(theta_r))
            if p_on:
                P = 0.5*(1+math.cos(2*theta_r)**2)
            if t_on:
                B = 0.8
                T = math.exp(-B*(math.sin(theta_r)/wavelength)**2)
            

            I = calculate_I(diffraction_type, form_factor_array, crystal, hkl, wavelength, theta_r, j0_coeffs, j2_coeffs, moments, partial_occupancy, occupancies, L, P, T)
            #hkl
            to_write+=str(h) + "," + str(k) + "," + str(l) + ","
            #|G|
            to_write+=str(mag_G) + ","
            #intensity
            to_write+=str(I) + ","
            #interplanar spacing
            to_write+=str(2*math.pi/np.dot(G, G)) + ","
            #2theta
            to_write+=str(2*theta) + "\n"

            #for plotting sharp peaks
            if (2*theta in two_theta):
                I_G[two_theta.index(2*theta)] += I
            elif 2*theta > 0:
                two_theta.append(2*theta)
                I_G.append(I)
            
            #for plotting with gaussian
            mu = mag_G
            sigma = 0.01
            x = np.linspace(mu - 5*sigma, mu + 5*sigma, 101)
            for i in range(len(x)):
                if (0.5*x[i]/mag_k > 1 or 0.5*x[i]/mag_k<-1):
                    continue
                theta_2 = math.degrees(math.asin(0.5*x[i]/mag_k))
                if 2*theta_2 in two_theta_g:
                    I_G_g[two_theta_g.index(2*theta_2)] += I * math.exp(-(mag_G-x[i])**2/(2*sigma**2))
                elif 2*theta_2 > 0:
                    two_theta_g.append(2*theta_2)
                    I_G_g.append(I * math.exp(-(mag_G-x[i])**2/(2*(sigma**2)))) 
            
#save to file
with open("crystal_data/" + name+"_"+ diffraction_type +"_data.csv", "w") as f:
    f.write(to_write)

max = 0
for i in range(len(I_G_g)):
    if I_G_g[i] > max and two_theta_g[i] > 0 and two_theta_g[i] < 120:
        max = I_G_g[i]

for i in range(len(I_G_g)):
    I_G_g[i] = I_G_g[i]/max

plt.plot(two_theta_g, I_G_g)
plt.xlim(0, 120)
plt.ylim(0, 1.1)
plt.xlabel(r"2${\Theta}$ [deg]")
plt.ylabel(r"Intensity")
plt.title("Neutron diffraction for type-II AlCrVTi")
plt.show()