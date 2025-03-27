import math
import numpy as np
import matplotlib.pyplot as plt
import sys
from diffraction_module import get_reciprocal_vectors, calculate_I, get_crystal, get_form_factor_array, get_j0_coeffs, get_j2_coeffs

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
else:
    raise ValueError("Need valid diffraction type")

lp_on = True
t_on = True
if "-nlp" in sys.argv:
    lp_on = False
if "-nt" in sys.argv:
    t_on = False

crystal = get_crystal(name, file_type)
form_factor_array = get_form_factor_array(crystal, diffraction_type)
rvectors = get_reciprocal_vectors(crystal)
j0_coeffs = 0
j2_coeffs = 0
if diffraction_type == "nm":
    j0_coeffs = get_j0_coeffs(crystal)
    j2_coeffs = get_j2_coeffs(crystal)

to_write = "h,k,l,mag(G),S,I,d,2theta\n"
two_theta = []
I_G = []
labels = []
I_G_g = []
two_theta_g = []

#diffraction calculations
lower = -5
upper = 5
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
            LP = 1
            T = 1
            if lp_on:
                LP = (1 + (math.cos(2*theta_r)**2))/(8*(math.sin(theta_r)**2)*math.cos(theta_r))
            if t_on:
                B = 1.2
                T = math.exp(-B*(math.sin(theta_r)/wavelength)**2)
            

            I = calculate_I(diffraction_type, form_factor_array, crystal, hkl, wavelength, theta_r, j0_coeffs, j2_coeffs)
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
                I_G[two_theta.index(2*theta)] += I*LP*T
            else:
                two_theta.append(2*theta)
                I_G.append(I*LP*T)
            
            #for plotting with gaussian
            mu = mag_G
            sigma = 0.01
            x = np.linspace(mu - 5*sigma, mu + 5*sigma, 101)
            for i in range(len(x)):
                if (0.5*x[i]/mag_k> 1 or 0.5*x[i]/mag_k<-1):
                    continue
                theta_2 = math.degrees(math.asin(0.5*x[i]/mag_k))
                if 2*theta_2 in two_theta_g:
                    I_G_g[two_theta_g.index(2*theta_2)] += I * math.exp(-(mag_G-x[i])**2/(2*sigma**2))*LP*T
                else:
                    two_theta_g.append(2*theta_2)
                    I_G_g.append(I * math.exp(-(mag_G-x[i])**2/(2*(sigma**2)))*LP*T) 
            
#save to file
with open("/crystal_data/" + name+"_"+ diffraction_type +"_data.csv", "w") as f:
    f.write(to_write)

plt.plot(two_theta_g, I_G_g)
plt.xlim(0, 120)
plt.ylim(0, None)
plt.xlabel(r"2${\Theta}$ [deg]")
plt.ylabel(r"Intensity")
plt.show()