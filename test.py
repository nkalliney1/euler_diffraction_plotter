#scp nkalline@euler.phys.cmu.edu:/home/widom/FilesForNedjma/Fe.cI2/Fe.cif ./crystals/
#scp nkalline@euler.phys.cmu.edu:/home/widom/FilesForNedjma/B2-I-AlCr-TiV/B2.cif ./crystals/
#P1-I.Al-V-Cr-Ti
#P1-II.Al-Cr-V-Ti
#scp -r ./XNDC nkalline@euler.phys.cmu.edu:./XNDC
#/Users/nedjmakalliney/Desktop/Programs/XNDC

#conda activate venv
#conda deactivate
#ssh nkalline@euler.phys.cmu.edu
#/B2-I-AlCr-TiV

import matplotlib.pyplot as plt
import numpy as np
import math

'''
x = []
y = []
y1 = []
y2 = []

name = sys.argv[1]
type = sys.argv[2]

with open("crystals/"+name+"-mag.dat", "r") as f:
    lines = f.readlines()
    for i in range(1,len(lines)):
        line = lines[i].split(" ")
        x.append(float(line[0]))
        y1.append(float(line[1]))

with open("crystals/"+name+"-nonmag.dat", "r") as f:
    lines = f.readlines()
    for i in range(1,len(lines)):
        line = lines[i].split(" ")
        y2.append(float(line[1]))

for i in range(len(y1)):
    y.append(y1[i] - y2[i])

plt.title("GSAS")
if type == "-m":    
    plt.plot(x, y)
if type == "-n":    
    plt.plot(x, y2)
if type == "-nm":    
    plt.plot(x, y1)
plt.xlim(0, 120)
plt.ylim(0, None)
plt.xlabel(r"2${\Theta}$ [deg]")
plt.ylabel(r"Intensity")
plt.show()

'''

x = np.linspace(0, 8, num=50)
y = []

symbol = "V3"

j0 = {}
with open("params/j0.csv", newline='') as f:
    lines = f.readlines()
    for i in range(0,len(lines)):
        line = lines[i]
        if line[0:3] == symbol:
            nums = lines[i][5:-1].split("\t")
            print(nums)
            j0[line[0:3]]=[float(z) for z in nums]

        if line[0:2] == symbol:
            nums = lines[i][4:-1].split("\t")
            print(nums)
            j0[line[0:2]]=[float(z) for z in nums]

print(j0)

for val in x:
    j0_atom = j0[symbol]
    y_val = 0
    for i in range(3):
        y_val += j0_atom[2*i] * math.exp(-j0_atom[2*i+1]*val*val/(16*math.pi*math.pi))
    y_val += float(j0_atom[6])
    y.append(y_val)

plt.plot(x, y)
plt.grid()
plt.show()