#scp nkalline@euler.phys.cmu.edu:/home/widom/FilesForNedjma/Fe.cI2/Fe.cif ./crystals/
#P1-I.Al-V-Cr-Ti
#P1-II.Al-Cr-V-Ti

#conda activate venv
#conda deactivate
#ssh nkalline@euler.phys.cmu.edu
#/B2-I-AlCr-TiV

import matplotlib.pyplot as plt

x = []
y = []

with open("to_plot.csv", "r") as f:
    lines = f.readlines()
    for i in range(1,len(lines)):
        line = lines[i].split("\t")
        x.append(float(line[0]))
        y.append(float(line[1]))

        
plt.plot(x[0:7132], y[0:7132], label = "P1-I")
plt.plot(x[7132:], y[7132:], label = "PI-II")
plt.xlim(0, 120)
plt.legend()
plt.xlim(0, 120)
plt.ylim(0, None)
plt.xlabel(r"2${\Theta}$ [deg]")
plt.ylabel(r"Intensity")
plt.show()