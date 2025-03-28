#scp nkalline@euler.phys.cmu.edu:/home/widom/FilesForNedjma/Fe.cI2/Fe.cif ./crystals/
#conda activate venv
#conda deactivate
#ssh nkalline@euler.phys.cmu.edu

import ase.io.cif

crystal = ase.io.cif.read_cif("crystals/B2-I.cif")
print(crystal.symbols)