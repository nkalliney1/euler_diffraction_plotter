#scp nkalline@euler.phys.cmu.edu:/home/widom/FilesForNedjma/Fe.cI2/Fe.cif ./crystals/
#conda activate venv
#conda deactivate
#ssh nkalline@euler.phys.cmu.edu
#/B2-I-AlCr-TiV

import ase.io.cif
from openeye import oechem

crystal = ase.io.cif.read_cif("crystals/B2-I-po.cif")
print(crystal.symbols)
print(crystal.get_tags())