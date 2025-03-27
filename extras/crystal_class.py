import math
import cmath
import numpy as np

class Crystal:
    def __init__(self, name):
        with open(name+".vasp", newline='') as crystalfile:
            lines = crystalfile.readlines()
            
            #get multiplier
            self.multiplier = float(lines[1][3:19])

            #get translation vectors
            translation_vectors = []
            for i in range(2,5):
                line = lines[i]
                v = np.array([float(line[5:21]), float(line[26:42]), float(line[47:63])])
                translation_vectors.append(v)

            self.tvector = translation_vectors

            #calculate reciprocal lattice vectors
            reciprocal_vectors = []
            volume = np.dot(self.tvector[0], np.cross(self.tvector[1], self.tvector[2]))
            reciprocal_vectors.append(2*math.pi*np.cross(self.tvector[1],self.tvector[2])/volume)
            reciprocal_vectors.append(2*math.pi*np.cross(self.tvector[2],self.tvector[0])/volume)
            reciprocal_vectors.append(2*math.pi*np.cross(self.tvector[0],self.tvector[1])/volume)
            self.rvectors = reciprocal_vectors

            #get atoms and number
            atom_list = []
            num_list = []
            cur = ""
            for i in range(3,len(lines[5])):
                if lines[5][i] == " ":
                    if cur != "":
                        atom_list.append(cur)
                    cur = ""
                elif lines[5][i] == "\n":
                    if cur != "":
                        atom_list.append(cur)
                else:
                    cur += lines[5][i]

            cur = ""
            for i in range(3,len(lines[6])):
                if lines[6][i] == " ":
                    if cur != "":
                        num_list.append(int(cur))
                    cur = ""
                elif lines[6][i] == "\n":
                    if cur != "":
                        num_list.append(int(cur))
                else:
                    cur += lines[6][i]
            
            atoms_dict = {}
            current_line = 8
            for i in range(len(atom_list)):
                position_list = []
                for j in range(num_list[i]):
                    line = lines[current_line]
                    v = float(line[2:20])*self.tvector[0]/np.linalg.norm(self.tvector[0])+float(line[22:40])*self.tvector[1]/np.linalg.norm(self.tvector[1])+\
                        float(line[42:60])*self.tvector[2]/np.linalg.norm(self.tvector[2])
                    position_list.append(v)
                    current_line += 1
                atoms_dict[atom_list[i]] = position_list

            self.atoms = atoms_dict    

            #set structure factors
            structure_factors= {}
            for i in range(len(atom_list)):    
                structure_factors[atom_list[i]] = 1
            
            self.sfactors = structure_factors