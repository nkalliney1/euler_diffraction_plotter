'''for element in crystal.atoms:
        for atom in crystal.atoms[element]:
            e = cmath.exp(complex(0,-2*math.pi*np.dot(v, atom)))
            SG += crystal.sfactors[element]*e'''
#label points
#for i, txt in enumerate(labels):
#    plt.annotate(txt, (two_theta[i], I_G[i]))
#point labels

marker=".", markersize=2

labels.append(str(h) + "," + str(k) + "," + str(l) + ",")
            
#plot:
#print out nonzero values
#for i in range(len(two_theta)):
#    if (I_G[i]>1 and two_theta[i]>0):   
#        print(labels[i])
#        print(two_theta[i])

#plot:
#print out nonzero values
#plot


#I(k)
sigma = 0.1
I_k += abs(SG)**2 * cmath.exp(-(mag_k-mag_G)**2/(2*sigma))
