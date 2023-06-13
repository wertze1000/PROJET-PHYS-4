#Default values (to get a Lennard Jones potential close to the project's Figure 1)
#Eo = 0.45 eV
#sigma = 0.45 nm

def lennardJones(r, Eo, sigma):
        lj = (4*Eo*((sigma/r)**12 - (sigma/r)**6))
        return lj 