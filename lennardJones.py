def lennardJones(r, Eo, sigma):
        lj = (4*Eo*((sigma/r)**12 - (sigma/r)**6))
        return lj 