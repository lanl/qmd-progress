#!/usr/bin/env python3

import numpy as np
from sdc_system import *
from sdc_ptable import ptable
import ctypes as ct
import os
import scipy.linalg as sp
from gpmd import *

        
def get_displ_born_charges(latticeVectors,symbols,atomTypes,coordsIn,verb):
    
    nats = len(atomTypes)
    bornDispl = np.zeros((nats,3,3))
    dspl = 1.0E-5
    field = np.zeros((3))
    coords = np.copy(coordsIn)
    err,charges_out,forces_out,dipole_saved, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
    for j in range(nats):
        coords[j,0] = coordsIn[j,0] + dspl
        err,charges_out,forces_out,dipole_out_p, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
        coords[j,0] = coordsIn[j,0] - 2*dspl
        err,charges_out,forces_out,dipole_out_m, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
        bornDispl[j,:,0] = (dipole_out_p - dipole_out_m)/(2*dspl)
        coords[j,0] = coordsIn[j,0] + dspl

    for j in range(nats):
        coords[j,1] = coordsIn[j,1] + dspl
        err,charges_out,forces_out,dipole_out_p,energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
        coords[j,1] = coordsIn[j,1] - 2*dspl
        err,charges_out,forces_out,dipole_out_m, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
        coords[j,1] = coordsIn[j,1] + dspl
        bornDispl[j,:,1] = (dipole_out_p - dipole_out_m)/(2*dspl)
    for j in range(nats):
        coords[j,2] = coordsIn[j,2] + dspl
        err,charges_out,forces_out,dipole_out_p, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
        coords[j,2] = coordsIn[j,2] - 2*dspl
        err,charges_out,forces_out,dipole_out_m, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
        coords[j,2] = coordsIn[j,2] + dspl
        bornDispl[j,:,2] = (dipole_out_p - dipole_out_m)/(2*dspl)

    
    bchMat = np.zeros((3,3))
    bornCharges = np.zeros((nats))
    for i in range(len(bornDispl[:,0,0])):
        bchMat[:,:] = bornDispl[i,:,:]
        E,Q = sp.eigh(bchMat)
        bornCharges[i] = sum(E)/3.0

    return bornCharges, bornDispl
    


###############################
## Main starts here
###############################

#******************************
#Check forces with no Field
#******************************
verb = 0
#Load the coordinates and atom types
latticeVectors,symbols,atomTypes,coords0 = read_pdb_file("H2O.pdb",lib="None",verb=True)

nats = len(atomTypes)
dspl = 0.00001
field = np.zeros((3))
field[0] = 0.0

#Get forces at coords0
err,charges_out,forces0,dipole_out, energyp = gpmd(latticeVectors,symbols,atomTypes,coords0,field,verb)

coordsP = np.copy(coords0)
coordsP[0,0] = coords0[0,0] + dspl
err,charges_out_pl,forces_out_pl,dipole_out, energyp = gpmd(latticeVectors,symbols,atomTypes,coordsP,field,verb)

coordsM = np.copy(coords0)
coordsM[0,0] = coords0[0,0] - dspl
err,charges_out_pl,forces_out_pl,dipole_out, energym = gpmd(latticeVectors,symbols,atomTypes,coordsM,field,verb)

print("Check: dE/dx,F",(energyp - energym)/(2*dspl),forces0[0,0])
exit(0)
#******************************
#Compute Born Charges using displacements 
#******************************

bornChDiag,bornChargesDspl = get_displ_born_charges(latticeVectors,symbols,atomTypes,coords0,verb)

for j in range(len(atomTypes)):
    print("Atom index",j,"","Mulliquen Charges:",charges_out[j])
    print("Born Charges:\n",bornChargesDspl[j,:,:])


exit(0)
