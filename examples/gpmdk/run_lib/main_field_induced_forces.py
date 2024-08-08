#!/usr/bin/env python3

import numpy as np
from sdc_system import *
from sdc_ptable import ptable
import ctypes as ct
import os
import scipy.linalg as sp
from gpmd import *

def get_dipole(coords,charges):
    dipole = np.zeros(3)

    for i in range(len(charges)):
        dipole[0] = dipole[0] + charges[i]*coords[i,0]
        dipole[1] = dipole[1] + charges[i]*coords[i,1]
        dipole[2] = dipole[2] + charges[i]*coords[i,2]


def get_field_born_charges(latticeVectors,symbols,atomTypes,coords,verb):

    field = np.zeros(3) 
    err,charges_out,forces_saved,dipole_out, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
    nats = len(atomTypes)
    bornField = np.zeros((nats,3,3))
    for i in range(3):

        field = np.zeros(3)
        field[i] = -1.0E-4  
        err,charges_out,forces_out,dipole_out, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
        bornField[:,:,i] = (forces_out - forces_saved)/np.linalg.norm(field)

        field = np.zeros(3)
        field[i] = 1.0E-4
        err,charges_out,forces_out,dipole_out, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
        bornField[:,:,i] = (forces_out - forces_saved)/np.linalg.norm(field)


    bchMat = np.zeros((3,3))
    bornCharges = np.zeros((nats))
    for i in range(len(bornField[:,0,0])):
        bchMat[:,:] = bornField[i,:,:]
        E,Q = sp.eigh(bchMat)
        bornCharges[i] = sum(E)/3.0

    return bornCharges

        
def get_displ_born_charges(latticeVectors,symbols,atomTypes,coordsIn,verb):
    
    nats = len(atomTypes)
    bornDispl = np.zeros((nats,3,3))
    dspl = 1.0E-7
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
field[0] = 1.0

#Get forces at coords0
err,charges_out_pl,forces0,dipole_out, energyp = gpmd(latticeVectors,symbols,atomTypes,coords0,field,verb)

coordsP = np.copy(coords0)
coordsP[0,0] = coords0[0,0] + dspl
err,charges_out_pl,forces_out_pl,dipole_out, energyp = gpmd(latticeVectors,symbols,atomTypes,coordsP,field,verb)

coordsM = np.copy(coords0)
coordsM[0,0] = coords0[0,0] - dspl
err,charges_out_pl,forces_out_pl,dipole_out, energym = gpmd(latticeVectors,symbols,atomTypes,coordsM,field,verb)

print("Check: dE/dx,F",(energyp - energym)/(2*dspl),forces0[0,0])

#******************************
#Compute forces with field 
#******************************

field = np.zeros((3))
#lambdaF = 1000000.0
lambdaF = 0.000001
field[0] = -lambdaF
coordsM = np.copy(coords0)
err,charges,forcesm,dipole, energy = gpmd(latticeVectors,symbols,atomTypes,coordsM,field,verb)

field[0] = lambdaF
coordsP = np.copy(coords0)
err,charges,forcesp,dipole, energy = gpmd(latticeVectors,symbols,atomTypes,coordsP,field,verb)


bornCh = (forcesp - forcesm)/(2*lambdaF) 


#bornChDiag,bornChargesDspl = get_displ_born_charges(latticeVectors,symbols,atomTypes,coords0,verb)

for j in range(len(atomTypes)):
    #print("Ch,BChD,BChF",j,bornChargesDspl[j,0,0],bornCh[j,0])
    print("Ch,BChD,BChF",j,bornCh[j,0])
    print("Fp,Fm",j,forcesp[0,0]-forcesm[0,0])


exit(0)
#latticeVectors,symbols,atomTypes,coords = read_pdb_file("H2O.pdb",lib="None",verb=True)

#Call the general gpmd funtion
field = np.zeros((3))
dipole_sved = np.zeros((3))
charges_saved = np.zeros((nats))
field[0] = 0.0 ; field[1] = 0.001 ; field[2] = 0.0 
err,charges_out,forces_out,dipole_out, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
charges_saved[:] = charges_out[:]

forces_saved = np.zeros((nats,3))
forces_saved[:] = forces_out[:] 
print("Charge on atom 1",charges_out[0])
print("Dipole",dipole_out)
#dipole_out = get_dipole(coords,charges_out)
#print("Dipole",dipole_out)
print("Forces on atom 1",forces_out[0,:])
dipole_saved = np.zeros((3))
dipole_saved[:] = dipole_out[:] 

field[0] = 0.0 ; field[1] = -0.001 ; field[2] = 0.0 
err,charges_out,forces_out,dipole_out, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
print("charges_out",charges_out)

nats = len(coords[:,0])
bornField = np.zeros((nats,3))
bornField = forces_out - forces_saved
bornField = bornField/np.linalg.norm(field)
print("bornField",bornField[9,1],charges_out[9],charges_saved[9])
dx = 0.001
field[:] = 0.0
bornDispl = np.zeros((3))
deltaQ = np.zeros((nats))
sumBF = 0.0

for j in range(nats):
    coords[j,0] = coords[j,0] + dx 
    err,charges_out,forces_out,dipole_out, energy = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb) 
    mydip = np.zeros((3))
    for k in range(nats):
        mydip[:] = mydip[:] + (charges_out[k]-charges_saved[k])*coords[k,:]

    
    print("mydip",mydip)
    print("dipole",dipole_out)
    print("dipole_saved",dipole_saved)
    print("dipole_diff",dipole_out - dipole_saved)
    print("charges_diff",charges_out - charges_saved)
    bornDispl = (dipole_out - dipole_saved)/dx
    print("BornDispl",bornDispl[0])
    print("BornField",bornField[j,0])
    sumBF = sumBF + bornField[j,0]
    deltaQ = charges_out - charges_saved
    print("deltaQ",charges_out[j],deltaQ[j],deltaQ[j]*(coords[j,0] - dx) + charges_out[j])
    print("BornBoth",charges_out[j],bornDispl[0],bornField[j,0],bornDispl[0]/bornField[j,0])
    exit(0)
    coords[j,0] = coords[j,0] - dx 


print("SumBF",sumBF)
  
 #   tch = 0.0
 #   bchMat = np.zeros((3,3))
 #   for i in range(len(bornch_out[:,0])):
 #       bchMat[0,:] = bornch_out[i,0:3]
 #       bchMat[1,:] = bornch_out[i,3:6]
 #       bchMat[2,:] = bornch_out[i,6:9]
 #       E,Q = sp.eigh(bchMat) 
 #       print(i,"\nbc",bchMat)
 #       print("E",E)
 #       print("resultant charge",np.sum(E))
 #       print("Muliken charge",charges_out[i])


