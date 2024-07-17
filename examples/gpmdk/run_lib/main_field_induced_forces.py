#!/usr/bin/env python3

import numpy as np
from sdc_system import *
from sdc_ptable import ptable
import ctypes as ct
import os
import scipy.linalg as sp
from gpmd import *


for i in range(1):
    #GPMD verbosity level.
    verb = 0
    #Load the coordinates and atom types
    #latticeVectors,symbols,atomTypes,coords = read_xyz_file("H2.xyz",lib="None",verb=True)
    latticeVectors,symbols,atomTypes,coords = read_pdb_file("benzene-CN.pdb",lib="None",verb=True)
    #latticeVectors,symbols,atomTypes,coords = read_pdb_file("H2O.pdb",lib="None",verb=True)
    print("Processing molecule H2O")
    
    #Call the general gpmd funtion
    field = np.zeros((3))
    field[0] = 0.0 ; field[1] = 0.0 ; field[2] = 0.0 
    err,charges_out,forces_out,dipole_out = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)

    print("Charge on atom 1",charges_out[0])
    print("Dipole",dipole_out)
    print("Forces on atom 1",forces_out[0,:])

    field[0] = 0.01 ; field[1] = 0.0 ; field[2] = 0.0 
    err,charges_out,forces_out,dipole_out = gpmd(latticeVectors,symbols,atomTypes,coords,field,verb)
    print("Forces on atom 1",forces_out[0,:])

  
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


