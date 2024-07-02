#!/usr/bin/env python3

import numpy as np
from sdc_system import *
from sdc_ptable import ptable
import ctypes as ct
import os

## General gpmd call function 
# @brief This will run gpmd and get back forces and charges, ect. 
# from a full SCF optimization. Other output could be added.
# @param LatticeVectors. 3x3 matrix containing the lattice vectors for the simulation box.
# latticeVectors[1,:] = first lattice vector.
# @param symbols Symbols for each atom type, e.g., the element symbol of the first atom is symbols[types[0]] 
# @param atomTypes Type for each atom, e.g., the first atom is of type "types[0]"
# @param coords for each atom, e.g., z-coordinate of the frist atom is coords[0,2]
# @param verb Verbosity level for gpmd output
#
def gpmd(latticeVectors,symbols,atomTypes,coords,verb):
 
    # Import the shared library
    gpmdLibFileName = os.environ['GPMD_PATH'] + '/libgpmd.so'

    try:
        gpmdLib = ct.CDLL(gpmdLibFileName)
        f = gpmdLib.gpmd_compute
    except:
        print("ERROR loading the gpmd library. Please the check that the \
            GPMD_PATH is set correctly and that there is a compiled libgpmd.so")
        exit(0)

    #Periodic table: We use this to pass the chemical atom types as integer instead of characters.
    pt = ptable()

    nats = len(coords[:,1])

    #Getting atomic numbers
    nTypes = len(symbols)
    atomicNumbers = np.zeros((nTypes),dtype=np.int32)
    atomTypes32 = np.zeros((nats),dtype=np.int32)
    atomTypes32[:] = atomTypes
    for i in range(len(symbols)):
        atomicNumbers[i] = pt.get_atomic_number(symbols[i])


    #Vectorizing 2D arrays for C-Fortran interoperability
    coordsFlat_in = np.zeros(3*nats) #Vectorized coordinates
    forcesFlat_out = np.zeros(3*nats) #Vectorized forces
    chargesFlat_out = np.zeros(nats) #We call this one Flat just for consistency
    dipoleFlat_out = np.zeros(3) #Same here 
    bornchFlat_out = np.zeros(9*nats) 

    for i in range(nats):
        coordsFlat_in[3*i] = coords[i,0]
        coordsFlat_in[3*i+1] = coords[i,1]
        coordsFlat_in[3*i+2] = coords[i,2]

    latticeVectorsFlat = np.zeros((9))
    latticeVectorsFlat[0] = latticeVectors[0,0]
    latticeVectorsFlat[1] = latticeVectors[0,1]
    latticeVectorsFlat[2] = latticeVectors[0,2]

    latticeVectorsFlat[3] = latticeVectors[1,0]
    latticeVectorsFlat[4] = latticeVectors[1,1]
    latticeVectorsFlat[5] = latticeVectors[1,2]

    latticeVectorsFlat[6] = latticeVectors[2,0]
    latticeVectorsFlat[7] = latticeVectors[2,1]
    latticeVectorsFlat[8] = latticeVectors[2,2]

    #Specify arguments as a pointers to pass to Fortran
    f.argtypes=[ct.c_int,ct.c_int,ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),\
            ct.POINTER(ct.c_int),ct.POINTER(ct.c_int),ct.POINTER(ct.c_double),\
            ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.POINTER(ct.c_double),ct.c_int]

    #Inputs
    coords_ptr = coordsFlat_in.ctypes.data_as(ct.POINTER(ct.c_double))
    atomTypes_ptr = atomTypes32.ctypes.data_as(ct.POINTER(ct.c_int))
    atomicNumbers_ptr = atomicNumbers.ctypes.data_as(ct.POINTER(ct.c_int))
    latticeVectors_ptr = latticeVectorsFlat.ctypes.data_as(ct.POINTER(ct.c_double))

    #Outputs
    charges_ptr = chargesFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))
    forces_ptr = forcesFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))
    dipole_ptr = dipoleFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))
    bornch_ptr = bornchFlat_out.ctypes.data_as(ct.POINTER(ct.c_double))

    #Call to the fortran funtion
    err = f(ct.c_int(nats),ct.c_int(nTypes),coords_ptr,latticeVectors_ptr,\
            atomTypes_ptr,atomicNumbers_ptr,charges_ptr,forces_ptr,\
            dipole_ptr,bornch_ptr,ct.c_int(verb))

    charges_out = np.zeros(nats)
    charges_out[:] = chargesFlat_out[:]
    dipole_out = np.zeros(3)
    dipole_out[:] = dipoleFlat_out[:]
    
    #Back to a 2D array for the forces
    forces_out = np.zeros((nats,3))
    bornch_out = np.zeros((nats,9))
    for i in range(nats):
        forces_out[i,0] = forcesFlat_out[i*3 + 0]
        forces_out[i,1] = forcesFlat_out[i*3 + 1]
        forces_out[i,2] = forcesFlat_out[i*3 + 2]

        bornch_out[i,0] = bornchFlat_out[i*9 + 0]
        bornch_out[i,1] = bornchFlat_out[i*9 + 1]
        bornch_out[i,2] = bornchFlat_out[i*9 + 2]
        
        bornch_out[i,3] = bornchFlat_out[i*9 + 3]
        bornch_out[i,4] = bornchFlat_out[i*9 + 4]
        bornch_out[i,5] = bornchFlat_out[i*9 + 5]
        
        bornch_out[i,6] = bornchFlat_out[i*9 + 6]
        bornch_out[i,7] = bornchFlat_out[i*9 + 7]
        bornch_out[i,8] = bornchFlat_out[i*9 + 8]


    return err, charges_out, forces_out, dipole_out, bornch_out



