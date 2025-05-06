.. _input_file_choices:

Input file choices
==================


In this section, every input file keyword is described. Every valid keyword should be formatted using camel-case syntax and must directly precede an equal ``=`` sign. For illustration, the syntax ``JobName= MyJob`` represents a valid keyword. Comments are required to begin with a hash ``#`` symbol, positioned immediately adjacent to the text intended for the comment. An example of such a comment is: ``#My comment``.
The input file is structured into various sections, each initiated by a capitalized keyword followed by a left-curly brace, and concluded with a right-curly brace. An example is provided for the LATTE section::

        LATTE{
                AKeyWord= value
         }


LATTE section
-------------
This section encompasses all the keywords pertinent to the LATTE Hamiltonian, along with several parameters associated with electronic structure calculations.
In order to open this section use the following systax: ``LATTE{ KeyWord= Value   }``.

`JobName=`
***********
This variable will indicate the name of the job we are sunning. 
It is just a tag to distinguish different outputs. 
As we mentioned before and example use should be: ``JobName= MyJob``

`Verbose= 0`
*************
Controls the verbosity level of the output. If set to ``0`` no 
output is printed out. If set to ``1``, only basic messages of 
the current execution point of the code will be printed. 
If set to ``2``, information about basic quantities are also 
printed. If set to ``3``, all relevant possible info is printed.

`BMLType= Dense`
****************
The code will employ a specific sparse matrix storage format when interfacing with the BML library. The default format utilized is ``Dense``, which is a two-dimensional Fortran array. Alternatively, the formats ELLPACK and CSR are also available for selection.


`CoordsFile= "coords.xyz"`
***************************
This file will include spatial coordinates and atom types, with compatibility for both xyz and pdb file formats.

`Method= DiagEfFull`
********************
This option determines the overall methodology used in the calculation of the electronic structure. Among the methodologies available are Diag, DiagEfFull, and SP2.

`Mdim= -1`
**********
This parameter defines the maximum row length for storing BML ELLPACK matrices. In cases where the system size is substantial, it is advisable to decrease this parameter. Typically, the value of MDim corresponds to the maximum number of neighboring elements present. Should MDim be assigned a value of -1, it will be configured to match the total number of atoms.

`StopAt=`
***********
he specified keyword serves as a critical marker to introduce designated halts within the code execution process. By assigning the value ``gpmdcov_Energ``, the program is instructed to carry out a single-point full self-consistent field (SCF) calculation. Conversely, when set to ``gpmdcov_DM_Min``, the execution is limited to the computation of the Density matrix.

`MPulay= 10`
************
Dimension of the Pulay matrix employed in the mixing algorithm to enhance the convergence rate of the self-consistent field (SCF) calculations.

`ZMat= Diag`
************
The methodology for calculating the inverse overlap depends on the selected strategy. When utilizing the ``Diag`` option, the process involves diagonalization of the overlap matrix. In contrast, selecting the ``ZSP`` option involves employing a recursive procedure backed by matrix multiplication techniques.

`PulayCoeff= 0.1`
*******************
A mixing coefficient employed within the framework of the Pulay Mixing methodology.

`MixCoeff= 0.2`
***************
Coefficient of mixing implemented within the framework of the Extended Lagrangian Molecular Dynamics approach.

`SCFTol= 1.0d-5`
****************
The tolerance criterion for the self-consistent field (SCF) calculation is defined with values expressed in electron units.

`MaxSCFIter= 100`
******************
Maximum SCF iterations.

`CoulAcc= 1.0d-5`
*****************
Precision in Coulombic Computation.

`TimeRatio= 10.0`
*****************
The proportion of the Real and Reciprocal Ewald summations expressed as a function of the length cutoff.

`TimeStep= 0.5`
****************
The temporal increment utilized for the Molecular Dynamics simulation is expressed in femtoseconds.

`MDSteps= 100`
***************
Total number of steps for the MD simulation

`ParamPath= "./TBparams"`
*************************
Location where the parameters of the Hamiltonian are found.

`NlistEach= 10`
****************
The quantity of molecular dynamics steps executed prior to the reevaluation of the neighbor list

`MuCalcType= FromParts`
***********************
Methodology for Calculating the Chemical Potential: The two strategies that have been implemented are FromPart and Dyn.

`EFermi= -0.0`
**************
Preliminary value of the electronic chemical potential (Fermi level)

`kBT= 0.025`
************
The electronic temperature expressed in electron volts

`Entropy= T`
************
Should the calculation of entropy attributable to the electronic temperature be necessary. This is set to T when kBT is larger than 0.

`DoKernel= T`
*************
In the event that the Extended Lagrangian method is implemented utilizing the kernel approach.


GPMD Section 
-------------
The GPMD section encompasses keywords that are designated for high-level procedural types and calculations. This section is initiated by the command ``GPMD{ MyKeyWord= Value }``.

`ApplyVoltage= F`
*****************
This keyword serves to facilitate the application of electrical potential across certain components within the system.

`VoltageFile= "voltage.dat"`
****************************
This file comprises the atoms influenced by the application of voltage. It is structured with the initial row denoting the total number of affected atoms, followed by a column of atom indices coupled with their respective applied voltage values::

    5 #Number of biased atoms
    1 1.0
    2 1.0
    3 -1.0
    4 -1.0
    5 -1.0
