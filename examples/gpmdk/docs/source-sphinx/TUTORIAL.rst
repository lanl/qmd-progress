Tutorial
==========


Steered MD 
##############

Ad tut here


  .. code-block:: bash
 
    STMR0=


Applying bias 
##############

In this tutorial we will explain the steps to apply a voltage to specific atoms
in the system and look at how the DOS changes acordingly.

The coordinates are given in this pdb file: 

  .. include:: ./_static/water_solvent.pdb
    :literal:

This `pdb` file corresponds to bulk water. 
The file can be visualized using `vmd` as follows::

  vmd -f coords.pdb

A vmd spnpshot can be seen in the figure below:

   .. image:: ./_static/figures/water.png

We will run the code using an `input.in` file which has the following content:

  .. include:: ./_static/input_voltage.in
    :literal:

Make sure to have `ApplyVoltage=  T`. The definition of the input variables can be see in 
:ref:`input file section <input_file_choices>`. The DFTB parameters can be dounloaded 
from the `LATTE <https://github.com/lanl/LATTE/tree/master/TBparam>`_ github repository.
One can also use wget as follows::
        
        wget https://raw.github.com/lanl/LATTE/master/TBparam/electrons.dat
        wget https://raw.github.com/lanl/LATTE/master/TBparam/ppots.nonortho
        wget https://raw.github.com/lanl/LATTE/master/TBparam/bondints.nonortho
        mkdir TBparams
        cp raw.github.com/lanl/LATTE/master/TBparam/electrons.dat ./TBparams
        cp raw.github.com/lanl/LATTE/master/TBparam/ppots.nonortho ./TBparams
        cp raw.github.com/lanl/LATTE/master/TBparam/bondints.nonortho ./TBparams


Next, we provide a voltage.in file indicating which are the atom indices with applied bias, 
together with a bias column. This file content is as follows::

        3 
        1 1.0 
        2 1.0 
        3 1.0 

Which means we will shift the energy levels of the first water molecule by 1.0 eV. 
Finally, we run the code as::

  OMP_NUM_THREADS=1 mpirun -n 1 /qmd-progress/examples/gpmdk/build/gpmdk input_voltage.in

Because we have use the `StopAt= DMMin` this will just do a Density matrix energy minimization or full
self-consistency (SCF) calculation. 
The code will produe the `tdos_output` containing the Density Of State (DOS). If we do the 
calculation with and without bias and plot the DOS we will get the following: 


   .. image:: ./_static/figures/bias.png




