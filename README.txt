=======
FourState
=======
   FourState is a simple class to explore the dynamics in a kinetic network 
   involving 4 states with 2 intermediates. We use the method developed
   by Berezhkovskii, Hummer and Szabo to investigate this system and estimate
   the effects of mutations in the network.
   
   For details on the method, check: Berezhkovskii, Hummer and Szabo,
   Reactive flux and folding pathways in network models of coarse-grained 
   protein dynamics, J. Chem. Phys. 130, 205102 (2009),  
   http://dx.doi.org/10.1063/1.3139063.

Contributors
============
   This code has been written by David De Sancho.

Installation
============
   python setup.py install --user

Basic usage
============

   import FourState 

   bhs = fourstate.FourState("random") # initialize the four-state model randomly

   bhs.run_simul("output_name") # run simulation to see time-evolution of populations

   bhs.run_commit() # calculate committor

   bhs.fold_rate() # calculate folding rate

   bhs.do_dot() # plot the network in dot format for graphviz processing
