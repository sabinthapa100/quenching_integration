============================================================
 quenching_cuba_integration
 Version 0.1a 
 Copyright (c) Michael Strickland and Sabin Thapa
 Oct 28 2023
============================================================

------------------------------------------------------------
 REQUIREMENT
------------------------------------------------------------

This code requires Cuba library (https://feynarts.de/cuba/)
installed in your machine. This code performs multidimen-
-sional integration of our function using Cuhre, a 
deterministic algorithm based in cubature rules available in 
Cuba library.

You can change the global parameters in the input/params.txt 
file.

------------------------------------------------------------
 COMPILING
------------------------------------------------------------

To compile, simply type 

   $make 

from the main code directory.  This should generate an
executable called "quenching".

------------------------------------------------------------
 RUNNING  
------------------------------------------------------------

To run with the code, simply type

  $./quenching


------------------------------------------------------------
 OUTPUT
------------------------------------------------------------

The code will generate outputs into the "output" folder.

There is a Mathematica notebook in the "mathematica"
directory that can read in the output and generate the 
necessary plots.  These are then saved in the 
"mathematica" folder in Mathematica's ".m" format so that
they can be easily imported into other notebooks. 

------------------------------------------------------------
 LICENSE
------------------------------------------------------------

GNU General Public License (GPLv3)
See detailed text in license directory 
