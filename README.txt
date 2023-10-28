============================================================
 quenching_cuba_integration
 Version 0.1a 
 Copyright (c) Michael Strickland and Sabin Thapa
 Oct 28 2023
============================================================

------------------------------------------------------------
 REQUIREMENTS
------------------------------------------------------------

This code requires Cuba library (https://feynarts.de/cuba/)
installed in your machine. This code performs multidimen-
-sional integration of our function using Cuhre, a 
deterministic algorithm based on cubature rules available in 
Cuba library.


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

You can change the parameters in the input/params.txt file.

There is a file "input/params.txt" that contains all
parameters that can be adjusted at runtime.  It includes
comments describing the various options.  To run with the
params.txt parameters simply type

  ./quenching

If you would like to override some parameters in the
params.txt file from the commandline the syntax is

  ./quenching -<NAME1> <value1> ... -<NAMEn> <valuen>

------------------------------------------------------------
 OUTPUT
------------------------------------------------------------

The code will generate outputs into the "output" folder.

There is a Mathematica notebook in the "mathematica"
directory that can read in the output and generate the 
necessary plots.  These are then saved in the 
"mathematica" folder in both pdf and Mathematica's ".m" 
format so that they can be easily imported into other 
notebooks. 

------------------------------------------------------------
 LICENSE
------------------------------------------------------------

GNU General Public License (GPLv3)
See detailed text in license directory 
