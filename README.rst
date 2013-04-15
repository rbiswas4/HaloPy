HaloPy 
==========

Author: Suman bhattacharya
=========

HaloPy is a package  built using Python that calculates the calculates the concentration- mass relation, the mass function and its universality (i.e. fsigma vs 1/sigma), given a halo output from the numerical simulations.


HaloPy is released under the MIT software liscense (see LISCENSE).

Prerequisites
=============
Python, SciPy, NumPy

Description
=============
The code assumes two sets of files- halo properties files which contain the integrated information (e.g. total particle in a halo) and halo profile file which contain the profile information of each halo as a function of radius of each halo.

The code also reads in the transfer function from a file outputted by CAMB.

If you have different input file formats you have to change the driver.py and cM.py accordingly

MF.py: reads in the halo particle counts(both SO and FOF) and computes the mass function for both the halo definition

fsigma.py: reads in the mass function and computes the universality of the mass function i.e. 1/sigma vs. f(sigma)  

growth_ode.py: calculates the growth factor and its log derivative wrt to the scale factor for a particular redshift (needed for the fsigma calc.)

cM.py: reads in the halo counts and the halo profiles, fit each halo using the NFW form, computes concentration for each halo, bin them up to compute c-M relation

params.py: contains the column number which contains fof, so halo counts in the halo properties file and radius, profile particle counts and overdensity in profile files.

driver.py: calls all the  modules to compute the mass function and the c-M relation.

Running
=========

download the package and place it anywhere you like. 
to run do:
python driver.py input.par
An example input.par file is given.

