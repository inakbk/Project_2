# Project_2
Project 2 in FYS-3150

Here you can find all the programs developed in project 2. I have used classes instead of functions for the jacobi algorithm, although they do the same (in my case). I copied and used, with only small changes, the same header files for all the projects in QT creator. I realized at a later stage that it is not usual to have the class inside the header file, this is my first time writing a class. To save time I decided to not use time on implementing the classes in a cpp file and only decleare the class in the header file. 

Did also not have time to implement unit testing, but did make my own test program for the jacobi algorithm.

About the folders:
(1e=one electron, 2e=two electrons without coulomc interaction, 2e_w_col=two electrons with coulomb interaction)

-The folder project_2_test_jacobi includes all files to test the implementation of the Jacobi algorithm with a random symetric matrix. All the other cpp files at built on this one.

-The folders starting with find_p_max run the algirithm for different p_max and plots to decide which to use

-The folder plotting_vs_nStep_1e takes the datafiles from the first runs with one electron and plots different parameters against n_step

-The folders project_2_1e, project_2_2e, project_2_2e_w_col are pretty much identical except the potential V_i. The programs solve the systems for the 1e and 2e case and plots the probability distrubution. 
