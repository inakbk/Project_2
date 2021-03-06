#Program to plot eigenvectors for the two electron case

from pylab import *
import os as os

"""
------------------------------------------------------------------------------------------
"""
#leser filen, lager arrays med dataen:
def read_file(filename):
    infile = open(filename, "r")
    all_lines = infile.readlines()

    i = 0
    for line in all_lines:
        if line.startswith('Execution time'):
            time = float(line.split()[2])
        if line.startswith('Value'):
            p_max = float(line.split()[3])        
        if line.startswith('Eigenvalues'):
            index_eigval = i
        if line.startswith('Eigenvector'):
            index_eigvec = i
        i += 1
    infile.close()

    eigenvalues = []
    for k in range(index_eigval+1,index_eigvec-1):
		eigenvalues.append(float(all_lines[k]))

    eigenvector = []
    for m in range(index_eigvec+1,i):
		eigenvector.append(float(all_lines[m]))

    return array(eigenvalues), eigenvector, time, p_max

"""
------------------------------------------------------------------------------------------
"""
N = [150]#, 10, 50, 100]

n_step = N[0]

eigval_arma, eigvec_arma_gs, time_arma, p_max = read_file("files/EigenValVecSolver_arma_pMax9_nStep%s.txt" %n_step)
eigvec_arma_gs = array([0] + eigvec_arma_gs + [0])
abs_sq_psi_arma = eigvec_arma_gs*eigvec_arma_gs

eigval_jacobi, eigvec_jacobi_gs, time_jacobi, p_max = read_file("files/EigenValVecSolver_jacobi_pMax9_nStep%s.txt" %n_step)
eigvec_jacobi_gs = array([0] + eigvec_jacobi_gs + [0])
abs_sq_psi_jacobi = eigvec_jacobi_gs*eigvec_jacobi_gs

p_min = 0;
h = (p_max - p_min)/n_step
p = linspace(p_min, p_max, n_step+1)

plot(p,abs_sq_psi_jacobi)
#hold('on')
#plot(p,abs_sq_psi_arma) #they are identical??? (arma and jacobi)
title('The probability distrubution 2e, n_step = %s' %n_step, fontsize=16)
xlabel("rho", fontsize=16)
ylabel("probability", fontsize=16)
show()





