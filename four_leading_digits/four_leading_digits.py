#This code is to check how large n_step needs to be to get the 3 first eigenvalues with four leading digits

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
		if line.startswith('Eigenvalues'):
			index_eigval = i
		if line.startswith('Number'):
			number_of_iterations = float(line.split()[-1])
		if line.startswith('Jacobi'):
			return [], [], 0

		i += 1
    infile.close()

    eigenvalues = []
    for k in range(index_eigval+1,len(all_lines)):
		eigenvalues.append(float(all_lines[k]))

    eigenvector = []

    return array(eigenvalues), eigenvector, time, number_of_iterations

"""
------------------------------------------------------------------------------------------
"""
N = [50, 100, 150, 200] 
p_max = 4.7

analytic_eigval = [3,7,11]

k = 3
diff_eigval = zeros(len(N))
number_of_iterations = zeros(len(N))
for eigval_nr in [0,1,2]:
	i = 0
	for n_step in N:
		eigval_jacobi, eigvec_jacobi_gs, time_jacobi, number_of_iterations_list[i] = read_file("figures_p_max/EigenValVecSolver_jacobi_pMax0_nStep%s.txt" %(n_step))
		diff_eigval[i] = abs(eigval_jacobi[eigval_nr] - analytic_eigval[eigval_nr])
		i += 1
	
	figure(k)	
	plot(N, diff_eigval)
	title('p_max = %s' %p_max)
	xlabel("n_step")
	ylabel("diff_eigval")
	hold('on')
	legend(["eigval 1", "eigval 2", "eigval 3"])
	
	figure(eigval_nr)
	plot(N, number_of_iterations_list)
	title('p_max = %s' %p_max)
	xlabel("n_step")
	ylabel("number_of_iterations")
	hold('on')
	legend(["eigval 1", "eigval 2", "eigval 3"])
	k += 1

show()

"""
#must do for arma to compare time
	figure(k+3)
	plot(N, time_jacobi)
	title('p_max = %s' %p_max)
	xlabel("n_step")
	ylabel("number_of_iterations")
	hold('on')
	legend(["eigval 1", "eigval 2", "eigval 3"])
"""
