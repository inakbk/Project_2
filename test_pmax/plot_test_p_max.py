#Program to run project 2 c++ code and plot datafiles

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
		if line.startswith('Jacobi'):
			return [], [], 0

		i += 1
    infile.close()

    eigenvalues = []
    for k in range(index_eigval+1,len(all_lines)):
		eigenvalues.append(float(all_lines[k]))

    eigenvector = []

    return array(eigenvalues), eigenvector, time

"""
------------------------------------------------------------------------------------------
"""
P = 14
p_max_list = linspace(3,6,P+1)
n_step = 100
N = [50, 100, 150, 200]
analytic_eigval = [3,7,11]

#P = 24
#p_max_list = linspace(3,8,P+1)
#n_step = 150

k = 0
for n_step in N:
	diff_eigval = zeros(P+1)
	for eigval_nr in [0,1,2]:
		i = 0
		for p_max in range(P+1):
			eigval_jacobi, eigvec_jacobi_gs, time_jacobi = read_file("figures_p_max/EigenValVecSolver_jacobi_pMax%s_nStep%s.txt" %(p_max,n_step))
			diff_eigval[i] = abs(eigval_jacobi[eigval_nr] - analytic_eigval[eigval_nr])
			i += 1
		figure(k)	
		plot(p_max_list, diff_eigval)
		title('n_step = %s' %n_step)
		xlabel("p_max")
		ylabel("diff_eigval")
		hold('on')

	plot(4.71428571,0.0005, "kx")
	legend(["eigval 1", "eigval 2", "eigval 3", "chosen p_max"])
	k += 1

show()


