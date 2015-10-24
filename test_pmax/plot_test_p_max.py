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
		if line.startswith('Eigenvector'):
			index_eigvec = i
		if line.startswith('Jacobi'):
			return [], [], 0

		i += 1
    infile.close()

    eigenvalues = []
    for k in range(index_eigval+1,index_eigvec-1):
		eigenvalues.append(float(all_lines[k]))

    eigenvector = []
    for m in range(index_eigvec+1,i):
		eigenvector.append(float(all_lines[m]))

    return array(eigenvalues), eigenvector, time

"""
------------------------------------------------------------------------------------------
"""
P = 24
p_max_list = linspace(3,8,P+1)
print len(p_max_list)
eigval_nr = 2
n_step = 150

diff_eigval = zeros(P)
analytic_eigval = [3,7,11]
i = 0

for p_max in range(P):
	eigval_jacobi, eigvec_jacobi_gs, time_jacobi = read_file("EigenValVecSolver_jacobi_pMax%s_nStep%s.txt" %(p_max,n_step))
	diff_eigval[i] = abs(eigval_jacobi[eigval_nr] - analytic_eigval[eigval_nr])
	i += 1
	print p_max

plot(diff_eigval,p_max_list)
show()


