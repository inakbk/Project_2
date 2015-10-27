#For two electrons with coulomb interaction: This code is for checking the choice of p_max with 4 different n_step, plotting 

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
#n_step = 100
N = [50, 100, 150, 200]
#analytic_eigval = [3,7,11] #is unknown for the 2e case

for n_step in N:
	#diff_eigval = zeros(P+1)
	eigval_jacobi_list = zeros(P+1)
	for eigval_nr in [0,1,2]:
		i = 0
		for p_max in range(P+1):
			eigval_jacobi, eigvec_jacobi_gs, time_jacobi = read_file("plot_files/EigenValVecSolver_jacobi_pMax%s_nStep%s.txt" %(p_max,n_step))
			eigval_jacobi_list[i] = eigval_jacobi[eigval_nr]
			#diff_eigval[i] = abs(eigval_jacobi[eigval_nr] - analytic_eigval[eigval_nr])
			i += 1	
		figure(eigval_nr)
		plot(p_max_list, eigval_jacobi_list)
		title('Plot of the %s eigenvalue against p_max for different n_step' %(eigval_nr+1), fontsize=16)
		xlabel("p_max", fontsize=16)
		ylabel("eigenvalue", fontsize=16)
		hold('on')

	#plot(2.215,0.1, "ko")
legend(["n_step=%s" %N[0], "n_step=%s" %N[1], "n_step=%s" %N[2], "n_step=%s" %N[3]], fontsize=14)

show()



