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
    number_of_iterations = 0
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
N = linspace(10,200,20) 
p_max = 4.7
analytic_eigval = [3,7,11]

number_of_iterations = zeros(len(N))
time_jacobi = zeros(len(N))
time_arma = zeros(len(N))

diff_eigval1_jacobi = zeros(len(N))
#diff_eigval1_arma = zeros(len(N))
diff_eigval2_jacobi = zeros(len(N))
#diff_eigval2_arma = zeros(len(N))
diff_eigval3_jacobi = zeros(len(N))
#diff_eigval3_arma = zeros(len(N))

i = 0
for n_step in N:
	eigval_jacobi, eigvec_jacobi_gs, time_jacobi[i], number_of_iterations[i] = read_file("files_with_jacobi_arma_n10_to_200_pmax4_7/EigenValVecSolver_jacobi_pMax0_nStep%s.txt" %(int(n_step)))
	eigval_arma, eigvec_arma_gs, time_arma[i], tull_ball = read_file("files_with_jacobi_arma_n10_to_200_pmax4_7/EigenValVecSolver_arma_pMax0_nStep%s.txt" %(int(n_step)))	

	diff_eigval1_jacobi[i] = abs(eigval_jacobi[0] - analytic_eigval[0])
	#diff_eigval1_arma[i] = abs(eigval_arma[0] - analytic_eigval[0])
	diff_eigval2_jacobi[i] = abs(eigval_jacobi[1] - analytic_eigval[1])
	#diff_eigval2_arma[i] = abs(eigval_arma[1] - analytic_eigval[1])
	diff_eigval3_jacobi[i] = abs(eigval_jacobi[2] - analytic_eigval[2])
	#diff_eigval3_arma[i] = abs(eigval_arma[2] - analytic_eigval[2])
	i += 1

figure(1)
plot(N, number_of_iterations, "g")
title('Number of iterations versus n_step with p_max = %s' %p_max)
xlabel("n_step")
ylabel("number_of_iterations")

figure(2)
plot(N, time_jacobi, "b")
hold('on')
plot(N, time_arma, "r")
title('Plot of the execution time versus n_step with p_max = %s' %p_max)
xlabel("n_step")
ylabel("time")
legend(["jacobi", "arma"], loc='upper left')

figure(3)	
plot(N, diff_eigval1_jacobi, "r")
hold("on")
plot(N, diff_eigval2_jacobi, "k")
plot(N, diff_eigval3_jacobi, "b")
title('Diff eigval with p_max = %s' %p_max)
xlabel("n_step")
ylabel("diff_eigval")
legend(["eigval 1", "eigval 2", "eigval 3"])

show()























