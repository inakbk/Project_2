#Program to plot eigenvectors

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
N = [150]#, 10, 50, 100]
max_number_of_iterations = 100000
p_max = 4.7

#Running code:
n_step = N[0]
# not working
#os.system('g++ -o project_2_1e/main project_2_1e/main.cpp -larmadillo -llapack -lblas')
#os.system('./project_2_1e/main %s %s %s' %(n_step, max_number_of_iterations, p_max))

eigval_arma, eigvec_arma_gs, time_arma = read_file("files/EigenValVecSolver_arma_pMax4_nStep%s.txt" %n_step)
eigvec_arma_gs = array([0] + eigvec_arma_gs + [0])
abs_sq_psi_arma = eigvec_arma_gs*eigvec_arma_gs

eigval_jacobi, eigvec_jacobi_gs, time_jacobi = read_file("files/EigenValVecSolver_jacobi_pMax4_nStep%s.txt" %n_step)
eigvec_jacobi_gs = array([0] + eigvec_jacobi_gs + [0])
abs_sq_psi_jacobi = eigvec_jacobi_gs*eigvec_jacobi_gs

p_min = 0;
h = (p_max - p_min)/n_step;
p = linspace(p_min, p_max, n_step+1); #p_i = p_min + i*h

plot(p,abs_sq_psi_jacobi)
#hold('on')
#plot(p,abs_sq_psi_arma)
show()
#they are identical???


"""
N = [5]#, 10, 50, 100]
max_number_of_iterations = 10000
p_max = 5.

#Running code:
#n_step = N[2]

for n_step in N:

	#Running code
	#os.system('g++ -o project_2/main project_2/main.cpp -larmadillo -llapack -lblas')
	#os.system('./project_2/main %s %s %s' %(n_step, max_number_of_iterations, p_max))

	FileName = 'EigvalSolver_jacobi_pMax%s_nStep%s.txt' %(int(p_max), n_step)

	eigval, t = read_file(FileName)

	#FirstEigenvalues = eigval[0:3]
	FirstEigenvalues = [2.999999, 6.499999999, 11.44444445]
	print FirstEigenvalues

	AnalyticEigval = [3, 7, 11]


	for i in range(3):
		print FirstEigenvalues[i]
		if FirstEigenvalues[i] > AnalyticEigval[i]:
			if (round(FirstEigenvalues[i],4) - AnalyticEigval[i]) < 0.5:
				print "yay"
				#print "---"
			else:
				print "not good enough precition"
				#print FirstEigenvalues[i] - AnalyticEigval[i]
				#print "---"
		if FirstEigenvalues[i] < AnalyticEigval[i]:
			if (AnalyticEigval[i] - round(FirstEigenvalues[i],4)) <= 0.5:
				#print AnalyticEigval[i] - FirstEigenvalues[i]
				print "yay"
				#print "---"
			else:
				print "not good enough precition"
				#print AnalyticEigval[i] - FirstEigenvalues[i]
				#print "---"

"""

