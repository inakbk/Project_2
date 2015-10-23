#Program to run project 2 c++ code and plot datafiles

from pylab import *
import os as os

"""
------------------------------------------------------------------------------------------
"""
#leser filen, lager arrays med dataen:
def read_file(filename):
    infile = open(filename, "r")
    infile_ = infile.readlines()

    eigenvalues = []
    eigenvectors = []
    
    for line in infile_:
		if line.startswith('Execution time'):
			time = float(line.split()[2])
		if line.startswith('Eigenvalues'):
			eigenvalues.append(float(line))
			if line.startswith('--'):
				break
		if line.startswith('Eigenvector'):
			eigenvectors.append(float(line))
			if line.startswith('--'):
				break
    infile.close()

    return array(eigenvalues), array(eigenvectors), time

"""
------------------------------------------------------------------------------------------
"""
N = [5]#, 10, 50, 100]
max_number_of_iterations = 10000
p_max = 5.

#Running code:
n_step = N[0]
os.system('g++ -o project_2_1e/main.cpp -larmadillo -llapack -lblas')
os.system('./project_2_1e/main %s %s %s' %(n_step, max_number_of_iterations, p_max))




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

