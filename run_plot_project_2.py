#Program to run project 2 c++ code


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
    
    for line in infile_:
		if len(line.split()) == 2:
			eigenvalues.append(float(line.split()[0]))
		if line.startswith('Execution time'):
			time = float(line.split()[2])
    infile.close()

    return array(eigenvalues), time

"""
------------------------------------------------------------------------------------------
"""

N = [5, 10, 50, 100]
max_number_of_iterations = 10000
p_max = 5.

#Running code:

n = N[2]

#os.system('g++ -o project_2/main project_2/main.cpp -larmadillo -llapack -lblas')
#os.system('./project_2/main %s %s %s' %(n, max_number_of_iterations, p_max))

FileName = 'EigvalSolver_jacobi_pMax%s_nStep%s.txt' %(int(p_max), n)

eigval, t = read_file(FileName)
#print t
#print "---"
#print eigval[0:-1]

FirstEigenvalues = eigval[0:3]
print FirstEigenvalues

AnalyticEigval = [3, 7, 11]

#print abs(AnalyticEigval - FirstEigenvalues)

#rounding off

print round(FirstEigenvalues[0],3)
print round(FirstEigenvalues[1],3)
print round(FirstEigenvalues[2],3)




