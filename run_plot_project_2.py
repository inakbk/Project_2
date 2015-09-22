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

    x = []
    v = []
    
    for line in infile_:
		if len(line.split()) == 2:
			x.append(float(line.split()[0]))
			v.append(float(line.split()[1]))
		if len(line.split()) >= 4:
			time = float(line.split()[4])
    infile.close()

    return array(x), array(v), time

def u(x):
	analytic_solution = 1 - (1 - exp(-10))*x - exp(-10*x)
	return analytic_solution

"""
------------------------------------------------------------------------------------------
"""

N = [10, 100, 1000, 10**4]#, 10**5, 10**6, 10**7, 10**8] #the highest n only works for the subst. method

os.system('g++ -o project_2/main project_2/main.cpp -larmadillo -llapack -lblas')
os.system('./project_2/main')


