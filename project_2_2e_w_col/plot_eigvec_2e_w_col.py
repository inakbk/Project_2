#Program to plot eigenvectors for the two electron case with the coulomb interaction

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
        if line.startswith('Strength'):
            w_r = float(line.split()[6])        
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
n_step = 150
p_max = 6
p_min = 0;
h = (p_max - p_min)/n_step
p = linspace(p_min, p_max, n_step+1)
w_r_list = [0.01, 0.05, 1, 5]

for w_r in w_r_list:
    eigval_jacobi, eigvec_jacobi_gs, time_jacobi = read_file("files/EigenValVecSolver_jacobi_pMax%s_nStep%s_wr%s.txt" %(int(p_max), n_step, int(w_r*100)))
    eigvec_jacobi_gs = array([0] + eigvec_jacobi_gs + [0])
    abs_sq_psi_jacobi = eigvec_jacobi_gs*eigvec_jacobi_gs

    figure(1)
    plot(p,abs_sq_psi_jacobi)
    hold('on')

title('Plot of the probability distribution against p_max for different strength \n of the coulomb interaction', fontsize=16)
xlabel("p_max", fontsize=16)
ylabel("probability", fontsize=16)

legend(["w_r=%s" %w_r_list[0], "w_r=%s" %w_r_list[1], "w_r=%s" %w_r_list[2], "w_r=%s" %w_r_list[3]], fontsize=14)
        
"""
for w_r in [0.01, 0.05, 1, 5]:
    eigval_arma, eigvec_arma_gs, time_arma = read_file("files/EigenValVecSolver_arma_pMax%s_nStep%s_wr%s.txt" %(int(p_max), n_step, int(w_r*100)))
    eigvec_arma_gs = array([0] + eigvec_arma_gs + [0])
    abs_sq_psi_arma = eigvec_arma_gs*eigvec_arma_gs

    figure(2)
    plot(p,abs_sq_psi_arma) #they are identical??? (arma and jacobi)
    hold('on')
"""

show()





