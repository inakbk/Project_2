#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <time.h>

using namespace std;
using namespace arma;

//Function to write data/output to file
void WriteToFile(const vec& eigenvalues, double p_max, const int n_step, const int number_of_iterations, const double time, string FileName, const int converge_test)
{
    ofstream myfile;
        string filename = "eigval_solver_" + FileName + "_n_step" + to_string(n_step) + ".txt";
        myfile.open (filename);
        myfile << "Solution of the eigenvalueproblem for the " << FileName << " algorithm." << endl;
        myfile << "Dimention of matrix, one less than n_step: " << n_step << endl;
        myfile << "Value of p_max: " << p_max << endl;
        myfile << "Execution time: " << time << endl;
        if(number_of_iterations != -1){
            myfile << "Number of iterations for jacobi algoritm: " << number_of_iterations << endl;
        }
        //Writing eigenvalues to file if the jacobi method did converge.
        if(converge_test == True)
        {
            myfile << "---------------------" << endl;
            myfile << "Eigenvalues (sorted)" << "  "<< "" << "     " << "Eigenvectors (soon)" << endl;
            myfile << "---------------------" << endl;
            int number_of_eigenvalues_printed = 10;
            for (int i=0; i < number_of_eigenvalues_printed; i++)
            {
                if(i == size(eigenvalues,0))
                {
                    cout << "Length of eigenval vec is shorter than number_of_eigenvalues_printed for " << FileName << " solver, exiting loop." << endl;
                    break;
                }
                myfile << eigenvalues[i] << "    " << "--" << endl;
            }
            myfile << "(Maximum writing " << number_of_eigenvalues_printed << " eigenvalues to file.)" << endl;
            myfile.close();

            cout << "Datafile done for n_step=" << n_step << " with the " << FileName <<" solver." << endl;
            cout << endl;
        }
        //Writing error message til file if jacobi method did not converge.
        if(converge_test == False)
        {
            myfile << endl; myfile << endl;
            myfile << "Jacobi method did not converge after " << number_of_iterations << " iterations!" << endl;
            myfile << endl; myfile << endl;

            cout << "Datafile done for n_step=" << n_step << " with the " << FileName <<" solver which did not converge." << endl;
            cout << endl;
        }


}

void find_max_elem_index(int& k, int& l, double& max_off_diagonal, const mat &B, const int n_step)
{
    max_off_diagonal = -1.0;
    //Checking all off-diagonal elements:
    for(int i=0; i<n_step-1; i++)  {
        for(int j=i+1; j<n_step-1; j++) {
            //storing value and index of max off-diag-element
            if(fabs(B(i,j)) > max_off_diagonal)
            {
                max_off_diagonal = fabs(B(i,j));
                k = i;
                l = j;
            }
        }
    }
    if(max_off_diagonal < 0) {
        cout << "Warning, find_max_elem_index did not find any offdiagonal element with absolute value larger than zero. Aborting!" << endl;
        exit(1);
    }
}

void transformation_matrix(double& c, double& s, const mat &B, const int k, const int l)
{
    if(B(k,l) != 0)
    {
        double t = 0;
        double tau = (B(l,l) - B(k,k))/(2.0*B(k,l));
        // ensuring that theta<=pi/4:
        if(tau>0)
        {
            t = -tau + sqrt(1 + (tau*tau));
        }
        else
        {
            t = -tau - sqrt(1 + (tau*tau));
        }
        c = 1.0/sqrt(1.0 + (t*t));
        s = t*c;
    }
    else
    {
        c = 1.0;
        s = 0.0;
        cout << "Element chosen is zero: B(" << k << "," << l << ") = 0 ==> c = 1 and s = 0" << endl;
    }
}

void jacobi_rotation(mat& B, const double c, const double s, const int k, const int l, const int n_step)
{
    double B_kk = B(k,k);
    double B_ll = B(l,l);
    B(k,k) = B_kk*c*c - 2*B(k,l)*c*s + B_ll*s*s;
    B(l,l) = B_ll*c*c + 2*B(k,l)*c*s + B_kk*s*s;
    B(k,l) = 0.0;//(B_kk - B_ll)*c*s + B(k,l)*(c*c - s*s); //or just set to 0
    B(l,k) = B(k,l); //symetrix matrix

    //changing the remaining elements:
    for(int i=0; i<n_step-1; ++i)
    {
        if( (i!=k) && (i!=l) )
        {
            double B_ik = B(i,k);
            double B_il = B(i,l);
            B(i,k) = B_ik*c - B_il*s;
            B(k,i) = B(i,k);
            B(i,l) = B_il*c + B_ik*s;
            B(l,i) = B(i,l);
        }
    }
}

int solve_eq_jacobi_rotation(mat& B, const int n_step, const int maxNumberOfIterations, int converge_test)
{
    double tolerance = 1.0e-08;
    double max_off_diagonal = 0;
    int k = 0;
    int l = 0;
    find_max_elem_index(k, l, max_off_diagonal, B, n_step); //initial value for max_off_diagonal to enter loop

    int numberOfIterations = 0;
    while(tolerance < max_off_diagonal)
    {
        if(++numberOfIterations > maxNumberOfIterations) {
            cout << "Jacobi algorithm did not converge after " << maxNumberOfIterations << " iterations for n_step= " << n_step << ". Aborting!" << endl;
            converge_test = 0;
            exit(1);
        }

        //finding the value and index(k,l) of the maximum element in B:
        find_max_elem_index(k, l, max_off_diagonal, B, n_step);

        //finding the values of c ans s (the S transformation matrix):
        double c = 0;
        double s = 0;
        transformation_matrix(c, s, B, k, l);
        //transformation of B:
        jacobi_rotation(B, c, s, k, l, n_step);

    }
    return numberOfIterations;
}

int main(int argc, char *argv[])
{
    if(argc == 1)
    {
        cout << "No arguments. Give 'n' on command line. Eks n=10: 10" << endl;
        exit(1);
    }
    else{
        //variables from command line:
        const int n_step = atof(argv[1]);

        //make these command line too
        int maxNumberOfIterations = 10000;
        const double p_max = 5; //writing p instead of rho

    //-------------------------------------------------------------
        const double p_min = 0;
        const double h = (p_max - p_min)/n_step;
        vec p = linspace(p_min, p_max, n_step+1); //p_i = p_min + i*h
        vec V = p%p;

        // Constructing matrix B:
        double e = -1/(h*h); // all elements of the e vec is the same
        vec d = 2/(h*h) + V;
        mat B = zeros<mat>(n_step-1,n_step-1);
        for(int i=0, j=1; (i<=n_step-2) && (j<=n_step-2); ++i, ++j)
        {
            B(i,j) = e;
            B(i,i) = d[i+1];
            B(j,i) = e;
        }
        B(n_step-2,n_step-2) = d[n_step-1];

    //-------------------------------------------------------------
        //clocking the operations (only solve, not making file):
        clock_t start_arma, finish_arma; //declaring start and finish time
        start_arma = clock();

        //solving equations with armadillo lib:
        vec eigval_arma = eig_sym(B);

        //stopping timer:
        finish_arma = clock();
        double time_arma = ( (finish_arma - start_arma)/((double)CLOCKS_PER_SEC ) );
        //cout << "Armadillo lib. eigenvalue solver: Time for n_step="
        //     << n_step << ":  " << time_arma << " seconds" << endl;

        WriteToFile(eigval_arma, p_max, n_step, -1, time_arma, "arma");
    //-------------------------------------------------------------
        //clocking the operations (only solve, not making file):
        clock_t start_jacobi, finish_jacobi; //declaring start and finish time
        start_jacobi = clock();

        //solving equations with jacobi rotation:
        converge_test = 1; //initializing to true (false if it does not converge)
        int numberOfIterations = solve_eq_jacobi_rotation(B, n_step, maxNumberOfIterations, converge_test);

        //stopping timer:
        finish_jacobi = clock();
        double time_jacobi = ( (finish_jacobi - start_jacobi)/((double)CLOCKS_PER_SEC ) );
    //    cout << "Jacobi rotation eigenvalue solver: Time for n_step="
    //         << n_step << ":  " << time_jacobi << " seconds" << endl;
        cout << "Total number of iterations with Jacobi rotation: " << numberOfIterations << endl;

        //retriving eigenvalues, interested in the three first:
        vec eigval_jacobi_rot = B.diag();
        eigval_jacobi_rot = sort(eigval_jacobi_rot);

        WriteToFile(eigval_jacobi_rot, p_max, n_step, numberOfIterations, time_jacobi, "jacobi_rot");
    //-------------------------------------------------------------


        //eigval_arma.print();
        //cout << "------" << endl;
        //eigval_jacobi_rot.print();


        /// fix write to file. write all data? or just the first ones?
        /// how big n_step need to get three lowest eigvals to 4 leading digits?
        /// dependency og p_max? save/write to file?
        /// how many transformations (before 0) as function of n_step
        /// also write time to file as function of n_step for both solvers
        /// if did not converge write that to file!!

    //    cout << a[0] << endl;
    //    cout << a[1] << endl;
    //    cout << a[2] << endl;
    //    cout << "------" << endl;
    //    cout << eigval_arma[0] << endl;
    //    cout << eigval_arma[1] << endl;
    //    cout << eigval_arma[2] << endl;

    }

    return 0;
}

