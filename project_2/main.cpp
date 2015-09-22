#include <iostream>
#include <cmath>
#include <armadillo>
#include <cstdlib>
#include <time.h>

using namespace std;
using namespace arma;

//Function to write data/output to file
void WriteToFile(const vec x, const vec solution, const int n, const double time_diag, string FileName)
{
    ofstream myfile;
        string filename = "linear_eq_solution_" + FileName + "_n" + to_string(n) + ".txt";
        myfile.open (filename);
        myfile << "Data:" << "  "<< "x" << "     " << "Solution" << endl;
        myfile << "Time calculating with " << FileName << ": "  << time_diag << " " << " seconds" << endl;
        myfile << "---------------------" << endl;
        for (int i=1; i<n+1; i++)
        {
            myfile << x[i] << "    " << solution[i] << endl;
        }
        myfile.close();
        cout << "Datafile done for n=" << n << endl;
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

int solve_eq_jacobi_rotation(mat& B, const int n_step, const int maxNumberOfIterations)
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
            cout << "Jacobi algorithm did not converge after " << maxNumberOfIterations << " iterations. Aborting!" << endl;
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

int main()
{
    //variables that can change for each run:
    const int n_step = 5;
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
    cout << "Armadillo lib. eigenvalue solver: Time for n_step="
         << n_step << ":  " << time_arma << " seconds" << endl;
//-------------------------------------------------------------
    //clocking the operations (only solve, not making file):
    clock_t start_jacobi, finish_jacobi; //declaring start and finish time
    start_jacobi = clock();

    //solving equations with jacobi rotation:
    int numberOfIterations = solve_eq_jacobi_rotation(B, n_step, maxNumberOfIterations);

    //stopping timer:
    finish_jacobi = clock();
    double time_jacobi = ( (finish_jacobi - start_jacobi)/((double)CLOCKS_PER_SEC ) );
    cout << "Jacobi rotation eigenvalue solver: Time for n_step="
         << n_step << ":  " << time_jacobi << " seconds" << endl;
    cout << "Total number of iterations with Jacobi rotation: " << numberOfIterations << endl;

//-------------------------------------------------------------
    //retriving eigenvalues, interested in the three first:
    vec eigval_jacobi_rot = B.diag();
    eigval_jacobi_rot = sort(eigval_jacobi_rot);

    eigval_arma.print();
    cout << "------" << endl;
    eigval_jacobi_rot.print();


    ///fix write to file. write all data? or just the first ones?
    /// how big n_step need to get three lowest eigvals to 4 leading digits?
    /// dependency og p_max? save/write to file?
    /// how many transformations (before 0) as function of n_step
    /// also write time to file as function of n_step for both solvers

//    cout << a[0] << endl;
//    cout << a[1] << endl;
//    cout << a[2] << endl;
//    cout << "------" << endl;
//    cout << eigval_arma[0] << endl;
//    cout << eigval_arma[1] << endl;
//    cout << eigval_arma[2] << endl;

    return 0;
}

