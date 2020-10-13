/**
 * @file   MultigridSolver.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 13 20:23:46 2020
 * 
 * @brief  solve 1D possion equation by multigrid method on interval [0,1].
 * 
 * 
 */

#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

class MultigridSolver
{
public:
    /// the total level of multigrid method
    int _n;
    /// the collection of right side
    std::vector<double> _f;
    /// the collection of approximate solution
    std::vector<double> _v;
    /// the Dirichlet condition
    double _u0, _u1;
    /// the tolerance
    double _tol;
    /// the upper limit of iteration steps
    int _maxstep;
    /// the iteration steps have done
    int _nowstep;


    MultigridSolver();
    MultigridSolver(int n, std::vector<double> f, std::vector<double> v, double u0 = 0,
		    double u1 = 0, double tol = 1e-6, int maxstep = 30);
    SetGridLevel(int n);
    SetRightSide(std::vector<double> f);
    SetInitialGuess(std::vector<double> v);
    SetBoundaryCond(double u0, double u1);
    SetTolerance(double tol);
    SetMaxStep(int maxstep);
    

    
};


