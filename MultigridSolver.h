/**
 * @file   MultigridSolver.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Tue Oct 13 20:23:46 2020
 * 
 * @brief  solve 1D possion equation by multigrid method on interval [0,1].
 * 
 * 
 */

#include "RestrictionOperator.h"
#include "InterpolationOperator.h"
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include <string>

class MultigridSolver
{
private:
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
    int _nowstep = 0;
    /// the pointer of restriction operator
    RestrictionOperator* _pRestrictOP;
    /// the pointer of interation operator
    InterpolationOperator* _pInterpolateOP;
    

public:
    MultigridSolver();
    MultigridSolver(int n, std::vector<double> f, std::vector<double> v, double u0 = 0,
		    double u1 = 0, double tol = 1e-6, int maxstep = 30,
		    std::string S1 = "FullWeighting", std::string S2 = "Linear");
    void SetGridLevel(int n);
    void SetRightSide(std::vector<double> f);
    void SetInitialGuess(std::vector<double> v);
    void SetBoundaryCond(double u0, double u1);
    void SetTolerance(double tol);
    void SetMaxStep(int maxstep);
    void SetRestrictionType(std::string S);
    void PrintInfo();
    RestrictionOperator* pRestrictOP();
    InterpolationOperator* pInterpolateOP();

    
};


