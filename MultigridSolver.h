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

#define PI (4.0 * atan(1.0))

class MultigridSolver
{
private:
    /// the total level of multigrid method
    int _n;
    /// the length of each unit
    double _h;
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
    /// the level where MG is proceeding now
    int _nowlevel = 1;
    /// the pointer of restriction operator
    RestrictionOperator* _pRestrictOP;
    /// the pointer of interation operator
    InterpolationOperator* _pInterpolateOP;
    /// the start and end index of f/v proceeding in this step
    std::vector<int> _Idx;
    /// the weight of WeightedJacobi
    double _w = 2.0/3;
    /// the times of Relaxtion
    int _RlxTimes = 10;

public:
    MultigridSolver();
    MultigridSolver(int n, std::vector<double> f, std::vector<double> v, double u0 = 0, double u1 = 0, double tol = 1e-6, int maxstep = 30, std::string S1 = "FullWeighting", std::string S2 = "Linear");
    void SetGridLevel(int n);
    void SetRightSide(std::vector<double> f);
    void SetInitialGuess(std::vector<double> v);
    void SetBoundaryCond(double u0, double u1);
    void SetTolerance(double tol);
    void SetMaxStep(int maxstep);
    void SetRestrictionType(std::string S);
    void SetNowLevel(int nowlevel);
    void PrintInfo();
    RestrictionOperator* pRestrictOP();
    InterpolationOperator* pInterpolateOP();
    void UpdateIndex();
    void PrintIdx();
    void WeightedJacobi();
    void BottomSolve();
    void VCycle();
    void Solve();
    std::vector<double> ReturnSolution();
    

    
};


