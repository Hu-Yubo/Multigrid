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
    /// the level where MG is proceeding now
    int _nowlevel = 1;
    /// the pointer of restriction operator
    RestrictionOperator* _pRestrictOP;
    /// the pointer of interation operator
    InterpolationOperator* _pInterpolateOP;
    /// the start and end index of f or v proceeding in this step
    std::vector<int> _Idx;
    /// the weight of WeightedJacobi
    double _w = 2.0/3;
    /// the times of Relaxtion
    int _RlxTimes = 10;
    /// the keywords of cycle type
    std::string _TypeofCycle;
    /// the real solution
    std::vector<double> _RS;
    

public:
    MultigridSolver();
    /** 
     * the constructor of MultigridSolver
     * 
     * @param n the total numbers of level
     * @param f the vector which stores the rightside of each level
     * @param v the vector which stores the approximate solution of each level
     * @param S1 the keywords of restriction operator type, 
     *           "FullWeighting" or "Injection"
     * @param S2 the keywords of interpolation operator type,
     *           "Linear" or "Quadratic"
     * @param S3 the keywords of cycle type
     *           "VC" or "FMG"
     */
    MultigridSolver(int n, std::vector<double> f, std::vector<double> v, std::string S1 = "FullWeighting", std::string S2 = "Linear", std::string S3 = "VC");
    void SetGridLevel(int n);
    void SetRightSide(std::vector<double> f);
    void SetInitialGuess(std::vector<double> v);
    void SetBoundaryCond(double u0, double u1);
    void SetRestrictionType(std::string S);
    void SetNowLevel(int nowlevel);
    void SetRS(std::vector<double> RS);
    void PrintInfo();
    RestrictionOperator* pRestrictOP();
    InterpolationOperator* pInterpolateOP();
    void UpdateIndex();
    void WeightedJacobi();
    void BottomSolve(double u0, double u1);
    void VCycle(int StartLevel);
    void FMG();
    void Solve();
    std::vector<double> ReturnSolution();
    double RE_2Norm();
    

    
};


