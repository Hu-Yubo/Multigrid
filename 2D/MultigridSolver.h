/**
 * @file   MultigridSolver.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Mon Oct 19 23:35:37 2020
 * 
 * @brief  solve 2D possion equation by multigrid method on domain [0,1]x[0,1].
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

#define PI (4.0*atan(1.0))

class MultigridSolver
{
private:
    /// the total level of multigrid method
    int _n;
    /// the length of each unit
    double _h;
    /// the unit count of each side
    int _SdLen;
    /// the collection of right side
    std::vector<std::vector<double> > _f;
    /// the collection of approximate solution
    std::vector<std::vector<double> > _v;
    /// the level where MG is proceeding now
    int _nowlevel = 1;
    /// the pointer of restriction operator
    RestrictionOperator* _pRestrictOP;
    /// the pointer of interation operator
    InterpolationOperator* _pInterpolateOP;
    /// the weight of WeightedJacobi
    double _w = 2.0/3;
    /// the times of Relaxtion
    int _RlxTimes = 10;
    /// the keywords of cycle type
    std::string _TypeofCycle;
    /// the mark whether this node is boundary
    std::vector<int> _BndMark;
    /// the bottom boundary condition
    std::vector<double> _BtmBnd;
    /// the real solution
    std::vector<double> _RS;

public:
    MultigridSolver();
    MultigridSolver(int n, std::vector<double> f, std::vector<double> v, std::string S1 = "FullWeighting", std::string S2 = "Linear", std::string S3 = "VC");
    std::vector<int> IsBoundary();
    void SetRS(std::vector<double> RS);
    void UpdateData();
    int CortoIdx(int i, int j);
    void WeightedJacobi();
    void BottomSolve(std::vector<double> BtmBnd);
    void VCycle(int StartLevel);
    void Solve();
    void FMG();
    std::vector<double> ReturnSolution();
    double RE_2norm();
};
