/**
 * @file   MultigridSolver.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Mon Oct 19 23:35:37 2020
 * 
 * @brief  solve 1D possion equation by multigrid method on domain [0,1]x[0,1].
 * 
 * 
 */

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
    /// the collection of right side
    std::vector<std::vector<double> > _f;
    /// the collection of approximate solution
    std::vector<std::vector<double> > _v;
    /// the Dirichlet condition
    /// while the nodes at corner belong to up and down boundary
    std::vector<double> _u_up;
    std::vector<double> _u_down;
    std::vector<double> _u_right;
    std::vector<double> _u_left;
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

public:
    MultigridSolver();
    MultigridSolver(int n, std::vector<double> f, std::vector<double> v, std::vector<double> u_up, std::vector<double> u_down, std::vector<double> u_left, std::vector<double> u_right, std::string S1 = "FullWeighting", std::string S2 = "Linear", std::string S3 = "VC");
    int IsBoundary(int nowlevel);
    std::vector<int> FindNeighbor(int i);
}
