/**
 * @file   Element.h
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Mon Oct 26 03:36:40 2020
 * 
 * @brief  Q1 element for solving 2D possion equation by Multigrid method
 * 
 * 
 */

#include <iostream>
#include <vector>
#include <math.h>

#define PI 4.0*atan(1.0)

class Node
{
private:
    int _GlbIdx;
    int _SdLen;
public:
    Node(int a, int n){
	_GlbIdx = a;
	_SdLen = n;
    }
    /// Suppose the domain of this problem is [0,1]x[0,1]
    double x();
    double y();
};


class Element
{
private:
    std::vector<Node> _Node;
    double _GaussPnt[4][2] = {{-1/sqrt(3), -1/sqrt(3)}, {1/sqrt(3), -1/sqrt(3)}, {1/sqrt(3), 1/sqrt(3)}, {-1/sqrt(3), 1/sqrt(3)}};
    double _w[4] = {1, 1, 1, 1};

public:
    Element(Node N1, Node N2, Node N3, Node N4);
    Element(std::vector<Node> NV){
	_Node = NV;
    };
    double phi(double xi, double eta, int i);
    double phi_xi(double xi, double eta, int i);
    double phi_eta(double xi, double eta, int i);
    double det_Jacobi(double xi, double eta);
    double xi_x(double xi, double eta);
    double xi_y(double xi, double eta);
    double eta_x(double xi, double eta);
    double eta_y(double xi, double eta);
    double phi_x(double xi, double eta, int i);
    double phi_y(double xi, double eta, int i);
    double a_ij(int i, int j);
    double f(double xi, double eta);
    double rhs(int i);
};
