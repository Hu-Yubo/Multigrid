/**
 * @file   Element.cpp
 * @author HirasawaYui <yui@Ubuntu18-04>
 * @date   Mon Oct 26 10:00:01 2020
 * 
 * @brief  
 * 
 * 
 */

#include "Element.h"

double Node::x()
{
    double h = 1.0 / _SdLen;
    int j = _GlbIdx % (_SdLen+1);
    return j*h;
}

double Node::y()
{
    double h = 1.0 / _SdLen;
    int i = _GlbIdx / (_SdLen+1);
    return i*h;
}

Element::Element(Node N1,  Node N2, Node N3, Node N4)
{
    _Node.push_back(N1);
    _Node.push_back(N2);
    _Node.push_back(N3);
    _Node.push_back(N4);
}

double Element::phi(double xi, double eta, int i)
{
    switch(i)
    {
    case 1:
	return (xi-1)*(eta-1)/4.0;
    case 2:
	return -(xi+1)*(eta-1)/4.0;
    case 3:
	return (xi+1)*(eta+1)/4.0;
    case 4:
	return -(xi-1)*(eta+1)/4.0;
    }
}

double Element::phi_xi(double xi, double eta, int i)
{
    switch(i)
    {
    case 1:
	return (eta-1)/4.0;
    case 2:
	return -(eta-1)/4.0;
    case 3:
	return (eta+1)/4.0;
    case 4:
	return -(eta+1)/4.0;
    }
}

double Element::phi_eta(double xi, double eta, int i)
{
    switch(i)
    {
    case 1:
	return (xi-1)/4.0;
    case 2:
	return -(xi+1)/4.0;
    case 3:
	return (xi+1)/4.0;
    case 4:
	return -(xi-1)/4.0;
    }
}

double Element::det_Jacobi(double xi, double eta)
{
    double J11 = 0;
    double J12 = 0;
    double J21 = 0;
    double J22 = 0;
    for (int i = 1; i <= 4; i++)
    {
	J11 = J11 + _Node[i-1].x()*phi_xi(xi, eta, i);
	J12 = J12 + _Node[i-1].y()*phi_xi(xi, eta, i);
	J21 = J21 + _Node[i-1].x()*phi_eta(xi, eta, i);
	J22 = J22 + _Node[i-1].y()*phi_eta(xi, eta, i);
    }
    return (J11*J22 - J12*J21);
}

double Element::xi_x(double xi, double eta)
{
    double a = 0;
    for (int i = 1; i <= 4; i++)
	a = a + _Node[i-1].y()*phi_eta(xi, eta, i);
    return a / det_Jacobi(xi, eta);
}

double Element::xi_y(double xi, double eta)
{
    double a = 0;
    for (int i = 1; i <= 4; i++)
	a = a + _Node[i-1].x()*phi_eta(xi, eta, i);
    return -a / det_Jacobi(xi, eta);
}

double Element::eta_x(double xi, double eta)
{
    double a = 0;
    for (int i = 1; i <= 4; i++)
	a = a + _Node[i-1].y()*phi_xi(xi, eta, i);
    return -a / det_Jacobi(xi, eta);
}

double Element::eta_y(double xi, double eta)
{
    double a = 0;
    for (int i = 1; i <= 4; i++)
	a = a + _Node[i-1].x()*phi_xi(xi, eta, i);
    return a / det_Jacobi(xi, eta);
}

double Element::phi_x(double xi, double eta, int i)
{
    return phi_xi(xi, eta, i)*xi_x(xi, eta) + phi_eta(xi, eta, i)*eta_x(xi, eta);
}

double Element::phi_y(double xi, double eta, int i)
{
    return phi_xi(xi, eta, i)*xi_y(xi, eta) + phi_eta(xi, eta, i)*eta_y(xi, eta);
}

double Element::a_ij(int i, int j)
{
    double a = 0;
    for (int k = 0; k < 4; k++)
    {
	double xi = _GaussPnt[k][0];
	double eta = _GaussPnt[k][1];
	a = a + det_Jacobi(xi, eta) * _w[k] * (phi_x(xi, eta, i)*phi_x(xi, eta, j)
					       + phi_y(xi, eta, i)*phi_y(xi, eta, j));
    }
    return a;
}

double Element::f(double xi, double eta)
{
    return 2*PI*PI*sin(PI*xi)*sin(PI*eta);
}

double Element::rhs(int i)
{
    double a = 0;
    for (int k = 0; k < 4; k++)
    {
	double xi = _GaussPnt[k][0];
	double eta = _GaussPnt[k][1];
	a = a + det_Jacobi(xi, eta) * _w[k] * f(xi, eta) * phi(xi, eta, i);
    }
    return a;
}

