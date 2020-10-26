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

class Node
{
private:
    int _GlbIdx;
    int _LocIdx;
    int _SdLen;
public:
    Node(int a, int b, int n){
	_GlbIdx = a;
	_LocIdx = b;
	_SdLen = n;
    }
    /// Suppose the domain of this problem is [0,1]x[0,1]
    double x();
    double y();
};


class Element
{
private:
    std::vector<Node> _Node(4);

public:
    Element(Node N1, Node N2, Node N3, Node N4);
    Element(std::vector<Node> NV){
	_Node = NV;
    };
    double phi(double xi, double eta, int i);
    double phi_xi(double xi, double eta, int i);
};
