# pragma once
#include <iostream>
#include <string>
#include <stack>
#include <cmath>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

class Token 
{
public:
    char op;
    double num;
};

class MyMath
{
public:
double sstod(string expr);    // change operations to double   Shunting Yard Algorithm
static Matrix<double, 6, 6> Initialize_Ai_nocavi(Matrix<double, 6, 6> M);
static Matrix<double, 6, 6> Canonized_A_nocavi(Matrix<double, 6, 6> Ak);
static void Check_Sympletic(Matrix<double, 6, 6> M);

};

