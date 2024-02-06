#pragma once
#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include<Eigen/Dense>
using namespace std;
using namespace Eigen;

class Element
{
private:
    /* data */
public:
    Element(/* args */);
    ~Element();


    int ID_;
    int TypeID_;
    string E_name_;

    double length_;
    double s_;   // position in a beamline
    void info();
    

    double x0_;    //fixed point information
    double px0_;
    double y0_;
    double py0_;
    double delta0_;
    double tau0_;

    virtual void integration(double &x,double &px, double &y, double &py, double &delta,double &tau);
    
    Eigen::Matrix<double, 6, 6> TM_;     // total TM
    vector< Eigen::Matrix <double, 6, 6> > TMsome_;  // TM of each slices
    vector< Matrix <double, 6, 6> > getdata() { return TMsome_; }
    double Lslice_;
    virtual void calc_TM();

    double AX_;
    double BX_;
    double EX_;
    double NX_;
    double AY_;
    double BY_;
    double EY_;
    double NY_;

    
};


