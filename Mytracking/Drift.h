#pragma once
#include "Element.h"

using namespace Eigen;


class Drift: public Element
{
private:
    /* data */
public:
    Drift(string DN, double DL);
    Drift(const Drift& obj);

    ~Drift();
    double D_step_;
    void calc(double &x,double &px, double &y, double &py, double &delta,double &tau);
    void integration(double &x,double &px, double &y, double &py, double &delta, double &tau);
    

    void calc_TM();

};
