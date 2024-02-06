#pragma once
#include "Element.h"
#include "Config.h"
using namespace Eigen;
class Bend:public Element
{
private:
    /* data */
public:
    Bend(string Bname, double BL, double theta, double K0, double K1, double E1, double E2);
    Bend(const Bend& obj);

    ~Bend();    
    double B_step_;
    double ANGLE_;
    double K0_;
    double K1_;
    double E1_;
    double E2_;
    Matrix <double, 6, 6> M1_;
    void calc_TM();


    void integration(double &x,double &px, double &y, double &py, double &delta,double &tau);

    void Calc_K_Drift(double &x,double &px, double &y, double &py, double &delta,double &tau,double L);
    void Calc_K_Kick(double &x,double &px, double &y, double &py, double &delta,double &tau,double L);
    void Calc_S_Drift(double &x,double &px, double &y, double &py, double &delta,double &tau,double L, double rho);
    void Calc_S_Kick(double &x,double &px, double &y, double &py, double &delta,double &tau,double L, double rho);

    void calc_SDrift_TM(double x, double px, double y, double py, double delta, double tau, double L, double rho);
    void calc_SKick_TM(double x, double px, double y, double py, double delta, double tau, double L, double rho);
    void calc_KDrift_TM(double x, double px, double y, double py, double delta, double tau, double L);
    void calc_KKick_TM(double x, double px, double y, double py, double delta, double tau, double L);

};
