#pragma once
#include "Element.h"

class Sext: public Element
{
private:
    /* data */
public:
    Sext(string SN, double SL, double SK);
    Sext(const Sext& obj);

    ~Sext();


    double K2_;
    double S_step_;


    void Calc_Drift(double &x,double &px, double &y, double &py, double &delta,double &tau,double L);
    void Calc_Kick(double &x,double &px, double &y, double &py, double &delta,double &tau,double L);
    void integration(double &x,double &px, double &y, double &py, double &delta,double &tau);
    
    Matrix <double, 6, 6> M1_;

    void calc_Drift_TM(double x, double px, double y, double py, double delta, double tau, double L);
    void calc_Kick_TM(double x, double px, double y, double py, double delta, double tau, double L);

    void calc_TM();

};
