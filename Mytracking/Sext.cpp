#include "Sext.h"

Sext::Sext(string SN, double SL, double SK)
{
    
    this->S_step_=10.0;
    this->E_name_=SN;
    this->length_=SL;
    this->K2_=SK;
    this->TypeID_=24;
    this->TM_.setZero();
    this->Lslice_ = this->length_ / this->S_step_;
}

Sext::Sext(const Sext& obj)
{
    this->S_step_ = 10.0;
    this->E_name_ = obj.E_name_;
    this->length_ = obj.length_;
    this->K2_ = obj.K2_;
    this->TypeID_ = 24;
    this->TM_.setZero();
    this->Lslice_ = this->length_ / this->S_step_;
}

Sext::~Sext()
{
}


void Sext::Calc_Drift(double &x,double &px, double &y, double &py, double &delta,double &tau, double L)
{
double pz;
double x1;
double x3;
pz=sqrt((1+delta)*(1+delta)-px*px-py*py);

x1=x+L*px/pz;
x3=y+L*py/pz;
//tau=tau+L*(1+delta)/pz-L;
tau=tau+L;
x=x1;
y=x3;
}


void Sext::Calc_Kick(double &x,double &px, double &y, double &py, double &delta,double &tau, double L)
{

if (L!=0)
{
   px=px-this->K2_/this->length_*(x*x-y*y)/2.0*L;
   py=py+this->K2_/this->length_*x*y*L;
}
else
{
   px=px-this->K2_*(x*x-y*y)/2.0;
   py=py+this->K2_*x*y;
}

}


void Sext::integration(double &x,double &px, double &y, double &py, double &delta,double &tau)
{
// a fourth order integrator

double LS1=this->length_/this->S_step_/(2.0-cbrt(2.0));
double LS2=this->length_/this->S_step_*(1.0-2.0/(2.0-cbrt(2.0)));

for (int i = 0; i < this->S_step_; i++)
    {
    Sext::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
    Sext::Calc_Kick(x,px,y,py,delta,tau,LS1);
    Sext::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);

    Sext::Calc_Drift(x,px,y,py,delta,tau,LS2/2.0);
    Sext::Calc_Kick(x,px,y,py,delta,tau,LS2);
    Sext::Calc_Drift(x,px,y,py,delta,tau,LS2/2.0);

    Sext::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
    Sext::Calc_Kick(x,px,y,py,delta,tau,LS1);
    Sext::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
   }

}


void Sext::calc_Drift_TM(double x, double px, double y, double py, double delta, double tau, double L)
{

    this->M1_.setZero();
    double pz = sqrt((1.0 + delta) * (1.0 + delta) - px * px - py * py);

    this->M1_(0, 0) = 1.0;
    this->M1_(0, 1) = L * (pz * pz + px * px) / (pow(pz, 3.0));
    this->M1_(0, 3) = L * px * py / pow(pz, 3.0);
    this->M1_(0, 4) = -L * px * (1 + delta) / pow(pz, 3.0);

    this->M1_(1, 1) = 1.0;

    this->M1_(2, 1) = L * px * py / pow(pz, 3.0);
    this->M1_(2, 2) = 1.0;
    this->M1_(2, 3) = L * (pz * pz + py * py) / (pow(pz, 3.0));
    this->M1_(2, 4) = -L * py * (1 + delta) / pow(pz, 3.0);

    this->M1_(3, 3) = 1.0;

    this->M1_(4, 4) = 1.0;


    this->M1_(5, 1) = L * px * (1 + delta) / pow(pz, 3.0);
    this->M1_(5, 3) = L * py * (1 + delta) / pow(pz, 3.0);
    this->M1_(5, 4) = L * (px * px + py * py) / pow(pz, 3.0);
    this->M1_(5, 5) = 1;

}


void Sext::calc_Kick_TM(double x, double px, double y, double py, double delta, double tau, double L)
{

    this->M1_.setZero();

    this->M1_(0, 0) = 1.0;
    this->M1_(1, 1) = 1.0;
    this->M1_(2, 2) = 1.0;
    this->M1_(3, 3) = 1.0;
    this->M1_(4, 4) = 1.0;
    this->M1_(5, 5) = 1.0;

    this->M1_(1, 0) = -this->K2_*x / this->length_ * L;
    this->M1_(1, 2) = this->K2_ * y / this->length_ * L;

    this->M1_(3, 0) = this->K2_*y / this->length_ * L;
    this->M1_(3, 2) = this->K2_*x / this->length_ * L;

}


void Sext::calc_TM()
{


    double x = this->x0_;
    double px = this->px0_;
    double y = this->y0_;
    double py = this->py0_;
    double delta = this->delta0_;
    double tau = this->tau0_;   //get the fixed point


    double LS1 = this->length_ / this->S_step_ / (2.0 - cbrt(2.0));
    double LS2 = this->length_ / this->S_step_ * (1.0 - 2.0 / (2.0 - cbrt(2.0)));
    this->TM_.setZero();
    this->TM_(0, 0) = 1.0;
    this->TM_(1, 1) = 1.0;
    this->TM_(2, 2) = 1.0;
    this->TM_(3, 3) = 1.0;
    this->TM_(4, 4) = 1.0;
    this->TM_(5, 5) = 1.0;



    for (int i = 0; i < this->S_step_; i++)
    {
        //cout<<"initial TM:"<<endl<<this->TM_<<endl;

        //cout<<i<<endl;
        //cout<<x<<"  "<<px<<"  "<<y<<"  "<<py<<"  "<<delta<<"  "<<tau<<"  "<<endl;
        Sext::calc_Drift_TM(x, px, y, py, delta, tau, LS1 / 2.0);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Drift(x, px, y, py, delta, tau, LS1 / 2.0);
        Sext::calc_Kick_TM(x, px, y, py, delta, tau, LS1);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Kick(x, px, y, py, delta, tau, LS1);
        Sext::calc_Drift_TM(x, px, y, py, delta, tau, LS1 / 2.0);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Drift(x, px, y, py, delta, tau, LS1 / 2.0);

        Sext::calc_Drift_TM(x, px, y, py, delta, tau, LS2 / 2.0);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Drift(x, px, y, py, delta, tau, LS2 / 2.0);
        Sext::calc_Kick_TM(x, px, y, py, delta, tau, LS2);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Kick(x, px, y, py, delta, tau, LS2);
        Sext::calc_Drift_TM(x, px, y, py, delta, tau, LS2 / 2.0);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Drift(x, px, y, py, delta, tau, LS2 / 2.0);


        Sext::calc_Drift_TM(x, px, y, py, delta, tau, LS1 / 2.0);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Drift(x, px, y, py, delta, tau, LS1 / 2.0);
        Sext::calc_Kick_TM(x, px, y, py, delta, tau, LS1);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Kick(x, px, y, py, delta, tau, LS1);
        Sext::calc_Drift_TM(x, px, y, py, delta, tau, LS1 / 2.0);
        this->TM_ = this->M1_ * this->TM_;
        Sext::Calc_Drift(x, px, y, py, delta, tau, LS1 / 2.0);


    }




}


