#include "Quad.h"

Quad::Quad(string QN, double QL, double QK)
{
    
    this->Q_step_=6;
    this->E_name_=QN;
    this->length_=QL;
    this->K1_=QK;
    this->TypeID_=23;
    this->TM_.setZero();
    this->Lslice_ = this->length_ / this->Q_step_;
}


Quad::Quad(const Quad& obj)
{
    this->Q_step_ = 6;
    this->E_name_ = obj.E_name_;
    this->length_ = obj.length_;
    this->K1_ = obj.K1_;
    this->TypeID_ = 23;
    this->TM_.setZero();
    this->Lslice_ = this->length_ / this->Q_step_;
}

Quad::~Quad()
{
}

void Quad::Calc_Drift(double &x,double &px, double &y, double &py, double &delta,double &tau, double L)
{
double pz;
double x1;
double x3;
pz=sqrt((1+delta)*(1+delta)-px*px-py*py);

x1=x+L*px/pz;
x3=y+L*py/pz;
tau = tau + L * delta / pz + L * (1.0 - pz * pz) / pz / (1.0 + pz);

x=x1;
y=x3;

}

void Quad::Calc_Kick(double &x,double &px, double &y, double &py, double &delta,double &tau, double L)
{

px=px-this->K1_/this->length_*x*L;
py=py+this->K1_/this->length_*y*L;

}



void Quad::integration(double &x,double &px, double &y, double &py, double &delta,double &tau)
{
// a fourth order integrator

double LS1=this->length_/this->Q_step_/(2.0-cbrt(2.0));
double LS2=this->length_/this->Q_step_*(1.0-2.0/(2.0-cbrt(2.0)));

for (int i = 0; i < this->Q_step_; i++)
    {
    Quad::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
    Quad::Calc_Kick(x,px,y,py,delta,tau,LS1);
    Quad::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);

    Quad::Calc_Drift(x,px,y,py,delta,tau,LS2/2.0);
    Quad::Calc_Kick(x,px,y,py,delta,tau,LS2);
    Quad::Calc_Drift(x,px,y,py,delta,tau,LS2/2.0);

    Quad::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
    Quad::Calc_Kick(x,px,y,py,delta,tau,LS1);
    Quad::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
   }


}


void Quad::calc_Drift_TM(double x,double px, double y, double py, double delta,double tau,double L)
{
    
    this->M1_.setZero();
    double pz=sqrt((1.0+delta)*(1.0+delta)-px*px-py*py);

    this->M1_(0,0)=1.0;  
    this->M1_(0,1)=L*(pz*pz+px*px)/(pow(pz, 3.0));
    this->M1_(0,3)=L*px*py/pow(pz, 3.0);
    this->M1_(0,4)=-L*px*(1+delta)/(pz*pz*pz);

    this->M1_(1,1)=1.0;

    this->M1_(2,1)=L*px*py/pow(pz, 3.0);
    this->M1_(2,2)=1.0;
    this->M1_(2,3)=L*(pz*pz+py*py)/(pow(pz, 3.0));
    this->M1_(2,4)=-L*py*(1+delta)/ (pz * pz * pz);

    this->M1_(3,3)=1.0;

    this->M1_(4,4)=1.0;


    this->M1_(5,1)=L * px*(1+delta)/ (pz * pz * pz);
    this->M1_(5,3)=L * py*(1+delta)/ (pz * pz * pz);
    this->M1_(5,4)=-L * (px*px+py*py)/ (pz * pz * pz);
    this->M1_(5,5)=1;
    
}


void Quad::calc_Kick_TM(double x,double px, double y, double py, double delta,double tau,double L)
{   
    
    this->M1_.setZero();

    this->M1_(0,0)=1.0; 
    this->M1_(1,1)=1.0; 
    this->M1_(2,2)=1.0; 
    this->M1_(3,3)=1.0; 
    this->M1_(4,4)=1.0; 
    this->M1_(5,5)=1.0; 

    this->M1_(1,0)=-this->K1_/this->length_*L;
    this->M1_(3,2)=this->K1_/this->length_*L;

}



void Quad::calc_TM()
{


    double x=this->x0_;
    double px=this->px0_;
    double y=this->y0_;
    double py=this->py0_;
    double delta=this->delta0_;
    double tau=this->tau0_;   //get the fixed point


    double LS1=this->length_/this->Q_step_/(2.0-cbrt(2.0));
    double LS2=this->length_/this->Q_step_*(1.0-2.0/(2.0-cbrt(2.0)));
    this->TM_.setZero();
    this->TM_(0,0)=1.0;
    this->TM_(1,1)=1.0;
    this->TM_(2,2)=1.0;
    this->TM_(3,3)=1.0;
    this->TM_(4,4)=1.0;
    this->TM_(5,5)=1.0;

    for (int i = 0; i < this->Q_step_; i++)
    {
        //cout<<"initial TM:"<<endl<<this->TM_<<endl;

        //cout<<i<<endl;
        //cout<<x<<"  "<<px<<"  "<<y<<"  "<<py<<"  "<<delta<<"  "<<tau<<"  "<<endl;
        Quad::calc_Drift_TM(x, px, y, py, delta, tau, LS1/2.0);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
        Quad::calc_Kick_TM(x, px, y, py, delta, tau, LS1);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Kick(x,px,y,py,delta,tau,LS1);
        Quad::calc_Drift_TM(x, px, y, py, delta, tau, LS1/2.0);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);

        Quad::calc_Drift_TM(x, px, y, py, delta, tau, LS2/2.0);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Drift(x,px,y,py,delta,tau,LS2/2.0);
        Quad::calc_Kick_TM(x, px, y, py, delta, tau, LS2);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Kick(x,px,y,py,delta,tau,LS2);
        Quad::calc_Drift_TM(x, px, y, py, delta, tau, LS2/2.0);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Drift(x,px,y,py,delta,tau,LS2/2.0);


        Quad::calc_Drift_TM(x, px, y, py, delta, tau, LS1/2.0);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
        Quad::calc_Kick_TM(x, px, y, py, delta, tau, LS1);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Kick(x,px,y,py,delta,tau,LS1);
        Quad::calc_Drift_TM(x, px, y, py, delta, tau, LS1/2.0);
        this->TM_=this->M1_*this->TM_;
        Quad::Calc_Drift(x,px,y,py,delta,tau,LS1/2.0);
        
        TMsome_.push_back(TM_);

    }




}



