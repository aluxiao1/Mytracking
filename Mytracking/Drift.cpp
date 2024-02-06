#include "Drift.h"

Drift::Drift(string DN, double DL)
{
    this->D_step_=6;
    this->E_name_=DN;
    this->length_=DL;
    this->TypeID_=20;
    this->TM_.setZero();
    this->Lslice_ = this->length_ / this->D_step_;
}

Drift::~Drift()
{
}

Drift::Drift(const Drift& obj)
{
    this->D_step_ = 6;
    this->E_name_ = obj.E_name_;
    this->length_ = obj.length_;
    this->TypeID_ = 20;
    this->TM_.setZero();
    this->Lslice_ = this->length_ / this->D_step_;
}

void Drift::calc(double &x,double &px, double &y, double &py, double &delta,double &tau)
{

double pz;
double x1;
double x3;
double l0;
pz=sqrt((1.0+delta)*(1.0+delta)-px*px-py*py);
l0=this->length_/this->D_step_;
x1=x+l0*px/pz;
x3=y+l0*py/pz;
//tau=tau+l0*(1+delta)/pz-l0;
tau = tau + l0 * delta / pz + l0 * (1.0 - pz * pz) / pz / (1.0 + pz) ;

x=x1;
y=x3;
}

void Drift::integration(double &x,double &px, double &y, double &py, double &delta,double &tau)
{
for (int i = 0; i < this->D_step_; i++)
    {
    Drift::calc(x,px,y,py,delta,tau);
   }
}

void Drift::calc_TM()
{
    
    this->TM_.setZero();
    this->TM_(0,0)=1.0;
    this->TM_(1,1)=1.0;
    this->TM_(2,2)=1.0;
    this->TM_(3,3)=1.0;
    this->TM_(4,4)=1.0;
    this->TM_(5,5)=1.0;

    Matrix<double,6,6> M1;
    M1.setZero();

    double x=this->x0_;
    double px=this->px0_;
    double y=this->y0_;
    double py=this->py0_;
    double delta=this->delta0_;
    double tau=this->tau0_;
    double L=this->length_/this->D_step_;

    cout << "Drift calc_TM" << endl;

    for (int i = 0; i < D_step_; i++)
    {
        double pz=sqrt((1.0+delta)*(1.0+delta)-px*px-py*py);

        M1(0,0)=1.0;  
        M1(0,1)=L*(pz*pz+px*px)/(pow(pz, 3.0));
        M1(0,3)=L*px*py/pow(pz, 3.0);
        M1(0,4)=-L*px*(1+delta)/pow(pz, 3.0);

        M1(1,1)=1.0;

        M1(2,1)=L*px*py/pow(pz, 3.0);
        M1(2,2)=1.0;
        M1(2,3)=L*(pz*pz+py*py)/(pow(pz, 3.0));
        M1(2,4)=-L*py*(1+delta)/pow(pz, 3.0);

        M1(3,3)=1.0;

        M1(4,4)=1.0;

        M1(5, 1) = L* px * (1 + delta) / (pz*pz*pz);
        M1(5,3)= L * py*(1+delta)/ (pz * pz * pz);
        M1(5,4)=-L * (px * px + py * py) / (pz * pz * pz);
        M1(5,5)=1 ;

        this->TM_=M1*this->TM_;

        Drift::calc(x,px,y,py,delta,tau);

        TMsome_.push_back(TM_);
    }
    

}

