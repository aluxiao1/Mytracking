#include "Bend.h"

Bend::Bend(string Bname, double BL, double theta, double K0, double K1, double E1, double E2)
{
    this->B_step_=10;
    this->E_name_=Bname;
    this->length_=BL;
    this->ANGLE_=theta;
    this->K0_=K0;
    this->K1_=K1;
    this->E1_=E1;
    this->E2_=E2;

    if (theta==0)
    {this->TypeID_=22;}   // a kicker
    else
    {this->TypeID_=21;}  // bending mag
    this->Lslice_ = this->length_ / this->B_step_;
}

Bend::Bend(const Bend& obj)
{
    this->B_step_ = 10;
    this->E_name_ = obj.E_name_;
    this->length_ = obj.length_;
    this->ANGLE_ = obj.ANGLE_;
    this->K0_ = obj.K0_;
    this->K1_ = obj.K1_;
    this->E1_ = obj.E1_;
    this->E2_ = obj.E2_;

    if (this->ANGLE_ == 0)
    {
        this->TypeID_ = 22;
    }   // a kicker
    else
    {
        this->TypeID_ = 21;
    }  // bending mag
    this->Lslice_ = this->length_ / this->B_step_;
}

Bend::~Bend()
{
}


void Bend::Calc_K_Drift(double &x,double &px, double &y, double &py, double &delta,double &tau,double L)
{
double pz;
double x1;
double x3;
pz=sqrt((1+delta)*(1+delta)-px*px-py*py);
x1=x+L*px/pz;
x3=y+L*py/pz;
tau=tau+L*(1.0+delta)/pz-L;
x=x1;
y=x3;
}

void Bend::Calc_K_Kick(double &x,double &px, double &y, double &py, double &delta,double &tau,double L)
{
if (L!=0)
{
   px=px-this->K0_/this->length_*L;
}
else
{
    px=px-this->K0_; 
}

}

void Bend::Calc_S_Drift(double &x,double &px, double &y, double &py, double &delta,double &tau,double L, double rho)
{

double A=L/rho;
double pz=sqrt((1.0+delta)*(1.0+delta)-px*px-py*py);

double x1;
double x2;
double x3;
double x6;

x1=(x+rho*(2*sin(A/2)*sin(A/2)+px*sin(A)/pz))/cos(A)/(1-px*tan(A)/pz);
x2=px*cos(A)+sin(A)*pz;
x3=y+py*(x+rho)*tan(A)/pz/(1-px*tan(A)/pz);
x6=tau+(1+delta)*(x+rho)*tan(A)/pz/(1-px*tan(A)/pz)-L;

x=x1;
px=x2;
y=x3;
tau=x6;

}

void Bend::Calc_S_Kick(double &x,double &px, double &y, double &py, double &delta,double &tau,double L, double rho)
{

double b=this->ANGLE_/this->length_;
px=px-(b+b*x/rho)*L-(1+x/rho)*L*(this->K1_/this->length_*x);
py=py+(1+x/rho)*L*(this->K1_/this->length_*y);


}



void Bend::integration(double &x,double &px, double &y, double &py, double &delta,double &tau)
{

    
    if (this->TypeID_==21)
    {
    /* Bending magnet */
    double r = this->length_/this->ANGLE_;
    double LS1=this->length_/this->B_step_/(2-cbrt(2));
    double LS2=this->length_/this->B_step_*(1-2/(2-cbrt(2)));

        for (int i = 0; i < this->B_step_; i++)
            {
            Bend::Calc_S_Drift(x,px,y,py,delta,tau,LS1/2,r);
            Bend::Calc_S_Kick(x,px,y,py,delta,tau,LS1,r);
            Bend::Calc_S_Drift(x,px,y,py,delta,tau,LS1/2,r);

            Bend::Calc_S_Drift(x,px,y,py,delta,tau,LS2/2,r);
            Bend::Calc_S_Kick(x,px,y,py,delta,tau,LS2,r);
            Bend::Calc_S_Drift(x,px,y,py,delta,tau,LS2/2,r);

            Bend::Calc_S_Drift(x,px,y,py,delta,tau,LS1/2,r);
            Bend::Calc_S_Kick(x,px,y,py,delta,tau,LS1,r);
            Bend::Calc_S_Drift(x,px,y,py,delta,tau,LS1/2,r);
            }
    }
    else if(this->TypeID_==22 && this->length_!=0)
    {   
        /* Kicker magnet */
        double LS1=this->length_/this->B_step_/(2-cbrt(2));
        double LS2=this->length_/this->B_step_*(1-2/(2-cbrt(2)));

        for (int i = 0; i < this->B_step_; i++)
            {
            Bend::Calc_K_Drift(x,px,y,py,delta,tau,LS1/2);
            Bend::Calc_K_Kick(x,px,y,py,delta,tau,LS1);
            Bend::Calc_K_Drift(x,px,y,py,delta,tau,LS1/2);

            Bend::Calc_K_Drift(x,px,y,py,delta,tau,LS2/2);
            Bend::Calc_K_Kick(x,px,y,py,delta,tau,LS2);
            Bend::Calc_K_Drift(x,px,y,py,delta,tau,LS2/2);

            Bend::Calc_K_Drift(x,px,y,py,delta,tau,LS1/2);
            Bend::Calc_K_Kick(x,px,y,py,delta,tau,LS1);
            Bend::Calc_K_Drift(x,px,y,py,delta,tau,LS1/2);
        }
    }
    else if(this->TypeID_==22 && this->length_==0)
    {   
        /* Kicker magnet---thin lens kick*/
        Bend::Calc_K_Kick(x,px,y,py,delta,tau,0);
    }

    


}


void Bend::calc_SDrift_TM(double x, double px, double y, double py, double delta, double tau, double L, double rho)
{
    double A = L / rho;
    this->M1_.setZero();
    double pz = sqrt((1.0 + delta) * (1.0 + delta) - px * px - py * py);

    this->M1_(0, 0) = 1 / cos(A) / (1 - px * tan(A) / pz);
    this->M1_(0, 1) = (rho * sin(A) + (x + 2 * rho * sin(A / 2) * sin(A / 2)) * tan(A)) * (pz + px * px / pz) / cos(A) / (pz - px * tan(A)) / (pz - px * tan(A));
    this->M1_(0, 3) = px * py / pz * (rho * sin(A) + (x + 2 * rho * sin(A / 2) * sin(A / 2)) * tan(A)) / cos(A) / (pz - px * tan(A)) / (pz - px * tan(A));
    double T1;
    T1 = ((1 + delta) / pz * x + 2 * rho * (1 + delta) / pz * sin(A / 2) * sin(A / 2)) * (pz - px * tan(A)) - (pz * x + rho * (2 * pz * sin(A / 2) * sin(A / 2) + px * sin(A))) * (1 + delta) / pz;
    this->M1_(0, 4) = T1 / cos(A) / (pz - px * tan(A)) / (pz - px * tan(A));


    this->M1_(1, 1) = cos(A) - sin(A) * px / pz;
    this->M1_(1, 3) = -py / pz;
    this->M1_(1, 4) = sin(A) * (1 + delta) / pz;

    this->M1_(2, 0) = py * tan(A) / (pz - px * tan(A));
    this->M1_(2, 1) = (x + rho) * py * tan(A) * (px / pz + tan(A)) / (pz - px * tan(A)) / (pz - px * tan(A));
    this->M1_(2, 2) = 1.0;
    this->M1_(2, 3) = (x + rho) * tan(A) / (pz - px * tan(A)) + py * (x + rho) * tan(A) * (py / pz) / (pz - px * tan(A)) / (pz - px * tan(A));
    this->M1_(2, 4) = -py * (x + rho) * tan(A) * (1 + delta) / pz / (pz - px * tan(A)) / (pz - px * tan(A));


    this->M1_(3, 3) = 1.0;
    this->M1_(4, 4) = 1.0;

    this->M1_(5, 0) = (1 + delta) * tan(A) / (pz - px * tan(A));
    this->M1_(5, 1) = (1 + delta) * (rho + x) * tan(A) * (px / pz + tan(A)) / (pz - px * tan(A)) / (pz - px * tan(A));
    this->M1_(5, 3) = (1 + delta) * (rho + x) * tan(A) * (py / pz) / (pz - px * tan(A)) / (pz - px * tan(A));
    this->M1_(5, 4) = (x + rho) * tan(A) / (pz - px * tan(A)) - (1 + delta) * (1 + delta) * (x + rho) * tan(A) / pz / (pz - px * tan(A)) / (pz - px * tan(A));
    this->M1_(5, 5) = 1.0;



}

void Bend::calc_SKick_TM(double x, double px, double y, double py, double delta, double tau, double L, double rho)
{

    double b = this->ANGLE_ / this->length_;

    this->M1_.setZero();

    this->M1_(0, 0) = 1.0;
    this->M1_(1, 1) = 1.0;
    this->M1_(2, 2) = 1.0;
    this->M1_(3, 3) = 1.0;
    this->M1_(4, 4) = 1.0;
    this->M1_(5, 5) = 1.0;

    this->M1_(1, 0) = -b / rho * L - this->K1_ * L / this->length_ * (1 + 2 * x / rho);

    this->M1_(3, 2) = (1 + x / rho) * L * this->K1_ / this->length_;

}

void Bend::calc_KDrift_TM(double x, double px, double y, double py, double delta, double tau, double L)
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
    this->M1_(5, 5) = 1.0;
}

void Bend::calc_KKick_TM(double x, double px, double y, double py, double delta, double tau, double L)
{
    this->M1_.setZero();

    this->M1_(0, 0) = 1.0;
    this->M1_(1, 1) = 1.0;
    this->M1_(2, 2) = 1.0;
    this->M1_(3, 3) = 1.0;
    this->M1_(4, 4) = 1.0;
    this->M1_(5, 5) = 1.0;

}





void Bend::calc_TM()
{
    double x = this->x0_;
    double px = this->px0_;
    double y = this->y0_;
    double py = this->py0_;
    double delta = this->delta0_;
    double tau = this->tau0_;   //get the fixed point

    this->TM_.setZero();
    this->TM_(0, 0) = 1.0;
    this->TM_(1, 1) = 1.0;
    this->TM_(2, 2) = 1.0;
    this->TM_(3, 3) = 1.0;
    this->TM_(4, 4) = 1.0;
    this->TM_(5, 5) = 1.0;

    if (this->TypeID_ == 21)
    {
        double r = this->length_ / this->ANGLE_;
        double LS1 = this->length_ / this->B_step_ / (2 - cbrt(2));
        double LS2 = this->length_ / this->B_step_ * (1 - 2 / (2 - cbrt(2)));


        for (int i = 0; i < this->B_step_; i++)
        {
            Bend::calc_SDrift_TM(x, px, y, py, delta, tau, LS1 / 2.0, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Drift(x, px, y, py, delta, tau, LS1 / 2.0, r);
            Bend::calc_SKick_TM(x, px, y, py, delta, tau, LS1, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Kick(x, px, y, py, delta, tau, LS1, r);
            Bend::calc_SDrift_TM(x, px, y, py, delta, tau, LS1 / 2.0, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Drift(x, px, y, py, delta, tau, LS1 / 2.0, r);

            Bend::calc_SDrift_TM(x, px, y, py, delta, tau, LS2 / 2.0, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Drift(x, px, y, py, delta, tau, LS2 / 2.0, r);
            Bend::calc_SKick_TM(x, px, y, py, delta, tau, LS2, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Kick(x, px, y, py, delta, tau, LS2, r);
            Bend::calc_SDrift_TM(x, px, y, py, delta, tau, LS2 / 2.0, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Drift(x, px, y, py, delta, tau, LS2 / 2.0, r);

            Bend::calc_SDrift_TM(x, px, y, py, delta, tau, LS1 / 2.0, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Drift(x, px, y, py, delta, tau, LS1 / 2.0, r);
            Bend::calc_SKick_TM(x, px, y, py, delta, tau, LS1, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Kick(x, px, y, py, delta, tau, LS1, r);
            Bend::calc_SDrift_TM(x, px, y, py, delta, tau, LS1 / 2.0, r);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_S_Drift(x, px, y, py, delta, tau, LS1 / 2.0, r);

            TMsome_.push_back(TM_);
        }
    }
    else if (this->TypeID_ == 21 && this->length_ != 0)
    {
        double LS1 = this->length_ / this->B_step_ / (2 - cbrt(2));
        double LS2 = this->length_ / this->B_step_ * (1 - 2 / (2 - cbrt(2)));

        for (int i = 0; i < this->B_step_; i++)
        {
            Bend::calc_KDrift_TM(x, px, y, py, delta, tau, LS1 / 2.0);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Drift(x, px, y, py, delta, tau, LS1 / 2.0);
            Bend::calc_KKick_TM(x, px, y, py, delta, tau, LS1);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Kick(x, px, y, py, delta, tau, LS1);
            Bend::calc_KDrift_TM(x, px, y, py, delta, tau, LS1 / 2.0);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Drift(x, px, y, py, delta, tau, LS1 / 2.0);

            Bend::calc_KDrift_TM(x, px, y, py, delta, tau, LS2 / 2.0);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Drift(x, px, y, py, delta, tau, LS2 / 2.0);
            Bend::calc_KKick_TM(x, px, y, py, delta, tau, LS2);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Kick(x, px, y, py, delta, tau, LS2);
            Bend::calc_KDrift_TM(x, px, y, py, delta, tau, LS2 / 2.0);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Drift(x, px, y, py, delta, tau, LS2 / 2.0);

            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Drift(x, px, y, py, delta, tau, LS1 / 2.0);
            Bend::calc_KKick_TM(x, px, y, py, delta, tau, LS1);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Kick(x, px, y, py, delta, tau, LS1);
            Bend::calc_KDrift_TM(x, px, y, py, delta, tau, LS1 / 2.0);
            this->TM_ = this->M1_ * this->TM_;
            Bend::Calc_K_Drift(x, px, y, py, delta, tau, LS1 / 2.0);

            TMsome_.push_back(TM_);
        }
    }
    else if (this->TypeID_ == 22 && this->length_ == 0)
    {
        //// none , Identity matrix
    }



}