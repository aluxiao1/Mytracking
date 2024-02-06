#include "Element.h"


Element::Element(/* args */)
{
    this->length_=0.0;
    this->s_=0.0;

    this->x0_=0.0;    
    this->px0_=0.0;
    this->y0_= 0.0;
    this->py0_= 0.0;
    this->delta0_=0.0;
    this->tau0_=0.0;

    this->TM_.setZero();

    this->AX_ = 0.0;
    this->BX_ = 0.0;
    this->EX_ = 0.0;
    this->NX_ = 0.0;
    this->AY_ = 0.0;
    this->BY_ = 0.0;
    this->EY_ = 0.0;
    this->NY_ = 0.0;



}

Element::~Element()
{
}

void Element::info()
{
    
}

void Element::integration(double &x,double &px, double &y, double &py, double &delta,double &tau)
{

}

void Element::calc_TM()
{
    cout << "element calc_TM" << endl;
}