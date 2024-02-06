#include "Marker.h"


Mark::Mark(string MN)
{
    this->E_name_ = MN;
    this->TypeID_ = 10;
    this->TM_.setZero();
    this->Lslice_ = 0.0;
}

Mark::Mark(const Mark& obj)
{
    this->E_name_ = obj.E_name_;
    this->TypeID_ = 10;
    this->TM_.setZero();
    this->Lslice_ = 0.0;
}

Mark::~Mark()
{
}

void Mark::integration(double& x, double& px, double& y, double& py, double& delta, double& tau)
{

}

void Mark::calc_TM()
{

    this->TM_.setZero();
    this->TM_(0, 0) = 1.0;
    this->TM_(1, 1) = 1.0;
    this->TM_(2, 2) = 1.0;
    this->TM_(3, 3) = 1.0;
    this->TM_(4, 4) = 1.0;
    this->TM_(5, 5) = 1.0;

    TMsome_.push_back(TM_);

}

