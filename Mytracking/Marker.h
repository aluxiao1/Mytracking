#pragma once
#include "Element.h"
using namespace Eigen;

class Mark : public Element
{
private:
    /* data */
public:
    Mark(string MN);
    Mark(const Mark& obj);
    ~Mark();

    void integration(double& x, double& px, double& y, double& py, double& delta, double& tau);
    void calc_TM();

};



