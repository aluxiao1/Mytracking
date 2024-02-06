#pragma once
#include "Drift.h"
#include "Bend.h"
#include "Quad.h"
#include "Sext.h"
#include "Marker.h"
#include "MyMath.h"
#include <regex>
#include<fstream>

using namespace std;
using namespace Eigen;
#define Filename "Beamline.dat"
//#define Filename "DBA_QF.dat"
//#define Filename "HiSOR2.dat"
#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

class BeamLine
{
private:
    /* data */
public:
    BeamLine(/* args */);
    ~BeamLine();
    vector<shared_ptr<Element>> Element_list_;
    //vector<Element> 
    string BL_name_;
};

class Beamlines
{
public:
    vector<BeamLine> BL_list_;
    vector<Drift> Drift_list_;
    vector<Bend> Bend_list_;
    vector<Quad> Quad_list_;
    vector<Sext> Sext_list_;
    vector<Mark> Mark_list_;
    BeamLine Main_List_;

    Matrix<double,6,6> Mat_;


    void ReadFile();
    void show();
    void Calc_L();

    void Find_Orbit();
    void One_turn_Matrix();
    void Calc_Twiss_nocavi();
    void Calc_Twiss();
};
