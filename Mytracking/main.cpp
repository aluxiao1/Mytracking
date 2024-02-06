#include "Beamline.h"
#include "MyMath.h"
//#include "matplotlibcpp.h"
#include<Eigen/Dense>
using namespace Eigen;
using namespace std;
//namespace plt = matplotlibcpp;


int main()
{

    Beamlines BL;
    BL.ReadFile();
    BL.show();

    //BL.Calc_L();
    //BL.Find_Orbit();
    //BL.One_turn_Matrix();

    //BL.Calc_Twiss_nocavi();  //BL.Calc_Twiss();

    //MyMath::Check_Sympletic(BL.Mat_);
 //   for (int i = 0; i < BL.Main_List_.Element_list_.size(); i++)   //BL.Main_List_.Element_list_.size()
 //   {
 //       cout << BL.Main_List_.Element_list_[i]->E_name_ << endl;
 //       //cout << BL.Main_List_.Element_list_[i]->TM_ << endl;
 //       BL.Main_List_.Element_list_[i]->calc_TM();
 //   }

     
 //double x1= 0.05;
 //double x2= 0.05;
 //double x3= 0.05;
 //double x4= 0.05;
 //double x5=0;
 //double x6=0;
 //double LC=0;
 //
 //
 //    for (int i = 0; i < BL.Main_List_.Element_list_.size(); i++)   //BL.Main_List_.Element_list_.size()
 //    {
 //        BL.Main_List_.Element_list_[i]->integration(x1,x2,x3,x4,x5,x6);
 //        LC=LC+BL.Main_List_.Element_list_[i]->length_;
 //        cout << BL.Main_List_.Element_list_[i]->E_name_ << endl;
 //        cout << x1 << endl;
 //        cout << x2 << endl;
 //        cout << x3 << endl;
 //        cout << x4 << endl;
 //        cout << x5 << endl;
 //        cout << x6 << endl;
 //     }
}