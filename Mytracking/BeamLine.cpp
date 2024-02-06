#include "BeamLine.h"

BeamLine::BeamLine(/* args */)
{

}

BeamLine::~BeamLine()
{
}


void Beamlines::ReadFile()
{
    MyMath m1;
    ifstream f1;
    f1.open(Filename,ios::in);
    if(!f1.is_open())
    {
        cout<<"File does not exist!"<<endl;
        f1.close();
        return;
    }

    string line;
    string result;  // 储存删除"!"的字符串
    vector<string> lines; // 用容器用于储存分割";"后的文本

    while (getline(f1, line))  
    {  
        int pos = line.find('!'); // 查找"！"的位置
        if (pos != string::npos)   // 没找到就返回npos给 pos
        {        
            line.erase(pos); // 移除"！"及其之后的内容
        }
        if(!line.empty())
        {
            result += line;
        }
    }
    
    f1.close();
//    cout<<result<<endl;
    istringstream newline(result);
    string Eline;
    while (getline(newline, Eline, ';')) 
    {
        Eline.push_back('\r');
        lines.push_back(Eline);
    }
    
    regex pattern("([A-Z0-9]+)\\s*=\\s*\\(([^)]+)\\)"); // this pattern is for all element except for marker and monitor
    regex pattern1("([A-Z0-9]+)\\s*=\\s*\\(.*?\\)");     // .* any string even null, but greedy match ,   .*?  non-greedy match    
    smatch match;
    string p1;
    string p2;
    double vp2;

    for (int i = 0; i < lines.size(); i++)
    {
        
        if(lines[i].find("BEND") != string::npos)
        {
            size_t found = lines[i].find("BEND");
            lines[i].erase(found,4);
            //cout<<lines[i]<<endl;

            // 正则表达式模式用于匹配变量和参数
            string::const_iterator searchStart(lines[i].cbegin());

            while (regex_search(searchStart, lines[i].cend(), match, pattern)) 
            {
                string BName = match[1];
                string parametersString = match[2];
                parametersString.erase(remove_if(parametersString.begin(), parametersString.end(), ::isspace), parametersString.end());
                //cout<<BName<<endl;
                //cout<<parametersString<<endl;

                stringstream sst(parametersString);      
                string Bpara;  

                map<string,double> Bmap;

                while (getline(sst, Bpara, ','))      //逐个提取字符串中的单词（以空格作为分隔符）
                {
                    Bpara.erase(Bpara.find_last_not_of(" \n\r\t") + 1);
                    auto pos = Bpara.find("=");
                    p1=Bpara.substr(0,pos);
                    p2=Bpara.substr(pos+1);
                    vp2=m1.sstod(p2);
                    Bmap[p1] = vp2;
                }

                double B_L=0;
                double ANGLE=0;
                double B_K0=0;
                double B_K1=0;
                double B_E1=0;
                double B_E2=0;

                for(auto it = Bmap.cbegin(); it != Bmap.cend(); ++it)
                {
                    if(it->first=="L")
                    {
                        B_L=it->second;
                    }
                    else if(it->first=="ANGLE")
                    {
                        ANGLE=it->second;
                    }
                    else if(it->first=="K0")
                    {
                        B_K0=it->second;
                    }
                    else if(it->first=="K1")
                    {
                        B_K1=it->second;
                    }
                    else if(it->first=="E1")
                    {
                        B_E1=it->second;
                    }
                    else if(it->first=="E2")
                    {
                        B_E2=it->second;
                    }
                 //cout << it->first << " " << it->second << endl;
                }
                Bend EB=Bend(BName,B_L,ANGLE,B_K0,B_K1,B_E1,B_E2);
                this->Bend_list_.push_back(EB);
                // for (int j = 0; j < Bparas.size(); j++)
                // {
                //     cout<<Bparas[j]<<endl;
                // }
                searchStart = match.suffix().first;   //来获取从上一个匹配项结束位置到整个字符串末尾的子字符,下一个 std::regex_search 将从这里开始继续查找。
            }
        }
        else if(lines[i].find("QUAD") != string::npos)
        {
            size_t found = lines[i].find("QUAD");
            lines[i].erase(found,4);
            //cout<<lines[i]<<endl;

            string::const_iterator searchStart(lines[i].cbegin());
            while (regex_search(searchStart, lines[i].cend(), match, pattern)) 
            {
                string QName = match[1];
                string QString = match[2];
                QString.erase(remove_if(QString.begin(), QString.end(), ::isspace), QString.end());
                //cout<<QName<<endl;
                //cout<<QString<<endl;

                stringstream sst(QString);      
                string Qpara;  

                map<string,double> Qmap;

                while (getline(sst, Qpara, ','))      //逐个提取字符串中的单词（以空格作为分隔符）
                {
                    Qpara.erase(Qpara.find_last_not_of(" \n\r\t") + 1);
                    auto pos = Qpara.find("=");
                    p1=Qpara.substr(0,pos);
                    p2=Qpara.substr(pos+1);
                    vp2=m1.sstod(p2);
                    Qmap[p1] = vp2;

                }

                double Q_L=0;
                double Q_K1;
                for(auto it = Qmap.cbegin(); it != Qmap.cend(); ++it)
                {
                    if(it->first=="L")
                    {
                        Q_L=it->second;
                    }
                    else if(it->first=="K1")
                    {
                        Q_K1=it->second;
                    }
                 //cout << it->first << " " << it->second << endl;
                }
                Quad EQ=Quad(QName,Q_L,Q_K1);
                this->Quad_list_.push_back(EQ);
                // for (int j = 0; j < Bparas.size(); j++)
                // {
                //     cout<<Bparas[j]<<endl;
                // }
                searchStart = match.suffix().first;   //来获取从上一个匹配项结束位置到整个字符串末尾的子字符,下一个 std::regex_search 将从这里开始继续查找。
            }

        }
        else if(lines[i].find("DRIFT") != string::npos)
        {
            size_t found = lines[i].find("DRIFT");
            lines[i].erase(found,4);
            //cout<<lines[i]<<endl;

            string::const_iterator searchStart(lines[i].cbegin());
            while (regex_search(searchStart, lines[i].cend(), match, pattern)) 
            {
                string DName = match[1];
                string DString = match[2];
                DString.erase(remove_if(DString.begin(), DString.end(), ::isspace), DString.end());
                //cout<< DName<<endl;
                //cout<<DString<<endl;

                stringstream sst(DString);      
                string Dpara;  

                map<string,double> Lmap;

                while (getline(sst, Dpara, ','))      //逐个提取字符串中的单词（以逗号作为分隔符）
                {
                    Dpara.erase(Dpara.find_last_not_of(" \n\r\t") + 1);
                    auto pos = Dpara.find("=");
                    p1=Dpara.substr(0,pos);
                    p2=Dpara.substr(pos+1);
                    vp2=m1.sstod(p2);
                    Lmap[p1] = vp2;
                }
                Drift ED = Drift(DName,vp2);
                this->Drift_list_.push_back(ED);
                // for(auto it = Lmap.cbegin(); it != Lmap.cend(); ++it)
                // {
                //  cout << it->first << " " << it->second << endl;
                // }
                // for (int j = 0; j < Bparas.size(); j++)
                // {
                //     cout<<Bparas[j]<<endl;
                // }
                searchStart = match.suffix().first;   //来获取从上一个匹配项结束位置到整个字符串末尾的子字符,下一个 std::regex_search 将从这里开始继续查找。
            }

        }
        else if (lines[i].find("MARK") != string::npos)
        {
           
            size_t found = lines[i].find("MARK");
            lines[i].erase(found, 4);

            string::const_iterator searchStart(lines[i].cbegin());
            while (regex_search(searchStart, lines[i].cend(), match, pattern1))
            {
                string MName = match[1];
                Mark EM = Mark(MName);
                this->Mark_list_.push_back(EM);
                searchStart = match.suffix().first;   //来获取从上一个匹配项结束位置到整个字符串末尾的子字符,下一个 std::regex_search 将从这里开始继续查找。
            }

        }
        
        else if (lines[i].find("SEXT") != string::npos)
        {
            size_t found = lines[i].find("SEXT");
            lines[i].erase(found, 4);
            //cout<<lines[i]<<endl;

            string::const_iterator searchStart(lines[i].cbegin());
            while (regex_search(searchStart, lines[i].cend(), match, pattern))
            {
                string SName = match[1];
                string SString = match[2];
                SString.erase(remove_if(SString.begin(), SString.end(), ::isspace), SString.end());
                //cout<<QName<<endl;
                //cout<<QString<<endl;

                stringstream sst(SString);
                string Spara;

                map<string, double> Smap;

                while (getline(sst, Spara, ','))      //逐个提取字符串中的单词（以空格作为分隔符）
                {
                    Spara.erase(Spara.find_last_not_of(" \n\r\t") + 1);
                    auto pos = Spara.find("=");
                    p1 = Spara.substr(0, pos);
                    p2 = Spara.substr(pos + 1);
                    vp2 = m1.sstod(p2);
                    Smap[p1] = vp2;

                }

                double S_L = 0.0;
                double S_K2 = 0.0;
                for (auto it = Smap.cbegin(); it != Smap.cend(); ++it)
                {
                    if (it->first == "L")
                    {
                        S_L = it->second;
                    }
                    else if (it->first == "K2")
                    {
                        S_K2 = it->second;
                    }
                    //cout << it->first << " " << it->second << endl;
                }
                Sext ES = Sext(SName, S_L, S_K2);
                this->Sext_list_.push_back(ES);
                // for (int j = 0; j < Bparas.size(); j++)
                // {
                //     cout<<Bparas[j]<<endl;
                // }
                searchStart = match.suffix().first;   //来获取从上一个匹配项结束位置到整个字符串末尾的子字符,下一个 std::regex_search 将从这里开始继续查找。
            }
    
        }
          
        else if(lines[i].find("LINE") != string::npos)
        {
            size_t found = lines[i].find("LINE");
            lines[i].erase(found,4);
            string::const_iterator searchStart(lines[i].cbegin());
            
            //map<string,string> BLmap;

            int size_BL=0;

            string Ename;

            while (regex_search(searchStart, lines[i].cend(), match, pattern)) 
            {
                
                string BLName = match[1];
                string BLString = match[2];
                BLString.erase(remove_if(BLString.begin(), BLString.end(), ::isspace), BLString.end());
                //cout<<BLName<<endl;
                //cout<<BLString<<endl;
                this->BL_list_.emplace_back(BeamLine());
                this->BL_list_[size_BL].BL_name_=BLName;
                istringstream BLiss(BLString);
                    while (getline(BLiss, Ename, ','))      //提出LINE中的字符串Ename，并从之前储存的list_里寻找对应原件，并放入该BL的Elementlist里面
                    {
                       for (int m = 0; m < this->Bend_list_.size(); m++)
                       {
                            if(this->Bend_list_[m].E_name_==Ename)
                            {
                                string Bname=this->Bend_list_[m].E_name_;
                                double BL= this->Bend_list_[m].length_;
                                double Angle= this->Bend_list_[m].ANGLE_;
                                double K0= this->Bend_list_[m].K0_;
                                double K1= this->Bend_list_[m].K1_;
                                double E1= this->Bend_list_[m].E1_;
                                double E2= this->Bend_list_[m].E2_;


                                this->BL_list_[size_BL].Element_list_.push_back(make_shared<Bend>(Bname,BL,Angle,K0,K1,E1,E2));
                                continue;
                            }
                        }    

                        for (int m = 0; m < this->Quad_list_.size(); m++)
                        {
                            if(this->Quad_list_[m].E_name_==Ename)
                            {
                                string QN=  this->Quad_list_[m].E_name_;
                                double QL=  this->Quad_list_[m].length_;
                                double QK=  this->Quad_list_[m].K1_;
            
                                this->BL_list_[size_BL].Element_list_.push_back(make_shared<Quad>(QN,QL,QK));
                                continue;
                            }
                        }
                        
                        for (int m = 0; m < this->Sext_list_.size(); m++)
                        {
                            if (this->Sext_list_[m].E_name_ == Ename)
                            {
                                string SN = this->Sext_list_[m].E_name_;
                                double SL = this->Sext_list_[m].length_;
                                double SK = this->Sext_list_[m].K2_;

                                this->BL_list_[size_BL].Element_list_.push_back(make_shared<Sext>(SN, SL, SK));
                                continue;
                            }
                        }
                       
                        for (int m = 0; m < this->Drift_list_.size(); m++)
                        {    
                            if(this->Drift_list_[m].E_name_==Ename)
                            {
                                string DN=this->Drift_list_[m].E_name_;
                                double DL=this->Drift_list_[m].length_;
                                this->BL_list_[size_BL].Element_list_.push_back(make_shared<Drift>(DN, DL));
                                continue;
                            }                            
                        }

                        for (int m = 0; m < this->Mark_list_.size(); m++)
                        {
                            if (this->Mark_list_[m].E_name_ == Ename)
                            {
                                string MN = this->Mark_list_[m].E_name_;
                                this->BL_list_[size_BL].Element_list_.push_back(make_shared<Mark>(MN));
                                continue;
                            }
                        }
                        
                        for (int m = 0; m < this->BL_list_.size()-1; m++)
                        {
                            int strfound = Ename.find(this->BL_list_[m].BL_name_);
                            if (strfound != string::npos) 
                            {
                                 if(this->BL_list_[m].BL_name_==Ename)
                                {
                                    for (int k = 0; k < BL_list_[m].Element_list_.size(); k++)
                                    {
                                        auto p = BL_list_[m].Element_list_[k];
                                        auto element = *(p.get());
                                       
                                        //auto node_knot.exit_to = std::make_shared<decltype(element)>();
                                        this->BL_list_[size_BL].Element_list_.push_back(std::make_shared<decltype(element)>(element));
                                        this->BL_list_[size_BL].Element_list_[k]->calc_TM();
                                    }
                                    continue;
                                }

                                else if (this->BL_list_[m].BL_name_!=Ename && this->BL_list_[m].BL_name_==Ename.substr(strfound) && strfound>1 && Ename[strfound-1]=='*')
                                {
                                    
                                    int cell_num=stoi(Ename.substr(0,strfound-1));
                                     for ( int j = 0; j < cell_num; j++)
                                     {
                                         for (int k = 0; k < BL_list_[m].Element_list_.size(); k++)
                                            {
                                            this->BL_list_[size_BL].Element_list_.push_back(BL_list_[m].Element_list_[k]);
                                            }
                                     }
                                    continue;
                                    cout<<m<<endl;
                                    cout<<this->BL_list_[m].BL_name_<<endl;
                                    cout<< Ename<<endl;
                                    cout<< cell_num<<endl;
                                }
                            }
                        }                                                                    
                    }
                searchStart = match.suffix().first; 
                size_BL=size_BL+1;
            }                      
        }
        else if(lines[i].find("USE") != string::npos)
        {
            size_t found = lines[i].find("USE");
            lines[i].erase(found,3);
            lines[i].erase(remove_if(lines[i].begin(), lines[i].end(), ::isspace), lines[i].end()); 
            string mainBL=lines[i];

            for (int m = 0; m < this->BL_list_.size(); m++)
            {
                if (mainBL==this->BL_list_[m].BL_name_)
                {
                    this->BL_list_[m].Element_list_.emplace(this->BL_list_[m].Element_list_.begin(), make_shared<Mark>("^^^"));
                    this->BL_list_[m].Element_list_.push_back(make_shared<Mark>("$$$"));
                    this->Main_List_=this->BL_list_[m];
                }
            }

            double c_position = 0.0;    // calculate the current position of each element in the beamline
            for (int i = 0; i < this->Main_List_.Element_list_.size(); i++)
            {
                this->Main_List_.Element_list_[i]->s_ = c_position;
                c_position= this->Main_List_.Element_list_[i]->length_+ c_position;
            }


        }
    }
}


void Beamlines::show()
{
    cout<<this->Main_List_.BL_name_<<endl;
    for (int i = 0; i < this->Main_List_.Element_list_.size(); i++)
    { 
        cout<<this->Main_List_.Element_list_[i]->E_name_<<" ";
    }
    cout<<endl;
}



 void Beamlines::Calc_L()
 {
    double length;
    length = 0.0;
    for (int i = 0; i < this->Main_List_.Element_list_.size(); i++)
    {
     /* code */
     length=length+this->Main_List_.Element_list_[i]->length_;
    }
    cout<<"Total length is "<<length<<endl;

 }





void Beamlines::Find_Orbit()
{
    double prec=1.0*pow(10, -8);

    double x=0.0;
    double px=0.0;
    double y=0.0;
    double py=0.0;
    double delta=0.0;
    double tau=0.0;   // used for saving closed orbit, inital point is 0,0,0,0,0,0
    
    Matrix<double,6,6> xd;
    xd.setZero();
    Matrix<double,4,4> Mat;
    Mat.setZero();
    Matrix<double,6,1> x0;
    Matrix<double,4,1> xdf;
    Matrix<double,4,1> Del;
    Matrix<double, 4, 4> MatI;
    for (int k = 0; k < 5; k++)
    {
        
    
        double x1=x;
        double px1=px;
        double y1=y;
        double py1=py;
        double delta1=delta;
        double tau1=tau;   // used as a probe

    // track for the intial point 

        for (int j = 0; j < this->Main_List_.Element_list_.size(); j++)   //  track one turn 
        {
            this->Main_List_.Element_list_[j]->integration(x1,px1,y1,py1,delta1,tau1);
        }    

        x0(0,0)=x1;        // first -> raw , second -> column
        x0(1,0)=px1;
        x0(2,0)=y1;
        x0(3,0)=py1;
        x0(4,0)=delta1;
        x0(5,0)=tau1;

        xdf(0,0)=x-x1;        // first -> raw , second -> column
        xdf(1,0)=px-px1;
        xdf(2,0)=y-y1;
        xdf(3,0)=py-py1;



    // track for the intial point with perturbation to get M components 
        for (int i = 0; i < 4; i++)
        {
            xd(i,0)=x;
            xd(i,1)=px;
            xd(i,2)=y;
            xd(i,3)=py;
            xd(i,4)=delta;
            xd(i,5)=tau;
            xd(i,i)=xd(i,i)+prec;    

        
            x1=xd(i,0);
            px1=xd(i,1);
            y1=xd(i,2);
            py1=xd(i,3);
            delta1=xd(i,4);
            tau1=xd(i,5);       // probe (initial point+ perturbation)


            for (int j = 0; j < this->Main_List_.Element_list_.size(); j++)   //  track one turn 
            {
                this->Main_List_.Element_list_[j]->integration(x1,px1,y1,py1,delta1,tau1);
            }
        
            xd(i,0)=x1;
            xd(i,1)=px1;
            xd(i,2)=y1;
            xd(i,3)=py1;
            xd(i,4)=delta1;
            xd(i,5)=tau1;        //save value  

            for (int j = 0; j < 4; j++)
            {
                Mat(j,i)=(xd(i,j)-x0(j,0))/prec;
            }
        }
    
    // get (M-I)
        for (int i = 0; i < 4; i++)
        {
            Mat(i,i)=Mat(i,i)-1.0;
        }

        MatI =Mat.inverse();  // invert matrix
       
        Del=MatI*xdf;

        x=Del(0,0)+x;       // get the new closed orbit information
        px=Del(1,0)+px;
        y=Del(2,0)+y;
        py=Del(3,0)+py;
        delta=0.0;
        tau=0.0;   
    
    }
    
   
    for (int i = 0; i < this->Main_List_.Element_list_.size(); i++)   //  track one turn 
    {
            
        this->Main_List_.Element_list_[i]->x0_=x;
        this->Main_List_.Element_list_[i]->px0_=px;
        this->Main_List_.Element_list_[i]->y0_=y;
        this->Main_List_.Element_list_[i]->py0_=py;
        this->Main_List_.Element_list_[i]->delta0_=delta;
        this->Main_List_.Element_list_[i]->tau0_=tau;


        cout << this->Main_List_.Element_list_[i]->E_name_ << endl;
        cout<<this->Main_List_.Element_list_[i]->x0_<<endl;
        cout<<this->Main_List_.Element_list_[i]->px0_<<endl;
        //cout << this->Main_List_.Element_list_[i]->delta0_ << endl;
        //cout << this->Main_List_.Element_list_[i]->tau0_ << endl;
        this->Main_List_.Element_list_[i]->integration(x,px,y,py,delta,tau);            
    }
   

    
}   




void Beamlines::One_turn_Matrix()
{
    this->Mat_.setZero();
    this->Mat_(0,0)=1.0;
    this->Mat_(1,1)=1.0;
    this->Mat_(2,2)=1.0;
    this->Mat_(3,3)=1.0;
    this->Mat_(4,4)=1.0;
    this->Mat_(5,5)=1.0;
    

    for (int i = 0; i < this->Main_List_.Element_list_.size(); i++)   //  track one turn  this->Main_List_.Element_list_.size()
        {
            //cout<<this->Main_List_.Element_list_[i]->E_name_<<endl;
            //cout<<this->Main_List_.Element_list_[i]->x0_<<endl;
            //cout<<this->Main_List_.Element_list_[i]->px0_<<endl;
            //cout<<this->Main_List_.Element_list_[i]->y0_<<endl;
            //cout<<this->Main_List_.Element_list_[i]->py0_<<endl;
            //cout<<this->Main_List_.Element_list_[i]->delta0_<<endl;
            //cout<<this->Main_List_.Element_list_[i]->tau0_<<endl;
            //cout<<i<<"th:"<<Main_List_.Element_list_[i]<<endl;
            this->Main_List_.Element_list_[i]->calc_TM();
            //cout << this->Main_List_.Element_list_[i]->E_name_ << endl;
            //cout<<this->Main_List_.Element_list_[i]->TM_<<endl;
            this->Mat_=this->Main_List_.Element_list_[i]->TM_*this->Mat_;
            //cout<<this->Mat_<<endl;
        }
   //cout << this->Main_List_.Element_list_[7]->E_name_ << endl;
   //cout<<this->Main_List_.Element_list_[7]->TM_;
   //cout << "One turn matrix:" << endl;
   //cout<<this->Mat_<<endl;
   //MyMath::Check_Sympletic(this->Main_List_.Element_list_[7]->TM_);
}


void Beamlines::Calc_Twiss_nocavi()
{
    vector <double> BL_s;   //  a vector for s of each slice.
    vector<double> BetaX;
    vector<double> BetaY;
    vector<double> AlphaX;
    vector<double> AlphaY;
    vector<double> EtaX;
    vector<double> EtaY;
    double AX;
    double BX;
    double EX;
    double AY;
    double BY;
    double EY;

    Matrix <double, 4, 4> N1 = this->Mat_.block(0, 0, 4, 4);
    //cout << N << endl;
    Matrix <double, 4, 1> v1 = this->Mat_.block(0, 4, 4, 1);
    //cout << v << endl;
    Matrix <double, 4, 4> I4 = MatrixXd::Identity(4, 4);
    auto eta1 = (I4 - N1).inverse() * v1;
    Matrix <double, 6, 6> F1 = MatrixXd::Identity(6, 6);

    for (int i = 0; i < 4; i++)
    {
        F1(i, 4) = eta1(i, 0);
    }

    F1(5, 0) = eta1(1, 0);
    F1(5, 1) = -eta1(0, 0);
    F1(5, 2) = eta1(3, 0);
    F1(5, 3) = -eta1(2, 0);

    // calculate slip factor:
    double alpha_s = this->Mat_(5, 4);
    for (int i = 0; i < 4; i++)
    {
        alpha_s = this->Mat_(5, i) * F1(i, 4) - F1(5, i) * this->Mat_(i, 4) + alpha_s;
    }
    cout << "Slip factor" << endl;
    cout << alpha_s << endl;

    MatrixXd Ai(6, 6);
    MatrixXd A(6, 6);
    MatrixXd Mr(6, 6);
    Mr = F1.inverse() * this->Mat_ * F1;
    Ai = MyMath::Initialize_Ai_nocavi(Mr);
    A =  MyMath::Canonized_A_nocavi(Ai.inverse());
    //cout << A << endl;
    //////   !!!This is R matrix   //////
    double pi = acos(-1);
    //Matrix <double, 6, 6> R = A * F1.inverse() * this->Mat_ * F1 * A.inverse();
    //cout << R << endl;
    // 
    //double mux, muy;
    //mux = acos(R(0, 0)) / 2 / pi;
    //cout << "Nx: " << mux << endl;
    //muy = acos(R(2, 2)) / 2 / pi;
    //cout << "Ny: " << muy << endl;

    // this is for the marker ^^^
    this->Main_List_.Element_list_[0]->BX_ = A(0, 0) * A(0, 0);
    this->Main_List_.Element_list_[0]->AX_ = -A(0, 0) * A(1, 0);
    this->Main_List_.Element_list_[0]->EX_ = eta1(0, 0);
    this->Main_List_.Element_list_[0]->BY_ = A(2, 2) * A(2, 2);
    this->Main_List_.Element_list_[0]->AY_ = eta1(2, 0);
    this->Main_List_.Element_list_[0]->EY_ = A(2, 4);
    // this is for the element after the marker ^^^
    this->Main_List_.Element_list_[1]->AX_ = this->Main_List_.Element_list_[0]->AX_;
    this->Main_List_.Element_list_[1]->BX_ = this->Main_List_.Element_list_[0]->BX_;
    this->Main_List_.Element_list_[1]->AY_ = this->Main_List_.Element_list_[0]->AY_;
    this->Main_List_.Element_list_[1]->BY_ = this->Main_List_.Element_list_[0]->BY_;
    this->Main_List_.Element_list_[1]->EX_ = this->Main_List_.Element_list_[0]->EX_;
    this->Main_List_.Element_list_[1]->EY_ = this->Main_List_.Element_list_[0]->EY_;


    BetaX.push_back(this->Main_List_.Element_list_[1]->BX_);
    BetaY.push_back(this->Main_List_.Element_list_[1]->BY_);
    AlphaX.push_back(this->Main_List_.Element_list_[1]->AX_);
    AlphaY.push_back(this->Main_List_.Element_list_[1]->AY_);
    EtaX.push_back(this->Main_List_.Element_list_[1]->EX_);
    EtaY.push_back(this->Main_List_.Element_list_[1]->EY_);
    BL_s.push_back(this->Main_List_.Element_list_[1]->s_);



    // calulated by slices after ^^^
    Matrix <double, 6, 6> Mn;
    Mn = this->Mat_;
    Matrix <double, 6, 6> Mi;

    Matrix <double, 6, 6> Mrij;
    Matrix <double, 6, 6> A1=A;
    Matrix <double, 6, 6> Rij;
    Matrix <double, 4, 4> N2 ;
    Matrix <double, 4, 1> v2 ;
    
    Matrix <double, 6, 6> F2 = MatrixXd::Identity(6, 6);
    Matrix <double, 6, 6> F2i;
    double Nx = 0.0;
    double Ny = 0.0;

    for (int i = 1; i < this->Main_List_.Element_list_.size() - 1; i++)  //this->Main_List_.Element_list_.size() - 1
    {
        auto Mlist = this->Main_List_.Element_list_[i]->getdata();
        int Count = 1;
        for (const auto& m : Mlist)
        {
            Mi = m * Mn * m.inverse();

            N2 = Mi.block(0, 0, 4, 4);
            v2 = Mi.block(0, 4, 4, 1);
            auto eta2 = (I4 - N2).inverse() * v2;

            for (int j = 0; j < 4; j++)
            {
                F2(j, 4) = eta2(j, 0);
            }

            F2(5, 0) = eta2(1, 0);
            F2(5, 1) = -eta2(0, 0);
            F2(5, 2) = eta2(3, 0);
            F2(5, 3) = -eta2(2, 0);
            F2i = F2.inverse();
            Mr = F2i * Mi * F2;
            Ai = MyMath::Initialize_Ai_nocavi(Mr);
            A = MyMath::Canonized_A_nocavi(Ai.inverse());

            BX = A(0, 0) * A(0, 0);
            AX = -A(0, 0) * A(1, 0);
            EX = eta2(0, 0);
            BY = A(2, 2) * A(2, 2);
            AY = -A(2, 2) * A(3, 2);
            EY = eta2(2, 0);
            BetaX.push_back(BX);
            BetaY.push_back(BY);
            AlphaX.push_back(AX);
            AlphaY.push_back(AY);
            EtaX.push_back(EX);
            EtaY.push_back(EY);

            BL_s.push_back(this->Main_List_.Element_list_[i]->s_ + this->Main_List_.Element_list_[i]->Lslice_ * Count);
            Count = Count + 1;
        }
        //cout << this->Main_List_.Element_list_[i]->E_name_ << endl;
        //cout << this->Main_List_.Element_list_[i]->TM_ << endl;

        Mrij = F2i * this->Main_List_.Element_list_[i]->TM_ * F1;
        //Mrij = Main_List_.Element_list_[i]->TM_;
        Rij = A.inverse() * Mrij* A1;

        if (this->Main_List_.Element_list_[i]->TypeID_ == 21)
        {

        }
        else
        {

        }


        this->Main_List_.Element_list_[i + 1]->BX_ = BX;
        this->Main_List_.Element_list_[i + 1]->AX_ = AX;
        this->Main_List_.Element_list_[i + 1]->EX_ = EX;
        this->Main_List_.Element_list_[i + 1]->BY_ = BY;
        this->Main_List_.Element_list_[i + 1]->AY_ = AY;
        this->Main_List_.Element_list_[i + 1]->EY_ = EY;
        Nx = Nx + acos(abs(Rij(0, 0))) / 2 / pi;
        Ny = Ny + acos(abs(Rij(2, 2))) / 2 / pi;
        this->Main_List_.Element_list_[i + 1]->NX_ = Nx;
        this->Main_List_.Element_list_[i + 1]->NY_ = Ny;

        Mn = Mi;
        A1 = A;
        F1 = F2;
    }

    for (int i = 0; i < this->Main_List_.Element_list_.size(); i++)
    {
        cout << this->Main_List_.Element_list_[i]->E_name_ << ",BX: " << this->Main_List_.Element_list_[i]->BX_
            << ",AX: " << this->Main_List_.Element_list_[i]->AX_ << ",EX: " << this->Main_List_.Element_list_[i]->EX_ << ",NX: " << this->Main_List_.Element_list_[i]->NX_
            << ",BY: " << this->Main_List_.Element_list_[i]->BY_ << ",AY: " << this->Main_List_.Element_list_[i]->AY_
            << ",EY: " << this->Main_List_.Element_list_[i]->EY_ << ",NY: " << this->Main_List_.Element_list_[i]->NY_ << endl;
    }

    plt::subplot(2, 1, 1);
    plt::named_plot("Beta_X", BL_s, BetaX);
    plt::named_plot("Beta_Y", BL_s, BetaY);
    plt::xlabel("s[m]");
    plt::ylabel("Beta[m]");
    plt::legend();
    plt::subplot(2, 1, 2);
    plt::named_plot("Eta_X", BL_s, EtaX);
    plt::xlabel("s[m]");
    plt::ylabel("Eta [m]");
    plt::legend();
    plt::show();


}



void Beamlines::Calc_Twiss()
{
    //MatrixXd Mt(6, 6);
    //Mt = this->Mat_.transpose();
    //cout << Mt << endl;
    //cout << "eigenvalue:" << endl << Mt.eigenvalues() << endl;
    //EigenSolver<MatrixXd> es(Mt);
    //    cout <<"eigenvectors:"<<endl<< es.eigenvectors() << endl;

    MatrixXd Ai(6, 6);
    Ai = MyMath::Initialize_Ai_nocavi(this->Mat_);
    //pb1 = Ai(0,0)*Ai(1,1)-Ai(0,1)*Ai(1,0)+Ai(0,2)*Ai(1,3)-Ai(0,3)*Ai(1,2)+Ai(0,4)*Ai(1,5)-Ai(0,5)*Ai(1,4);
    //pb2 = Ai(2,0)*Ai(3,1)-Ai(2,1)*Ai(3,0)+Ai(2,2)*Ai(3,3)-Ai(2,3)*Ai(3,2)+Ai(2,4)*Ai(3,5)-Ai(2,5)*Ai(3,4);
    //pb3 = Ai(4,0)*Ai(5,1)-Ai(4,1)*Ai(5,0)+Ai(4,2)*Ai(5,3)-Ai(4,3)*Ai(5,2)+Ai(4,4)*Ai(5,5)-Ai(4,5)*Ai(5,4);

    //rotation matrix, make A12, A34=0, A56=0 and A11>0 A33>0 A55>0
    MatrixXd Ri(6, 6);
    Ri.setZero();
    double cs1 = sqrt(Ai(1, 1) * Ai(1, 1) / (Ai(0, 1) * Ai(0, 1) + Ai(1, 1) * Ai(1, 1)));
    double sn1 = sqrt(Ai(0, 1) * Ai(0, 1) / (Ai(0, 1) * Ai(0, 1) + Ai(1, 1) * Ai(1, 1)));
    if (abs(cs1 * Ai(0, 1) - sn1 * Ai(1, 1)) > 0.00001)
    {
        cs1 = -cs1;
    }
    if (cs1 * Ai(0, 0) - sn1 * Ai(1, 0) < 0)
    {
        cs1 = -cs1;
        sn1 = -sn1;
    }


    double cs2 = sqrt(Ai(3, 3) * Ai(3, 3) / (Ai(2, 3) * Ai(2, 3) + Ai(3, 3) * Ai(3, 3)));
    double sn2 = sqrt(Ai(2, 3) * Ai(2, 3) / (Ai(2, 3) * Ai(2, 3) + Ai(3, 3) * Ai(3, 3)));
    if (abs(cs2 * Ai(2, 3) - sn2 * Ai(3, 3)) > 0.00001)
    {
        cs2 = -cs2;
    }
    if (cs2 * Ai(2, 2) - sn2 * Ai(3, 2) < 0)
    {
        cs2 = -cs2;
        sn2 = -sn2;
    }


    double cs3 = sqrt(Ai(5, 5) * Ai(5, 5) / (Ai(4, 5) * Ai(4, 5) + Ai(5, 5) * Ai(5, 5)));
    double sn3 = sqrt(Ai(4, 5) * Ai(4, 5) / (Ai(4, 5) * Ai(4, 5) + Ai(5, 5) * Ai(5, 5)));
    if (abs(cs3 * Ai(4, 5) - sn3 * Ai(5, 5)) > 0.00001)
    {
        cs3 = -cs3;
    }
    if (cs3 * Ai(4, 4) - sn3 * Ai(5, 4) < 0)
    {
        cs3 = -cs3;
        sn3 = -sn3;
    }

    Ri(0, 0) = cs1;
    Ri(0, 1) = -sn1;
    Ri(1, 0) = sn1;
    Ri(1, 1) = cs1;
    Ri(2, 2) = cs2;
    Ri(2, 3) = -sn2;
    Ri(3, 2) = sn2;
    Ri(3, 3) = cs2;
    Ri(4, 4) = cs3;
    Ri(4, 5) = -sn3;
    Ri(5, 4) = sn3;
    Ri(5, 5) = cs3;

    cout << "rotation Ri:" << endl;
    cout << Ri << endl;

    MatrixXd Ai_Cano(6, 6);
    Ai_Cano = Ri * Ai;

    cout << "Ai_Cano:" << endl;
    cout << Ai_Cano << endl;

    cout << "BetaX: " << 1 / (Ai_Cano(0, 0) * Ai_Cano(0, 0)) << endl;
    cout << "AlphaX: " << Ai_Cano(1, 0) / Ai_Cano(0, 0) << endl;
    cout << "BetaY: " << 1 / (Ai_Cano(2, 2) * Ai_Cano(2, 2)) << endl;
    cout << "AlphaY: " << Ai_Cano(3, 2) / Ai_Cano(2, 2) << endl;

    // this is for the marker ^^^
    this->Main_List_.Element_list_[0]->AX_ = Ai_Cano(1, 0) / Ai_Cano(0, 0);
    this->Main_List_.Element_list_[0]->BX_ = 1 / (Ai_Cano(0, 0) * Ai_Cano(0, 0));
    this->Main_List_.Element_list_[0]->AY_ = Ai_Cano(3, 2) / Ai_Cano(2, 2);
    this->Main_List_.Element_list_[0]->BY_ = 1 / (Ai_Cano(2, 2) * Ai_Cano(2, 2));
    MatrixXd R(6, 6);
    MatrixXd A(6, 6);
    R.setZero();
    A = Ai_Cano.inverse();
    cout << "A:" << endl;
    cout << A << endl;

    this->Main_List_.Element_list_[0]->EX_ = A(0, 4);
    this->Main_List_.Element_list_[0]->EY_ = A(2, 4);
    R = Ai_Cano * (this->Mat_ * A);

    this->Main_List_.Element_list_[1]->AX_ = this->Main_List_.Element_list_[0]->AX_;
    this->Main_List_.Element_list_[1]->BX_ = this->Main_List_.Element_list_[0]->BX_;
    this->Main_List_.Element_list_[1]->AY_ = this->Main_List_.Element_list_[0]->AY_;
    this->Main_List_.Element_list_[1]->BY_ = this->Main_List_.Element_list_[0]->BY_;
    this->Main_List_.Element_list_[1]->EX_ = this->Main_List_.Element_list_[0]->EX_;
    this->Main_List_.Element_list_[1]->EY_ = this->Main_List_.Element_list_[0]->EY_;

    cout << "rotation R:" << endl;
    cout << R << endl;
    // 
    //double mux, muy;
    double pi = acos(-1);
    //mux = acos(R(0, 0)) / 2 / pi;
    //cout << "Nx: " << mux << endl;
    //muy = acos(R(2, 2)) / 2 / pi;
    //cout << "Ny: " << muy << endl;

    MatrixXd Rij(6, 6);
    MatrixXd Ak(6, 6);
    MatrixXd Ac(6, 6);
    double Nx = 0.0;
    double Ny = 0.0;
    Rij.setZero();
    Ak.setZero();
    Ac.setZero();
    /*
    for (int i = 0; i < this->Main_List_.Element_list_.size() - 1; i++)
    {
        
        //cout << "Transfer Matrix" << endl;
        //cout << this->Main_List_.Element_list_[i]->TM_<<endl;
        Ak = this->Main_List_.Element_list_[i]->TM_ * A;
        A = MyMath::Canonized_A(Ak);
        Rij = Ak.inverse() * A;
        //cout << this->Main_List_.Element_list_[i + 1]->E_name_ << endl;
        //cout << A << endl;
        this->Main_List_.Element_list_[i + 1]->BX_ = A(0, 0) * A(0, 0);
        this->Main_List_.Element_list_[i + 1]->AX_ = -A(0, 0) * A(1, 0);
        this->Main_List_.Element_list_[i + 1]->EX_ = A(0, 4);
        this->Main_List_.Element_list_[i + 1]->BY_ = A(2, 2) * A(2, 2);
        this->Main_List_.Element_list_[i + 1]->AY_ = -A(2, 2) * A(3, 2);
        this->Main_List_.Element_list_[i + 1]->EY_ = A(2, 4);

        Nx = Nx + acos(abs(Rij(0, 0))) / 2 / pi;
        Ny = Ny + acos(abs(Rij(2, 2))) / 2 / pi;
        this->Main_List_.Element_list_[i + 1]->NX_ = Nx;
        this->Main_List_.Element_list_[i + 1]->NY_ = Ny;
    }
    */
    
    /// have a try.....
    vector <double> BL_s;   //  a vector for s of each slice.
    vector<double> BetaX;
    vector<double> BetaY;
    vector<double> AlphaX;
    vector<double> AlphaY;
    vector<double> EtaX;
    vector<double> EtaY;
    double AX;
    double BX;
    double EX;
    double AY;
    double BY;
    double EY;

    BetaX.push_back(this->Main_List_.Element_list_[1]->BX_);
    BetaY.push_back(this->Main_List_.Element_list_[1]->BY_);
    AlphaX.push_back(this->Main_List_.Element_list_[1]->AX_);
    AlphaY.push_back(this->Main_List_.Element_list_[1]->AY_);
    EtaX.push_back(this->Main_List_.Element_list_[1]->EX_);
    EtaY.push_back(this->Main_List_.Element_list_[1]->EY_);
    BL_s.push_back(this->Main_List_.Element_list_[1]->s_);

    for (int i = 1; i < this->Main_List_.Element_list_.size() - 1; i++)
    {
        auto Mlist=this->Main_List_.Element_list_[i]->getdata();        
        int Count = 1;
        for (const auto& m : Mlist)
            {
                Ak = m * A;
                Ac = MyMath::Canonized_A_nocavi(Ak);
                
                AX = -Ac(0, 0) * Ac(1, 0);
                BX = Ac(0, 0) * Ac(0, 0);
                EX = Ac(0, 4);
                AY = -Ac(2, 2) * Ac(3, 2);
                BY = Ac(2, 2) * Ac(2, 2);
                EY = Ac(2, 4);
                BetaX.push_back(BX);
                BetaY.push_back(BY);
                AlphaX.push_back(AX);
                AlphaY.push_back(AY);
                EtaX.push_back(EX);
                EtaY.push_back(EY);

                BL_s.push_back(this->Main_List_.Element_list_[i]->s_ + this->Main_List_.Element_list_[i]->Lslice_ * Count);
                Count = Count + 1;
            }
        A = Ac;
        Rij = Ak.inverse() * A;
        this->Main_List_.Element_list_[i + 1]->BX_ = BX;
        this->Main_List_.Element_list_[i + 1]->AX_ = AX;
        this->Main_List_.Element_list_[i + 1]->EX_ = EX;
        this->Main_List_.Element_list_[i + 1]->BY_ = BY;
        this->Main_List_.Element_list_[i + 1]->AY_ = AY;
        this->Main_List_.Element_list_[i + 1]->EY_ = EY;
        Nx = Nx + acos(abs(Rij(0, 0))) / 2 / pi;
        Ny = Ny + acos(abs(Rij(2, 2))) / 2 / pi;
        this->Main_List_.Element_list_[i + 1]->NX_ = Nx;
        this->Main_List_.Element_list_[i + 1]->NY_ = Ny;

    }
    

    for (int i = 0; i < this->Main_List_.Element_list_.size(); i++)
    {
        cout << this->Main_List_.Element_list_[i]->E_name_ << ",BX: " << this->Main_List_.Element_list_[i]->BX_
            << ",AX: " << this->Main_List_.Element_list_[i]->AX_ << ",EX: " << this->Main_List_.Element_list_[i]->EX_ << ",NX: " << this->Main_List_.Element_list_[i]->NX_
            << ",BY: " << this->Main_List_.Element_list_[i]->BY_ << ",AY: " << this->Main_List_.Element_list_[i]->AY_
            << ",EY: " << this->Main_List_.Element_list_[i]->EY_ << ",NY: " << this->Main_List_.Element_list_[i]->NY_ << endl;
    }
    //cout << "NX:" << endl;
    //for (int i = 0; i < this->Main_List_.Element_list_.size(); i++)
    //{
    //    cout << this->Main_List_.Element_list_[i]->NX_ << endl;
    //}
   
    //plt::plot(s, BetaX);
    //plt::plot(s, BetaY);
    //plt::show();
    plt::subplot(2, 1, 1);
    plt::named_plot("Beta_X", BL_s, BetaX);
    plt::named_plot("Beta_Y", BL_s, BetaY);
    plt::xlabel("s[m]");
    plt::ylabel("Beta[m]");
    plt::legend();
    plt::subplot(2, 1, 2);
    plt::named_plot("Eta_X", BL_s, EtaX);
    plt::xlabel("s[m]");
    plt::ylabel("Eta [m]");
    plt::legend();
    plt::show();
    

}







