
#include "MyMath.h"

map<char, int> pre{{'+',1},{'-',1},{'*',2},{'/',2},{'^',3}};



double MyMath::sstod(string expr) 
{
    string tempNum = "";
    vector<Token> infix, suffix;
    stack<Token> oper;
    //cin>>expr;
    
    expr.erase(remove_if(expr.begin(), expr.end(), ::isspace), expr.end()); 
    expr.push_back('\0');
    //cout<<expr<<endl;

    if (expr[0]=='-' || expr[0]=='.')     // 对于首个字符是 负号或者. 的情况，在其前面加上0
    {expr='0'+expr;}
    for (auto c : expr)  // 遍历expr中的每一个字符， 将字符串解析成若干Token存入vector
    {
        if(isdigit(c) || c=='.')    // isdigit()函数判断是否字符为数字，如果为数字，
            {tempNum+=c;}                 // 则将其添加到临时字符串中
        else if(c=='+' || c=='-' || c=='*' || c=='/' || c=='^' )
            {infix.push_back(Token{0,stod(tempNum)});  //则 op存为0,而num 存为该数字，利用char '0'做了一个计算。 
            infix.push_back(Token{c,0});  //如果不为数字，为符号，则op存为c,而num 设为0
            tempNum="";}
            
        if(c=='\0')  //最后一个字符，之前的tempNum是数字且需要保存
            {infix.push_back(Token{0,stod(tempNum)});}

    }

    for (auto t : infix)
    {
        if (!t.op)
            suffix.push_back(t);   // 如果 op是数字，则将改元素t整个存入 suffix中
        else if (t.op == '(')     // 如果 op是'('，则将改元素t整个存入 stack的栈中
            oper.push(t);
        else if (t.op == ')')   // 如果 op是')'
        {
            while (!oper.empty() && oper.top().op != '(')     // 如果 op是')'
            {
                suffix.push_back(oper.top());          //  pop operators from the oper stack and push them onto the suffix vector until you encounter an opening parenthesis '('
                oper.pop();      // 将oper末尾的token弹出
            } 
            oper.pop();
        }
        else
        {
            while (!oper.empty() && pre[oper.top().op] >= pre[t.op])       
        //If the token is an operator (+, -, *, /, ^), you need to check its precedence and associativity to determine its position in the postfix expression.
            {
                suffix.push_back(oper.top());
                oper.pop();
            }
            oper.push(t);
        }
    }
    while (!oper.empty())
    {
        suffix.push_back(oper.top());
        oper.pop();
    }

    //   for (int i = 0; i < suffix.size(); i++)
    //   {
    //  //     /* code */
    //       cout<<suffix[i].num<<"    "<<suffix[i].op<<endl;
    //   }
    


    ////这里对后缀表达式的计算
    while (suffix.size() > 1)
    {
        // for (auto t : suffix)
        //     if (t.op)
        //         cout << t.op << " ";
        //     else
        //         cout << t.num << " ";
        // cout << endl;
        for (int i = 2; i < suffix.size(); ++i)
        {
            double res;
            Token a = suffix[i - 2], b = suffix[i - 1];
            switch (suffix[i].op)
            {
            case 0:
                continue;
            case '+':
                res = a.num + b.num;
                break;
            case '-':
                res = a.num - b.num;
                break;
            case '*':
                res = a.num * b.num;
                break;
            case '/':
                res = a.num / b.num;
                break;
            case '^':
                res = pow(a.num, b.num);
                break;
            }
            suffix.insert(suffix.begin() + i + 1, Token{0, res});
            suffix.erase(suffix.begin() + i - 2, suffix.begin() + i + 1);
            break;
        }
    }
    //cout << suffix[0].num; // 输出最后结果
    return suffix[0].num;
}

Matrix<double, 6, 6> MyMath::Initialize_Ai_nocavi(Matrix<double, 6, 6> M)
{
    MatrixXd Mt(6, 6);
    Mt = M.transpose();
    EigenSolver<MatrixXd> es(Mt);
    //cout << Mt.eigenvalues() << endl;
    //cout << real(Mt.eigenvalues()(0)) << endl;

    MatrixXd Ai(6, 6);
    Ai.setZero();

    int a = 0;
    const double epsilon = 1e-7;

    for (int i = 0; i < 6; i++)
    {
        if (abs(real(Mt.eigenvalues()(i)) - 1.0) > epsilon)
        {
            Ai(0, 0) = real(es.eigenvectors()(0, i));
            Ai(0, 1) = real(es.eigenvectors()(1, i));
            Ai(0, 2) = real(es.eigenvectors()(2, i));
            Ai(0, 3) = real(es.eigenvectors()(3, i));
            Ai(0, 4) = real(es.eigenvectors()(4, i));
            Ai(0, 5) = real(es.eigenvectors()(5, i));
            Ai(1, 0) = imag(es.eigenvectors()(0, i));
            Ai(1, 1) = imag(es.eigenvectors()(1, i));
            Ai(1, 2) = imag(es.eigenvectors()(2, i));
            Ai(1, 3) = imag(es.eigenvectors()(3, i));
            Ai(1, 4) = imag(es.eigenvectors()(4, i));
            Ai(1, 5) = imag(es.eigenvectors()(5, i));
            break;
        }
        a = a + 1;
    }



    for (int i = 1; i < 6; i++)
    {
        if (  (real(Mt.eigenvalues()(i)) != real(Mt.eigenvalues()(a))) && (abs(real(Mt.eigenvalues()(i)) - 1.0) > epsilon)  )
        {
            //cout << Mt.eigenvalues()(i) << endl;
            //cout << real(Mt.eigenvalues()(i)) << endl;
            Ai(2, 0) = real(es.eigenvectors()(0, i));
            Ai(2, 1) = real(es.eigenvectors()(1, i));
            Ai(2, 2) = real(es.eigenvectors()(2, i));
            Ai(2, 3) = real(es.eigenvectors()(3, i));
            Ai(2, 4) = real(es.eigenvectors()(4, i));
            Ai(2, 5) = real(es.eigenvectors()(5, i));
            Ai(3, 0) = imag(es.eigenvectors()(0, i));
            Ai(3, 1) = imag(es.eigenvectors()(1, i));
            Ai(3, 2) = imag(es.eigenvectors()(2, i));
            Ai(3, 3) = imag(es.eigenvectors()(3, i));
            Ai(3, 4) = imag(es.eigenvectors()(4, i));
            Ai(3, 5) = imag(es.eigenvectors()(5, i));
            break;
        }
        a = a + 1;
    }
            Ai(4, 4) = 1.0;
            Ai(5, 5) = 1.0;

    //cout << "eigenvectors:" << endl << es.eigenvectors() << endl;
    //cout << "Ai before:" << endl;
    //cout << Ai << endl;

    double pb1, pb2, pb3;
    pb1 = Ai(0, 0) * Ai(1, 1) - Ai(0, 1) * Ai(1, 0) + Ai(0, 2) * Ai(1, 3) - Ai(0, 3) * Ai(1, 2) + Ai(0, 4) * Ai(1, 5) - Ai(0, 5) * Ai(1, 4);
    pb2 = Ai(2, 0) * Ai(3, 1) - Ai(2, 1) * Ai(3, 0) + Ai(2, 2) * Ai(3, 3) - Ai(2, 3) * Ai(3, 2) + Ai(2, 4) * Ai(3, 5) - Ai(2, 5) * Ai(3, 4);
    pb3 = Ai(4, 0) * Ai(5, 1) - Ai(4, 1) * Ai(5, 0) + Ai(4, 2) * Ai(5, 3) - Ai(4, 3) * Ai(5, 2) + Ai(4, 4) * Ai(5, 5) - Ai(4, 5) * Ai(5, 4);

    // manipulation
    double cp1, cp2, cp3;
    if (pb1 < 0)
    {
        cp1 = sqrt(abs(pb1)) * (-1);
    }
    else
    {
        cp1 = sqrt(pb1);
    }
    if (pb2 < 0)
    {
        cp2 = sqrt(abs(pb2)) * (-1);
    }
    else
    {
        cp2 = sqrt(pb2);
    }
    if (pb3 < 0)
    {
        cp3 = sqrt(abs(pb3)) * (-1);
    }
    else
    {
        cp3 = sqrt(pb3);
    }

    Ai.row(0) = Ai.row(0) / abs(cp1);
    Ai.row(1) = Ai.row(1) / cp1;
    Ai.row(2) = Ai.row(2) / abs(cp2);
    Ai.row(3) = Ai.row(3) / cp2;
    Ai.row(4) = Ai.row(4) / abs(cp3);
    Ai.row(5) = Ai.row(5) / cp3;

    return Ai;
}

Matrix<double, 6, 6> MyMath::Canonized_A_nocavi(Matrix<double, 6, 6> Ak)
{
    MatrixXd Rij(6, 6);
    MatrixXd A(6, 6);
    Rij.setZero();
    double cs1, sn1, cs2, sn2;

    cs1 = sqrt(Ak(0, 0) * Ak(0, 0) / (Ak(0, 0) * Ak(0, 0) + Ak(0, 1) * Ak(0, 1)));
    sn1 = sqrt(Ak(0, 1) * Ak(0, 1) / (Ak(0, 0) * Ak(0, 0) + Ak(0, 1) * Ak(0, 1)));
    if (abs(cs1 * Ak(0, 1) + sn1 * Ak(0, 0)) > 0.00001)
    {
        cs1 = -cs1;
    }
    if (cs1 * Ak(0, 1) + sn1 * Ak(0, 0) < 0)
    {
        cs1 = -cs1;
        sn1 = -sn1;
    }


    cs2 = sqrt(Ak(2, 2) * Ak(2, 2) / (Ak(2, 2) * Ak(2, 2) + Ak(2, 3) * Ak(2, 3)));
    sn2 = sqrt(Ak(2, 3) * Ak(2, 3) / (Ak(2, 2) * Ak(2, 2) + Ak(2, 3) * Ak(2, 3)));
    if (abs(cs2 * Ak(2, 3) + sn2 * Ak(2, 2)) > 0.00001)
    {
        cs2 = -cs2;
    }
    if (cs2 * Ak(2, 3) + sn2 * Ak(2, 2) < 0)
    {
        cs2 = -cs2;
        sn2 = -sn2;
    }
    /*
    cs3 = sqrt(Ak(4, 4) * Ak(4, 4) / (Ak(4, 4) * Ak(4, 4) + Ak(4, 5) * Ak(4, 5)));
    sn3 = sqrt(Ak(4, 5) * Ak(4, 5) / (Ak(4, 4) * Ak(4, 4) + Ak(4, 5) * Ak(4, 5)));
    if (abs(cs3 * Ak(4, 5) + sn3 * Ak(4, 4)) > 0.00001)
    {
        cs3 = -cs3;
    }
    if (cs3 * Ak(4, 5) + sn3 * Ak(4, 4) < 0)
    {
        cs3 = -cs3;
        sn3 = -sn3;
    }
    */
    Rij(0, 0) = cs1;
    Rij(0, 1) = sn1;
    Rij(1, 0) = -sn1;
    Rij(1, 1) = cs1;
    Rij(2, 2) = cs2;
    Rij(2, 3) = sn2;
    Rij(3, 2) = -sn2;
    Rij(3, 3) = cs2;
    Rij(4, 4) = 1;
    Rij(5, 5) = 1;
    //Rij(4, 4) = cs3;
    //Rij(4, 5) = sn3;
    //Rij(5, 4) = -sn3;
    //Rij(5, 5) = cs3;

    A = Ak * Rij;

    return A;
}


void MyMath::Check_Sympletic(Matrix<double, 6, 6> M)
{
    Matrix<double, 6, 6> S;
    S.setZero();
    for (int i = 0; i < 3; i++)
    {
        S(2 * i+1, 2 * i ) = -1.0;
        S(2 * i, 2 * i+1) = 1.0;
    }
    Matrix<double, 6, 6> Mt;
    Mt = M.transpose();
    Matrix<double, 6, 6> res;
    res = Mt*S * M;
    cout << "Mt*S*M:" << endl;
    cout << res << endl;
}
