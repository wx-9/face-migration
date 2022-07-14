#include <iostream>
#include<vector>
#include<fstream>
#include<cstring>
#include<math.h>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<time.h>
using namespace std;
using namespace Eigen;
struct  Point
{
    double num_x;
    double num_y;
   double num_z;
};

struct  nomal_vector
{
  double  num_x;
  double num_y;
  double num_z;  
};

struct  v_index
{
   int fir;
   int sed;
   int thr;  
};
struct mesh
{
    Point my_v;
    nomal_vector f_v;
};

vector<Point> s0_p;
vector<Point> s1_p;
vector<Point> t0_p;
vector<Point> t1_p;
vector<Point> h_p;

vector<v_index> index_1;
vector<v_index> index_2;
vector<v_index> index_3;

vector<MatrixXd> my_S;
vector<MatrixXd> S_d;
vector<MatrixXd> my_T;
vector<MatrixXd> T_d;
 
Point caculate_v4(double x1,double y1,double z1,double x2,double y2,double z2,double x3,double y3,double z3)
{
    Point v4;
    double X1=x2-x1; double Y1=y2-y1;double Z1=z2-z1;
    double X2=x3-x1; double Y2=y3-y1;double Z2=z3-z1;
    double X3=Y1*Z2-Z1*Y2; double Y3=-(X1*Z2-Z1*X2);double Z3=X1*Y2-Y1*X2;
    double div=sqrt(X3*X3+Y3*Y3+Z3*Z3);
    div=sqrt(div); 
    double X4=x1+X3/(div);
    double Y4=y1+Y3/(div);
    double Z4=z1+Z3/(div);
    v4.num_x=X4;v4.num_y=Y4;v4.num_z=Z4;
    return v4; 
}
int add_S()
{
    for (int i=0;i<index_1.size();i++)
    {
        MatrixXd mat0(3,3);
        MatrixXd mat1(3,3);
        MatrixXd d(3,1);
        MatrixXd v(3,1);
        MatrixXd V(3,1);
        MatrixXd S(3,3);

        int fri =index_1[i].fir-1;int sed=index_1[i].sed-1;int thr=index_1[i].thr-1;
        double x1=s0_p[fri].num_x;double y1=s0_p[fri].num_y;double z1=s0_p[fri].num_z;
        double x2=s0_p[sed].num_x;double y2=s0_p[sed].num_y;double z2=s0_p[sed].num_z;
        double x3=s0_p[thr].num_x;double y3=s0_p[thr].num_y;double z3=s0_p[thr].num_z;
        Point v4=caculate_v4(x1,y1,z1,x2,y2,z2,x3,y3,z3);
        double x4=v4.num_x;double y4=v4.num_y;double z4=v4.num_z;

        double X1=s1_p[fri].num_x;double Y1=s1_p[fri].num_y;double Z1=s1_p[fri].num_z;
        double X2=s1_p[sed].num_x;double Y2=s1_p[sed].num_y;double Z2=s1_p[sed].num_z;
        double X3=s1_p[thr].num_x;double Y3=s1_p[thr].num_y;double Z3=s1_p[thr].num_z;
        Point V4=caculate_v4(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3);
        double X4=V4.num_x;double Y4=V4.num_y;double Z4=V4.num_z;

        mat0<<x2-x1,x3-x1,x4-x1,y2-y1,y3-y1,y4-y1,z2-z1,z3-z1,z4-z1;
        mat1<<X2-X1,X3-X1,X4-X1,Y2-Y1,Y3-Y1,Y4-Y1,Z2-Z1,Z3-Z1,Z4-Z1;
        v<<x1,y1,z1;V<<X1,Y1,Z1; 
        S=mat1*mat0.inverse();
        d=V-S*v;
        my_S.push_back(S);
        S_d.push_back(d);

    }
}
int add_S2()
{
    for (int i=0;i<index_1.size();i++)
    {
        MatrixXd mat0(3,3);
        MatrixXd mat1(3,3);
        MatrixXd d(3,1);
        MatrixXd v(3,1);
        MatrixXd V(3,1);
        MatrixXd S(3,3);

        int fri =index_1[i].fir-1;int sed=index_1[i].sed-1;int thr=index_1[i].thr-1;
        double x1=s0_p[fri].num_x;double y1=s0_p[fri].num_y;double z1=s0_p[fri].num_z;
        double x2=s0_p[sed].num_x;double y2=s0_p[sed].num_y;double z2=s0_p[sed].num_z;
        double x3=s0_p[thr].num_x;double y3=s0_p[thr].num_y;double z3=s0_p[thr].num_z;
        Point v4=caculate_v4(x1,y1,z1,x2,y2,z2,x3,y3,z3);
        double x4=v4.num_x;double y4=v4.num_y;double z4=v4.num_z;

        double X1=s1_p[fri].num_x;double Y1=s1_p[fri].num_y;double Z1=s1_p[fri].num_z;
        double X2=s1_p[sed].num_x;double Y2=s1_p[sed].num_y;double Z2=s1_p[sed].num_z;
        double X3=s1_p[thr].num_x;double Y3=s1_p[thr].num_y;double Z3=s1_p[thr].num_z;
        Point V4=caculate_v4(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3);
        double X4=V4.num_x;double Y4=V4.num_y;double Z4=V4.num_z;
        // mat0<<x2-x1,y2-y1,z2-z1,x3-x1,y3-y1,z3-z1,x4-x1,y4-y1,z4-z1;
         mat0<<x2-x1,x3-x1,x4-x1,y2-y1,y3-y1,y4-y1,z2-z1,z3-z1,z4-z1;
        mat1<<X2-X1,X3-X1,X4-X1,Y2-Y1,Y3-Y1,Y4-Y1,Z2-Z1,Z3-Z1,Z4-Z1;
        // mat1<<X2-X1,Y2-Y1,Z2-Z1,X3-X1,Y3-Y1,Z3-Z1,X4-X1,Y4-Y1,Z4-Z1;
        v<<x1,y1,z1;V<<X1,Y1,Z1; 
        S=mat0.inverse()*mat1;
        MatrixXd ST=S.transpose();
         d=V-S*v;
        my_S.push_back(ST);
         S_d.push_back(d);

    }
}


vector <string>  split(const string &str,const string &delim )  
{
    vector<string> res;
    if (""==str) return res;
    char *cstr = new char[str.length()+1];
    char *cdelim=new char[delim.length()+1];
    strcpy(cstr,str.c_str());
    strcpy(cdelim,delim.c_str());
    char *p=strtok(cstr,cdelim);
    while(p)
    {
        string s=p;
        res.push_back(s);
        p=strtok(NULL,cdelim);
    }
return res;
}
int read_obj(string path,vector<Point> &p,vector<v_index> &v_inde)  //读取每个三角片面的顶点和法向量信息
{
    ifstream my_file;
    my_file.open(path);
    if(!my_file)
    {
        return -1; 
    }
    string line;
    int i=10;
    int num_p;
    char head;
    char head1;
    char head2;
    while (getline(my_file,line) && i>0)
    {
        int len_line = line.length();
        head=line[0];
        vector<string>res = split(line," ");
        Point my_p;
        v_index my_index;
        if (head=='v')
        {
             my_p.num_x=stof(res[1]);
             my_p.num_y=stof(res[2]);
             my_p.num_z=stof(res[3]);
             p.push_back(my_p);
        }

        if (head == 'f')
        {
            my_index.fir=stoi(res[1]);
            my_index.sed=stoi(res[2]);
            my_index.thr=stoi(res[3]);
            v_inde.push_back(my_index);
        }
    
    }
    return 0;
}

int write_obj(string path)
{
   ofstream file;
   file.open(path);
   file<<"# 1069 vertices, 1947 faces"<<endl;
   for(int i=0;i<t1_p.size();i++)
   {
     file<<"v"<<" "<<t1_p[i].num_x<<" "<<t1_p[i].num_y<<" "<<t1_p[i].num_z<<endl;
   } 
   for(int i=0;i<index_3.size();i++)
   {
       file<<"f"<<" "<<index_3[i].fir<<" "<<index_3[i].sed<<" "<<index_3[i].thr<<endl;
   }
return 0 ;
    
}


int caculate_target()
{
    SparseMatrix<double> A(3*index_1.size(),index_1.size()+t0_p.size());
    vector<Triplet<double>> tripletLIst;
    //MatrixXd A(3*index_1.size(),index_1.size()+t0_p.size());
    MatrixXd c(3*index_1.size(),3);

    MatrixXd x;
    MatrixXd mat0(3,3);
    Point my_t1;
   // A=A*0;
    for(int i =0;i<index_1.size();i++)
    { //do    transpose  
        c(0+3*i,0)=my_S[i](0,0);c(0+3*i,1)=my_S[i](1,0);c(0+3*i,2)=my_S[i](2,0);
        c(1+3*i,0)=my_S[i](0,1);c(1+3*i,1)=my_S[i](1,1);c(1+3*i,2)=my_S[i](2,1);
        c(2+3*i,0)=my_S[i](0,2);c(2+3*i,1)=my_S[i](1,2);c(2+3*i,2)=my_S[i](2,2);
    }

   for(int j =0;j<index_1.size();j++)
   {
        int fri =index_1[j].fir-1;int sed=index_1[j].sed-1;int thr=index_1[j].thr-1;
        double x1=t0_p[fri].num_x;double y1=t0_p[fri].num_y;double z1=t0_p[fri].num_z;
        double x2=t0_p[sed].num_x;double y2=t0_p[sed].num_y;double z2=t0_p[sed].num_z;
        double x3=t0_p[thr].num_x;double y3=t0_p[thr].num_y;double z3=t0_p[thr].num_z;     
        Point v4=caculate_v4(x1,y1,z1,x2,y2,z2,x3,y3,z3);
        double x4=v4.num_x;double y4=v4.num_y;double z4=v4.num_z;
        mat0<<x2-x1,x3-x1,x4-x1,y2-y1,y3-y1,y4-y1,z2-z1,z3-z1,z4-z1;
        MatrixXd mat0i=mat0.inverse();
        tripletLIst.push_back(Triplet<double>(3*j,fri,-mat0i(0,0)-mat0i(1,0)-mat0i(2,0)));tripletLIst.push_back(Triplet<double>(3*j,sed,mat0i(0,0)));tripletLIst.push_back(Triplet<double>(3*j,thr,mat0i(1,0)));tripletLIst.push_back(Triplet<double>(3*j,t0_p.size()+j,mat0i(2,0)));
        tripletLIst.push_back(Triplet<double>(1+3*j,fri,-mat0i(0,1)-mat0i(1,1)-mat0i(2,1)));tripletLIst.push_back(Triplet<double>(1+3*j,sed,mat0i(0,1)));tripletLIst.push_back(Triplet<double>(1+3*j,thr,mat0i(1,1)));tripletLIst.push_back(Triplet<double>(1+3*j,t0_p.size()+j,mat0i(2,1)));
        tripletLIst.push_back(Triplet<double>(2+3*j,fri,-mat0i(0,2)-mat0i(1,2)-mat0i(2,2)));tripletLIst.push_back(Triplet<double>(2+3*j,sed,mat0i(0,2)));tripletLIst.push_back(Triplet<double>(2+3*j,thr,mat0i(1,2)));tripletLIst.push_back(Triplet<double>(2+3*j,t0_p.size()+j,mat0i(2,2)));
   } 
   A.setFromTriplets(tripletLIst.begin(),tripletLIst.end());
   SparseMatrix<double> A_transpose=A.transpose();
   SparseMatrix<double> ATA=A.transpose()*A;
   MatrixXd ls_b=A_transpose*c;
//    SimplicialCholesky<SparseMatrix<double>>MatricesCholesky(ATA);
SparseLU<SparseMatrix<double>> solver;
    // x=MatricesCholesky.solve(ls_b);
    solver.analyzePattern(ATA);
    solver.factorize(ATA);
    x=solver.solve(ls_b);

   
   for(int k=0;k<t0_p.size();k++)
   {
        my_t1.num_x=x(0+k,0);
        my_t1.num_y=x(0+k,1);
        my_t1.num_z=x(0+k,2);
        t1_p.push_back(my_t1);
   }
   return 0;
}

int caculate_target_2()
{
    // MatrixXd V(3*index_1.size(),3*index_1.size());
    // MatrixXd A(3*index_1.size(),index_1.size()+t0_p.size());
        SparseMatrix<double> A(3*index_1.size(),index_1.size()+t0_p.size());
        SparseMatrix<double> V(3*index_1.size(),3*index_1.size());
           SparseMatrix<double> V1(3*index_1.size(),3*index_1.size());
    vector<Triplet<double>> tripletLIst;
     vector<Triplet<double>> tripletLIst1;
     vector<Triplet<double>> tripletLIst2;



    MatrixXd c(3*index_1.size(),3);
    c.setZero();
     MatrixXd S(3*index_1.size(),3);
     MatrixXd s0(index_1.size()+t0_p.size(),3);
    SparseMatrix<double> D;
    MatrixXd x(index_1.size()+t0_p.size(),3);
    MatrixXd mat0(3,3);
    MatrixXd mat1(3,3);
    Point my_t1;
    for(int i =0;i<index_1.size();i++)
    { //do    transpose  
        c(0+3*i,0)=my_S[i](0,0);c(0+3*i,1)=my_S[i](0,1);c(0+3*i,2)=my_S[i](0,2);
        c(1+3*i,0)=my_S[i](1,0);c(1+3*i,1)=my_S[i](1,1);c(1+3*i,2)=my_S[i](1,2);
        c(2+3*i,0)=my_S[i](2,0);c(2+3*i,1)=my_S[i](2,1);c(2+3*i,2)=my_S[i](2,2);
    }
    //////////////////////////////////////////////////////////////////////////
    for(int i =0;i<t0_p.size();i++)
    { //do    transpose  
        s0(i,0)=s1_p[i].num_x;s0(i,1)=s1_p[i].num_y;s0(i,2)=s1_p[i].num_z;
    }
   
      for(int i =0; i<index_1.size();i++)
    { //do    transpose  
        int fri =index_1[i].fir-1;int sed=index_1[i].sed-1;int thr=index_1[i].thr-1;
        double x1=s1_p[fri].num_x;double y1=s1_p[fri].num_y;double z1=s1_p[fri].num_z;
        double x2=s1_p[sed].num_x;double y2=s1_p[sed].num_y;double z2=s1_p[sed].num_z;
        double x3=s1_p[thr].num_x;double y3=s1_p[thr].num_y;double z3=s1_p[thr].num_z;     
        Point v4=caculate_v4(x1,y1,z1,x2,y2,z2,x3,y3,z3); 
        double x4=v4.num_x;double y4=v4.num_y;double z4=v4.num_z;
        s0(i+s0_p.size(),0)=x4;s0(i+s0_p.size(),1)=y4;s0(i+s0_p.size(),2)=z4;
    }
    //////////////////////////////////////////////////////////////////////////////// 
   for(int j =0;j<index_1.size();j++)
   {
        int fri =index_1[j].fir-1;int sed=index_1[j].sed-1;int thr=index_1[j].thr-1;
        double x1=t0_p[fri].num_x;double y1=t0_p[fri].num_y;double z1=t0_p[fri].num_z;
        double x2=t0_p[sed].num_x;double y2=t0_p[sed].num_y;double z2=t0_p[sed].num_z;
        double x3=t0_p[thr].num_x;double y3=t0_p[thr].num_y;double z3=t0_p[thr].num_z;     
        Point v4=caculate_v4(x1,y1,z1,x2,y2,z2,x3,y3,z3); 
        double x4=v4.num_x;double y4=v4.num_y;double z4=v4.num_z;
         mat0<<x2-x1,x3-x1,x4-x1,y2-y1,y3-y1,y4-y1,z2-z1,z3-z1,z4-z1;

         MatrixXd mat0_inverse=mat0.inverse();
         MatrixXd mat0inT=mat0_inverse.transpose(); 
         //*******************//
        double X1=s0_p[fri].num_x;double Y1=s0_p[fri].num_y;double Z1=s0_p[fri].num_z;
        double X2=s0_p[sed].num_x;double Y2=s0_p[sed].num_y;double Z2=s0_p[sed].num_z;
        double X3=s0_p[thr].num_x;double Y3=s0_p[thr].num_y;double Z3=s0_p[thr].num_z;     
        Point V4=caculate_v4(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3); 
        double X4=V4.num_x;double Y4=V4.num_y;double Z4=V4.num_z;
        mat1<<X2-X1,X3-X1,X4-X1,Y2-Y1,Y3-Y1,Y4-Y1,Z2-Z1,Z3-Z1,Z4-Z1;
        MatrixXd mat1_inverse=mat1.inverse();
         MatrixXd mat1inT=mat1_inverse.transpose(); 
         tripletLIst2.push_back(Triplet<double>(3*j,3*j,mat1inT(0,0)));tripletLIst2.push_back(Triplet<double>(3*j,1+3*j,mat1inT(0,1)));tripletLIst2.push_back(Triplet<double>(3*j,2+3*j,mat1inT(0,2)));
        tripletLIst2.push_back(Triplet<double>(1+3*j,3*j,mat1inT(1.0)));tripletLIst2.push_back(Triplet<double>(1+3*j,1+3*j,mat1inT(1,1)));tripletLIst2.push_back(Triplet<double>(1+3*j,2+3*j,mat1inT(1,2)));
        tripletLIst2.push_back(Triplet<double>(2+3*j,3*j,mat1inT(2,0)));tripletLIst2.push_back(Triplet<double>(2+3*j,1+3*j,mat1inT(2,1)));tripletLIst2.push_back(Triplet<double>(2+3*j,2+3*j,mat1inT(2,2)));
         //*******************//
    
        tripletLIst.push_back(Triplet<double>(3*j,3*j,mat0inT(0,0)));tripletLIst.push_back(Triplet<double>(3*j,1+3*j,mat0inT(0,1)));tripletLIst.push_back(Triplet<double>(3*j,2+3*j,mat0inT(0,2)));
        tripletLIst.push_back(Triplet<double>(1+3*j,3*j,mat0inT(1.0)));tripletLIst.push_back(Triplet<double>(1+3*j,1+3*j,mat0inT(1,1)));tripletLIst.push_back(Triplet<double>(1+3*j,2+3*j,mat0inT(1,2)));
        tripletLIst.push_back(Triplet<double>(2+3*j,3*j,mat0inT(2,0)));tripletLIst.push_back(Triplet<double>(2+3*j,1+3*j,mat0inT(2,1)));tripletLIst.push_back(Triplet<double>(2+3*j,2+3*j,mat0inT(2,2)));
        tripletLIst1.push_back(Triplet<double>(0+3*j,fri,-1.0));tripletLIst1.push_back(Triplet<double>(0+3*j,sed,1.0));
        tripletLIst1.push_back(Triplet<double>(1+3*j,fri,-1.0));tripletLIst1.push_back(Triplet<double>(1+3*j,thr,1.0));
        tripletLIst1.push_back(Triplet<double>(2+3*j,fri,-1.0));tripletLIst1.push_back(Triplet<double>(2+3*j,t0_p.size()+j,1.0));
        
   }    
  V.setFromTriplets(tripletLIst.begin(),tripletLIst.end());
   A.setFromTriplets(tripletLIst1.begin(),tripletLIst1.end());
  V1.setFromTriplets(tripletLIst2.begin(),tripletLIst2.end());

   S=V1*A*s0;
    //  cout<<V(0,0)<<endl;
    //  cout<<V(1,0)<<endl;

     D=V*A;
  
     SparseMatrix<double> D_transpose=D.transpose();
     SparseMatrix<double> DTD=D.transpose()*D;
     MatrixXd ls_b=D_transpose*S;
     
     SparseLU<SparseMatrix<double>> solver;
    solver.analyzePattern(DTD);
    solver.factorize(DTD);
  
    x=solver.solve(ls_b);
   
   
   for(int k=0;k<t0_p.size();k++)
   {
        my_t1.num_x=x(0+k,0);
        my_t1.num_y=x(0+k,1);
        my_t1.num_z=x(0+k,2);
        t1_p.push_back(my_t1);
   }
   return 0;
}
int main()
{
    double t_start = time(NULL);
     t_start=t_start*1000;

   int res1=read_obj("mesh/small/s0.obj",s0_p,index_1);
   int res2=read_obj("mesh/small/s1.obj",s1_p,index_2);
   int res3=read_obj("mesh/small/t0.obj",t0_p,index_3);

   if (res1<0) cout<<"re1 read worried"<<endl; 
   if (res2<0) cout<<"re2 read worried"<<endl; 
   if (res3<0) cout<<"re3 read worried"<<endl; 

    add_S2();
    caculate_target_2();

int res4= write_obj("mesh/small/t1.obj");
 if (res4<0) cout<<"re4 read worried"<<endl; 

    double t_end = time(NULL);
     t_end=t_end*1000;
cout<<difftime(t_end,t_start)<<"ms"<<endl;


   return 0;



}
