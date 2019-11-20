#include <iostream>
#include<math.h>
#include<fstream>
#include<cmath>
using namespace std;
double mass(double pm,double r0,double ya);
double fmass(double sigma,double x);
int main()
{
   double pmm,pm,r00,r0,dt,tmax,t,ya,va,vv,yy,rm,a,aux;
   const double c=6.671*(pow(10,-11));
   //gravitational constant
   cout<<"Please insert the following parameters:"<<endl;
   cout << "(Planet's mass)x(10^24) in kg=" << endl;
   cin>>pmm;
   cout<<"Planet's radius in km="<<endl;
   cin>>r00;
   cout<<"Time step in seconds=?"<<endl;
   cin>>dt;
   cout<<"Total time in seconds=?"<<endl;
   cin>>tmax;
   ofstream myfile("mass_distribution.dat");
   myfile<<"time(s)"<<"\t"<<"speed(km/s)"<<"\t"<<"position(km)"<<"\t"<<"M(r)(kg)"<<endl;
   pm=pmm*(pow(10,24));
   r0=r00*1000;
   ya=r0;
   va=0;
   for(double t=0; t<tmax; t=t+dt) //calculating mass distribution and writing to file
   {
       rm=mass(pm,r0,ya);
       aux=t;
       vv=va*3.6;
       yy=ya/1000;
       //conversion
       myfile<<aux<<"\t"<<vv<<"\t"<<yy<<"\t"<<rm<<endl;
       if(ya>=0)
       a=-c*rm/ya/ya;
       else a=c*rm/ya/ya;
       va=va+a*dt;
       ya=ya+va*dt+a*dt*dt/2;
   }
   myfile.close();
   cout<<"The file was created"<<endl;

   return 0;
}
double mass(double pm,double r0,double ya)
{
    int n;
    double d0,a,b,h,sp,si,ff,x,sigma,fa,fb;
    //d0=mean density
    const double pi=3.14159265;
    const double e=2.71828;
    const double d=1270*pow(10,8);
    d0=pm/(4*pi*(pow(r0,3))/3);
    a=0;
    b=abs(ya);
    n=300;
    sigma=r0/2;
    //sigma=dispersion
    h=(b-a)/n;
    sp=0;
    si=0;
    for(int i=0; i<=n; i++)
    {
        x=a+i*h;
        ff=fmass(sigma,x);
        if (i%2==0)
        sp=sp+4*ff;
        //i is even
        else si=si+2*ff;
    }
    fa=0;
    fb=fmass(sigma,b);
    double rm=(h/3)*(fa+si+sp+fb);
    //numerical integration using Simson's formula to determine planet mass at the given radius
    return rm;
}
double fmass(double sigma,double x) //Gaussian distribution function
{
    const double pi=3.14159265;
    const double d=1270*pow(10,8);
    const double e=2.71828;
    double ff=(4*pi*d/(sigma*(sqrt(2*pi))))*x*x*(pow(e,(-(x*x)/(2*sigma*sigma))));
    return ff;
}
