#include <iostream>
#include <cmath>
#include <fstream>


using namespace std;


void f(double* k, double* r, const double a,const double b,const double c){  


k[0] = a*(r[1] - r[0]);

k[1] = r[0]*(b-r[2]) - r[1];

k[2] = r[0]*r[1] - c*r[2];

}

int main(){
    double a = 10.0;
    double b = 28.0;
    double c = 8.0/3.0;
    ofstream out("data.txt");
    
    double r[3];
    r[0] = 1;
    r[1] = 1;
    r[2] = 1;

    
    double k1[3],k2[3],k3[3],k4[3],temp[3] = {0};
    double dt = 0.01;
    double N = 100/dt;
    
    out <<0 << "\t" << r[0] << "\t" << r[1] << "\t" << r[2] << endl;
    for(int i=1;i < N; i++){
    
        f(k1,r,a,b,c);
        
        int j;
        for(j=0; j < 3; j++){
            temp[j] = r[j] + (dt/2)* k1[j];
        }
        f(k2,temp,a,b,c);
        cout << k1[0] << k1[2] << endl;
        
        for(j=0; j < 3; j++){
            temp[j] = r[j] + (dt/2)* k2[j];
        }
        f(k3,temp,a,b,c);
        
        
        for(j=0; j < 3; j++){
            temp[j] = r[j] + dt* k3[j];
         }
         f(k4,temp,a,b,c);
         
        for(j=0; j < 3; j++){
            r[j]+= (dt)* (k1[j]/6 + k2[j]/3 + k3[j]/3 + k4[j]/6);
          }
          out <<i*dt << "\t" << r[0] << "\t" << r[1] << "\t" << r[2] << endl;
    }
    out.close();
    
            
    /*k2x[0] = x[0] + (dx/2)*k1x[0]/2;
    k2y[0] = y[0] + (dx/2)*k1y[0]/2;
    k2z[0] = z[0] + (dx/2)*k1z[0]/2;    
    
    k3x[0] = x[0] + (dx/2) * k2x[0]/2;
    k3y[0] = y[0] + (dx/2) * k2y[0]/2;
    k3z[0] = z[0] + (dx/2) * k2z[0]/2;
    
    k4x[0] = x[0] + dx * k3x[0];
    k4x[0] = x[0] + dx * k3x[0];
    k4x[0] = x[0] + dx * k3x[0];
    
    x[0] = (dx/6) * (k1x[0]/6 + k2x[0]/3 + k3x[0]/3 + k4x[0]/6);
    y[0] = (dx/6) * (k1y[0]/6 + k2y[0]/3 + k3y[0]/3 + k4y[0]/6);
    z[0] = (dx/6) * (k1z[0]/6 + k2z[0]/3 + k3z[0]/3 + k4z[0]/6);*/
    
    return 0;
    }

    