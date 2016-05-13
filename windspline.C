#include <TString.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
void windspline( double aspl[], double xspl[], int n, int integrate, double *ap, double *app, double *Sa, double *SSa){
    //Taken from H. Wind paper.
    double p[75*2],q[75*2],d[75*2],ta,tb,r;
    int i;
    
    //first find p and q
    p[0]=1;
    q[0]=0;
    d[0]=xspl[1]-xspl[0];
    for(i=0;i<n-2;i++){
        d[i+1]=xspl[i+2]-xspl[i+1];
        ta=d[i]/(2*(d[i]+d[i+1]));
        tb=d[i+1]/(2*(d[i]+d[i+1]));
        r=3*((aspl[i+2]-aspl[i+1])/d[i+1]-(aspl[i+1]-aspl[i])/d[i])/(d[i]+d[i+1]);
        
        p[i+1]=-tb/(1+ta*p[i]);
        q[i+1]=(r-ta*q[i])/(1+ta*p[i]);
        
    }
    
    //next find app
    app[n-1]=q[n-1]/(1-p[n-1]);
    app[n-2]=app[n-1];
    for(i=n-3;i>0;i--){
        app[i]=p[i]*app[i+2]+q[i];
    }
    app[0]=app[1];
    
    //next find ap, Sa, Saa
    Sa[0]=0;
    SSa[0]=0;
    for(i=0;i<n-1;i++){
        ap[i+1]=(aspl[i+1]-aspl[i])/d[i]-d[i]*(app[i]/3.0+app[i+1]/6.0-app[i]-(app[i+1]-app[i])/2.0);
        
        if(integrate){
            Sa[i+1]=Sa[i]+d[i]*(aspl[i+1]+aspl[i])/2.0-d[i]*d[i]*d[i]*(app[i]+app[i+1])/24.0;
            SSa[i+1]=SSa[i]+Sa[i]*d[i]+d[i]*d[i]*(2*aspl[i]+aspl[i+1])/6.0-d[i]*d[i]*d[i]*d[i]*(0.8*app[i]+0.7*app[i+1])/36.0;
        }
    }
    ap[0]=ap[1];
}