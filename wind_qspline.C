#include <iostream>
#include "TMatrixD.h"
#include "TVectorD.h"
void wind_qspline(double xalg[], double yalg[], double zalg[], double yc[], double zc[], int n, double *c,double *erralg) {
    /*This routine finds the coeeficients of the least squares track fit from cubic spline fits y and z.
    Taken from AN IMPROVEMENT TO ITERATIVE TRACKING FOR MOMENTUM DETERMINATION
, H. Wind, NUCLEAR INSTRUMENTS AND METHODS 153 (1978) 195-197*/

    double sx,sx2,sxy,sxz,sz,sz2,sy,sy2,syyc,szzc,sxyc,sxzc,szc,syc,szc2,syc2;
    sx=sx2=sxy=sxz=sz=sz2=sy=sy2=syyc=szzc=sxyc=sxzc=szc=syc=szc2=syc2=0.0;
    double t1,t2,t3,t4,t5,t7,t8,t10,t11,det,num;
    t1=t2=t3=t4=t5=t7=t8=t10=t11=det=num=0.0;
    TMatrixD e(2*n,5);
    double mx,my,mz,myc,mzc;
    mx=my=mz=myc=mzc=0.0;
    
    //need the center of masses to do the least squares fit
    for (int i=0; i<n ;i++) {
        mx +=xalg[i]/n;
        my +=yalg[i]/n;
        mz +=zalg[i]/n;
        myc +=yc[i]/n;
        mzc +=zc[i]/n;
    }
    for (int i=0; i<n ;i++) {
        sx +=(xalg[i]-mx);
        sx2 += (xalg[i]-mx)*(xalg[i]-mx);
        sxyc+= (xalg[i]-mx)*(yc[i]-myc);
        sxzc +=(xalg[i]-mx)*(zc[i]-mzc);
        szc += (zc[i]-mzc);
        szc2 += (zc[i]-mzc)*(zc[i]-mzc);
        syc += (yc[i]-myc);
        syc2 += (yc[i]-myc)*(yc[i]-myc);
        syyc +=(yalg[i]-my)*(yc[i]-myc);
        szzc +=(zalg[i]-mz)*(zc[i]-mzc);
        sxy += (xalg[i]-mx)*(yalg[i]-my);
        sxz+= (xalg[i]-mx)*(zalg[i]-mz);
        sz +=(zalg[i]-mz);
        sy += (yalg[i]-my);
        
        //This is for the error matrix
        e[i][0]=1;
        e[i][1]=xalg[i]-mx;
        e[i][2]=yc[i]-myc;
        e[i][3]=0;
        e[i][4]=0;
        
        e[n+i][0]=0;
        e[n+i][1]=0;
        e[n+i][3]=1;
        e[n+i][4]=xalg[i]-mx;
        e[n+i][2]=zc[i]-mzc;
    }
    //Belowed are as defined in the followup to the initial Wind paper
    t1=sx2*sy-sx*sxy;
    t2=syc*sx2-sxyc*sx;
    t3=n*sx2-sx*sx;
    t4=n*sxy-sx*sy;
    t5=n*sxyc-syc*sx;
    t7=sx2*sz-sx*sxz;
    t8=sx2*szc-sx*sxzc;
    t10=n*sxz-sx*sz;
    t11=n*sxzc-sx*szc;
    det= t3*(syc2+szc2)-t2*syc-t5*sxyc-t8*sz-t11*sxzc;
    num=t3*(syyc+szzc)-t1*syc-t4*sxyc-t7*sz-t10*sxzc;

    if(det!=0){
    c[2]=num/det;
    }
    else c[2]=0;
    if(t3!=0){
    c[0]=(t1-t2*c[2])/t3;
    c[1]=(t4-t5*c[2])/t3;
    c[3]=(t7-t8*c[2])/t3;
    c[4]=(t10-t11*c[2])/t3;
    }
    else{
        c[0]=0;
        c[1]=0;
        c[3]=0;
        c[4]=0;
    }
    
    //Caclulate the error matrix
    double cov=0;
    for(int i=0;i<n;i++){
        cov+=(pow(c[0]+c[1]*(xalg[i]-mx)+c[2]*(yc[i]-myc)-(yalg[i]-my),2.0)+pow(c[3]+c[4]*(xalg[i]-mx)+c[2]*(zc[i]-mzc)-(zalg[i]-mz),2.0))/(n-5);
    }
    //cov=(.5*.5+.07*.07); //Use this for error free tracks
    
    TMatrixD et(2*n,5);
    TMatrixD esq(5,5);
    et=e;
    et.T();
    esq=et*e;
    esq=esq.Invert();
    
    erralg[0]=esq[0][0]*cov;
    erralg[1]=esq[0][1]*cov;
    erralg[2]=esq[0][2]*cov;
    erralg[3]=esq[0][3]*cov;
    erralg[4]=esq[0][4]*cov;
    erralg[5]=esq[1][0]*cov;
    erralg[6]=esq[1][1]*cov;
    erralg[7]=esq[1][2]*cov;
    erralg[8]=esq[1][3]*cov;
    erralg[9]=esq[1][4]*cov;
    erralg[10]=esq[2][0]*cov;
    erralg[11]=esq[2][1]*cov;
    erralg[12]=esq[2][2]*cov;
    erralg[13]=esq[2][3]*cov;
    erralg[14]=esq[2][4]*cov;
    erralg[15]=esq[3][0]*cov;
    erralg[16]=esq[3][1]*cov;
    erralg[17]=esq[3][2]*cov;
    erralg[18]=esq[3][3]*cov;
    erralg[19]=esq[3][4]*cov;
    erralg[20]=esq[4][0]*cov;
    erralg[21]=esq[4][1]*cov;
    erralg[22]=esq[4][2]*cov;
    erralg[23]=esq[4][3]*cov;
    erralg[24]=esq[4][4]*cov;
}
