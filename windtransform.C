#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDEigen.h"
#include <TString.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
TMatrixD windtransform(double *xin, double *yin, double *zin,double *bxin, double *byin, double *bzin, int n, double *errin)
{
    double sx,sy,sz;
    TMatrixD amat(3,3);//3x3 dispersion matrix
    TMatrixD w(3,3);
    TMatrixD m(n,3);
    TMatrixD mt(n,3);
    amat[0][0]=amat[0][1]=amat[0][2]=amat[1][0]=amat[1][1]=amat[1][2]=amat[2][0]=amat[2][1]=amat[2][2]=0;
    sx=sy=sz=0.0;
    
    //Find centers of mass
    for(int i=0; i<n; i++)
    {
        sx+=xin[i]/n;
        sy+=yin[i]/n;
        sz+=zin[i]/n;
    }
    sx=sy=sz=0;
    //set up dispersion matrix
    for(int i=0; i<n; i++){
        amat[0][0]+=(xin[i]-sx-xin[0])*(xin[i]-sx-xin[0]);
        amat[0][1]+=(xin[i]-sx-xin[0])*(yin[i]-sy-yin[0]);
        amat[0][2]+=(xin[i]-sx-xin[0])*(zin[i]-sz-zin[0]);
        amat[1][0]+=(yin[i]-sy-yin[0])*(xin[i]-sx-xin[0]);
        amat[1][1]+=(yin[i]-sy-yin[0])*(yin[i]-sy-yin[0]);
        amat[1][2]+=(yin[i]-sy-yin[0])*(zin[i]-sz-zin[0]);
        amat[2][0]+=(zin[i]-sz-zin[0])*(xin[i]-sx-xin[0]);
        amat[2][1]+=(zin[i]-sz-zin[0])*(yin[i]-sy-yin[0]);
        amat[2][2]+=(zin[i]-sz-zin[0])*(zin[i]-sz-zin[0]);
    }
    //Find eigenvectors
    TMatrixDEigen b(amat);
    TVectorD w1(3);
    TVectorD wscaled(3);
    double wvals[3];
    w=b.GetEigenVectors();
    
    //Get the eigenvalues
    w1[0]=w[0][0];
    w1[1]=w[1][0];
    w1[2]=w[2][0];
    wscaled=amat*w1;
    sx=w1[0];
    if(sx==0) sx=w1[1];
    if(sx==0) sx=w1[2];
    wvals[0]=wscaled[1]/sx;
    
    w1[0]=w[0][1];
    w1[1]=w[1][1];
    w1[2]=w[2][1];
    wscaled=amat*w1;
    sx=w1[0];
    if(sx==0) sx=w1[1];
    if(sx==0) sx=w1[2];
    wvals[1]=wscaled[1]/sx;
    
    w1[0]=w[0][2];
    w1[1]=w[1][2];
    w1[2]=w[2][2];
    wscaled=amat*w1;
    sx=w1[0];
    if(sx==0) sx=w1[1];
    if(sx==0) sx=w1[2];
    wvals[2]=wscaled[1]/sx;

    //find the largest absolute value eigen vale
    double max,min;
    min=max=pow(wvals[0],2.0);
    for(int i=1; i<3; i++){
        if(pow(wvals[i],2.0)>max) max = pow(wvals[i],2.0);
        if(pow(wvals[i],2.0)<min) min = pow(wvals[i],2.0);
    }
    
    //Order eigenvectors by eigenvalue size. Largest one corresponds to best fit line, second largest to best fit plane, and third largest to plane normal vector.
    TMatrixD wtr(3,3);
    for(int i=0;i<3;i++){
        if(pow(wvals[i],2.0)==max){
            wtr[0][0]=w[0][i];
            wtr[1][0]=w[1][i];
            wtr[2][0]=w[2][i];
           // printf("max\n");
        }
        else if(pow(wvals[i],2.0)==min){
            wtr[0][2]=w[0][i];
            wtr[1][2]=w[1][i];
            wtr[2][2]=w[2][i];
           // printf("min\n");
        }
        else {
            wtr[0][1]=w[0][i];
            wtr[1][1]=w[1][i];
            wtr[2][1]=w[2][i];
            //printf("mid\n");
        }
    }
    
   /* max=pow(w[1][0],2.0);
    for(int i=1; i<3; i++){
        double m=pow(w[1][i],2.0);
        if(m>max) max=m;
    }
  //  printf("\na1: %f, a2: %f, a3: %f, max: %f, min %f\n",wvals[0],wvals[1],wvals[2],sqrt(max),min);
    for(int i=0;i<3;i++){
        if(pow(w[1][i],2.0)==max){
            int ind;
            for(int j=0;j<3;j++){
                if(pow(wtr[1][j],2.0)==max) ind=j;
            }
            sx=wtr[0][1];
            sy=wtr[1][1];
            sz=wtr[2][1];
            wtr[0][1]=w[0][i];
            wtr[1][1]=w[1][i];
            wtr[2][1]=w[2][i];
           // printf("i:  %d\n",ind);
            wtr[0][ind]=sx;
            wtr[1][ind]=sy;
            wtr[2][ind]=sz;
            break;
        }
    }*/
 
    //Now do the actual transformation
    sx=sy=sz=0;
    for(int i=0; i<n; i++){
        sx=xin[i];
        sy=yin[i];
        sz=zin[i];
        xin[i]=wtr[0][0]*sx+wtr[0][1]*sy+wtr[0][2]*sz;
        yin[i]=wtr[1][0]*sx+wtr[1][1]*sy+wtr[1][2]*sz;
        zin[i]=wtr[2][0]*sx+wtr[2][1]*sy+wtr[2][2]*sz;
        
        sx=bxin[i];
        sy=byin[i];
        sz=bzin[i];
        bxin[i]=wtr[0][0]*sx+wtr[0][1]*sy+wtr[0][2]*sz;
        byin[i]=wtr[1][0]*sx+wtr[1][1]*sy+wtr[1][2]*sz;
        bzin[i]=wtr[2][0]*sx+wtr[2][1]*sy+wtr[2][2]*sz;
    }
    
    //Transform the errors
    sx=errin[0];
    sy=errin[1];
    sz=errin[2];
    errin[0]=wtr[0][0]*sx+wtr[0][1]*sy+wtr[0][2]*sz;
    errin[1]=wtr[1][0]*sx+wtr[1][1]*sy+wtr[1][2]*sz;
    errin[2]=wtr[2][0]*sx+wtr[2][1]*sy+wtr[2][2]*sz;
    return wtr;
}

                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                
                                

