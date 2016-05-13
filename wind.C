#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TSystem.h"
#include "TRandom2.h"
#include <stdio.h>
#include <stdlib.h>
#include "wind_qspline.C"
#include "windspline.C"
#include "windtransform.C"
#include "math.h"

void wind();
void wind()
{
    //this might be needed for root 5
    //gSystem->CompileMacro("wind_qspline.C");
   // gSystem->CompileMacro("windspline.C");
   // gSystem->CompileMacro("windtransform.C");
    
    double clight=0.00002999792458;
   
   //TString ipfname = "TTest_AllAngXY_Fixed.root";
    TString ipfname = "TTest_AllAngXY_Real.root";
    //TString ipfname = "TTest_Centered_Fixed.root";
    //TString ipfname = "TTest_Centered_Real.root";
    TString bfname = "HELIXBmapFixed.root";

    TFile *ipf = new TFile(ipfname.Data()); //data file
    TFile *bf = new TFile(bfname.Data()); //B field file
    
    TH2D *bfrho = (TH2D*)bf->Get("hbrho"); //rho B field
    TH2D *bfx = (TH2D*)bf->Get("hbx"); //x b field
    if(!ipf || ! bf){
        printf("Couldn't open %s!!! Bailing out\n", ipfname.Data());
        exit(-1);
    }

    // Variables to read from the tree
    double *inG = new double;   //gamma
    double inB;   //beta
    double inA;   //nucleus mass
    double inQ;   //nucleus charge
    double *inE = new double;   //energy
    double inR;   //rigidity
    double inP;   //momentum
    double inPLC; //pathlength cosine

    double inLoc[3]; //starting position
    double inPhi;    //starting azimuth
    double inTheta;  //starting polar
    

    //Track truth values
    const int fDCTLayMax = 75*2; //bigger than actual number to cover rare multiple crossings
    int    inTrkCNT;                    //actual number of layer hits
    int    inTrkLay[fDCTLayMax];        //ID of the layer we are in
    double inTrkX[fDCTLayMax];
    double inTrkY[fDCTLayMax];
    double inTrkZ[fDCTLayMax];
    double inTrkT[fDCTLayMax];

    // make the tree and associate variable names with the variables
    TTree *t = (TTree*)ipf->Get("t");
    t->SetBranchAddress("inE",         &inE);
    t->SetBranchAddress("inQ",         &inQ);
    t->SetBranchAddress("inG",         &inG);
    t->SetBranchAddress("inB",         &inB);
    t->SetBranchAddress("inP",         &inP);
    t->SetBranchAddress("inR",         &inR);
    t->SetBranchAddress("inA",         &inA);
    t->SetBranchAddress("inLoc",       inLoc);
    t->SetBranchAddress("inTheta",     &inTheta);
    t->SetBranchAddress("inPhi",       &inPhi);

    t->SetBranchAddress("inTrkCNT",    &inTrkCNT);
    t->SetBranchAddress("inTrkLay",    inTrkLay);
    t->SetBranchAddress("inTrkX",      inTrkX);
    t->SetBranchAddress("inTrkY",      inTrkY);
    t->SetBranchAddress("inTrkZ",      inTrkZ);
    //t->SetBranchAddress("inTrkT",      inTrkT);

    bool displayOne = true;
    bool doFit      = true;

    // loop over all the events in the tree
    int CNT = t->GetEntries();


    printf("CNT: %d, \n",CNT);
    
    //TH2D *h2 = new TH2D("h2","",1000,-.5,.5,200,0,1000);//-4000000,10000000);
    TH1D *h2 = new TH1D("h2","",1000,-.5,.5);
    h2->GetXaxis()->SetTitle("(R-R_{in})/R_{in}");
    h2->SetTitle(ipfname.Data());
      //h12->GetYaxis()->SetTitle("Detector Z (mm)");
      //h12->GetYaxis()->SetTitleOffset(1.4);

      double xq,yq,zq,bxq,byq,bzq;
      double bx[fDCTLayMax],by[fDCTLayMax],bz[fDCTLayMax];
      double x[fDCTLayMax],y[fDCTLayMax],z[fDCTLayMax];
      double a[fDCTLayMax];
      double b[fDCTLayMax];
      double yfit[fDCTLayMax];
      double yp[fDCTLayMax];
      double yppp[fDCTLayMax];
      double ypppp[fDCTLayMax];
      double zfit[fDCTLayMax];
      double zp[fDCTLayMax];
      double zppp[fDCTLayMax];
      double zpppp[fDCTLayMax];
      double ex,ey;
      double chisq=0;
      double chisq2=0;
      double err[3];
      TMatrixD lsqerr(5,5);
      err[0]=0.0;//500 micron error in x
      err[1]=0.07;//70 micron error in y
      err[2]=0.5;//0 micron error in z
    double coef[5];

    for (int i=0;i<CNT;i++) {
              t->GetEntry(i);
        if(inTrkCNT>74)
        {
            //at this point all the variables above are available for analysis!!
            

            if(i%1==0){
                printf("n: %d \n",i);
                TRandom2 r;
                double xm=0;
                double ym=0;
                double zm=0;
                
                inTrkCNT=inTrkCNT-0;
                r.SetSeed(i+1);
                for(int g =0;g<inTrkCNT;g++) {
                    //spread the x and y measurements
                    inTrkY[g]+=r.Gaus(0,0.07);
                    inTrkX[g]+=r.Gaus(0,0.5);
                    
                    xm+=inTrkX[g]/inTrkCNT;
                    ym+=inTrkY[g]/inTrkCNT;
                    zm+=inTrkZ[g]/inTrkCNT;
                }
                TMatrixD w(3,3);
                double rho;
                //get the B fields
                for (int q=0;q<inTrkCNT;q++)
                {
                    rho=sqrt(inTrkY[q]*inTrkY[q]+inTrkZ[q]*inTrkZ[q]);
                    if(ipfname.Contains("Real")){
                        bzq=(bfrho->Interpolate(inTrkX[q],rho))*inTrkZ[q]/rho;
                        byq=(bfrho->Interpolate(inTrkX[q],rho))*inTrkY[q]/rho;
                        bxq=bfx->Interpolate(inTrkX[q],rho);
                    }
                    else{
                        bzq=0;
                        byq=0;
                        bxq=1;
                    }
                     bx[q]=bxq;
                     by[q]=byq;
                     bz[q]=bzq;
                }
                //transform into the algorithim coordinate system where dy/dx and dz/dx will be small and dx direction is the best fit line to the track
                w= windtransform(inTrkX,inTrkY,inTrkZ,bx,by,bz,inTrkCNT,err);
                
                //estimate dy/dx and dz/dx with a cubic spline fit
                windspline(inTrkY,inTrkX,inTrkCNT,0,yp,a,yppp,ypppp);
                windspline(inTrkZ,inTrkX,inTrkCNT,0,zp,b,zppp,zpppp);

        
                //Set up a=p*d^2y/dx^2 and b=p*d^2z/dx^2 as described in the Wind paper
                for(int n=0;n<inTrkCNT;n++)
                {
                    a[n] = sqrt(1+yp[n]*yp[n]+zp[n]*zp[n])*(bx[n]*zp[n]+by[n]*zp[n]*yp[n]-bz[n]*(1+yp[n]*yp[n]));
                    b[n] = sqrt(1+yp[n]*yp[n]+zp[n]*zp[n])*(-bx[n]*yp[n]-bz[n]*zp[n]*yp[n]+by[n]*(1+zp[n]*zp[n]));
                }
                //Transform the center of masses into the algorithm coordinate system
                bxq=xm;
                byq=ym;
                bzq=zm;
                xm=w[0][0]*bxq+w[0][1]*byq+w[0][2]*bzq;
                ym=w[1][0]*bxq+w[1][1]*byq+w[1][2]*bzq;
                zm=w[2][0]*bxq+w[2][1]*byq+w[2][2]*bzq;
                
                //we will also have to calculate the center of mass of the integrated spline fits to a and b
                double syf,szf,lsqerr[25];
                for(int it=0;it<5;it++){
                    windspline(a,inTrkX,inTrkCNT,1,ypppp,yppp,yp,yfit);
                    windspline(b,inTrkX,inTrkCNT,1,zpppp,zppp,zp,zfit);
                    wind_qspline(inTrkX,inTrkY,inTrkZ,yfit,zfit,inTrkCNT,coef,lsqerr);
                    //This is the problem, if I don't print here the code messes up and I get runnaway rigidities
                    printf("sdf: %f",yfit[4]);
                    syf=szf=0;
                    for(int g=0;g<inTrkCNT;g++){
                        yp[g]=coef[1]+coef[2]*yp[g];
                        zp[g]=coef[4]+coef[2]*zp[g];
                        a[g] = sqrt(1+yp[g]*yp[g]+zp[g]*zp[g])*(bx[g]*zp[g]+by[g]*zp[g]*yp[g]-bz[g]*(1+yp[g]*yp[g]));
                        b[g] = sqrt(1+yp[g]*yp[g]+zp[g]*zp[g])*(-bx[g]*yp[g]-bz[g]*zp[g]*yp[g]+by[g]*(1+zp[g]*zp[g]));
                        syf+=yfit[g]/inTrkCNT;
                        szf+=zfit[g]/inTrkCNT;
                    }
                }
                double dy,dz;
                double da1,db1;
                dy=dz=da1=db1=0.0;
                chisq2=0;
                chisq=0;
                for(int n=0;n<inTrkCNT;n++)
                {
                    da1=coef[1]*inTrkX[n]*sqrt(lsqerr[1+5*1]/(coef[1]*coef[1])+.25/(inTrkX[n]*inTrkX[n]));
                    db1=coef[4]*inTrkX[n]*sqrt(lsqerr[4+4*5]/(coef[4]*coef[4])+.25/(inTrkX[n]*inTrkX[n]));
                    dy=lsqerr[0]+da1*da1/pow(coef[1]*inTrkX[n],2.0)+pow(sqrt(lsqerr[2+2*5])*yfit[n],2.0);
                    dz=lsqerr[3+3*5]+db1*db1/pow(coef[4]*inTrkX[n],2.0)+pow(sqrt(lsqerr[2+2*5])*zfit[n],2.0);
                    chisq+=(dy+dz)/inTrkCNT;//what the least squares fit thinks its residuals are
                    
                    dy=(coef[0]+coef[1]*(inTrkX[n]-xm)+coef[2]*(yfit[n]-syf)-(inTrkY[n]-ym))/(inTrkY[n]-ym);
                    dz=(coef[3]+coef[4]*(inTrkX[n]-xm)+coef[2]*(zfit[n]-szf)-(inTrkZ[n]-zm))/(inTrkZ[n]-zm);
                    
                    chisq2+=(dy*dy+dz*dz)/inTrkCNT;//sum of squares of actual residuals
                }
                h2->Fill((sqrt(clight/coef[2]*clight/coef[2])-inP/inQ)/(inP/inQ)); //MDR=clight/sqrt(lsqerr[2+2*5]*10
                //h2->Fill(clight/sqrt(lsqerr[2+2*5])/100);
                if(1){
                printf("RMS= %.20f\n", chisq2);
                printf("rigidity resolution: %f \n",(clight/coef[2]-inP/inQ));
                printf("rigidty: %f, fit rigidity %f , P %f\n\n", inP/inQ, sqrt(clight/coef[2]*clight/coef[2]), inP);
                printf("a1: %f, a2: %f, q/p: %.15f\n",coef[0],coef[1],coef[2]);
                printf("b1: %f, b2: %f, q/p: %.15f\n\n",coef[3],coef[4],coef[2]);
                printf("11: %f, 22: %.15f, 33: %f, 44: %.15f, 55: %f, \n\n",lsqerr[0],lsqerr[1+1*5],clight/sqrt(lsqerr[2+2*5]),lsqerr[3+3*5],lsqerr[4+4*5]);
                printf("phi: %f, theta: %f\n",inPhi,inTheta);
                    
                    
        }

    }//track cut 1
  }//Track cut 2
    }//tree loop

  TCanvas *c = new TCanvas("c","",0,0,1.2*500,1.2*600);
  h2->Draw();
    
}//function
