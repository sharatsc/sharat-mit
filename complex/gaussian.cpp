#include "mex.h"
#include "matrix.h"
#include <math.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{


const int  *vettore_im;        /* vettore con componenti*/
int imx,imy;                   // componenti
double *imd;                   // ingresso
double sigma;                  // parametro sigma
int padding;                   // padding di zeri
double *outd;                  // uscita 
double *tempd;                 // immagine temporanea
int tempx,tempy;               // dimensione immagine paddata con zero
int ii,jj;                     // contatori
double *in_scand,*out_scand;   // puntatori che scandiscono immagine
double q,qq,qqq,b0,b1,b2,b3,B; // parametri legati a sigma
double *tempd1,*tempd2;        // matrici temporanee
double bc;                     // variabile ausiliaria
int riga;                      // boolean
          

imd=mxGetPr(prhs[0]);
vettore_im=mxGetDimensions(prhs[0]);
imx=vettore_im[0];
imy=vettore_im[1];
sigma=mxGetScalar(prhs[1]);

q = 1.31564 * (sqrt(1 + 0.490811 * sigma*sigma) - 1);
qq = q*q;
qqq = qq*q;
b0 = 1.0/(1.57825 + 2.44413*q + 1.4281*qq + 0.422205*qqq);
b1 = (2.44413*q + 2.85619*qq + 1.26661*qqq)*b0;
b2 = (-1.4281*qq - 1.26661*qqq)*b0;
b3 = 0.422205*qqq*b0;
B = 1.0 - (b1 + b2 + b3);

if ((imx>3)&&(imy>3))
{
if (nrhs==3)                  // 3 ingressi: padding
 {
 
  padding=mxGetScalar(prhs[2]);
   
  tempx=imx+2*padding;
  tempy=imy+2*padding;
  tempd=(double*)mxCalloc(tempx*tempy,sizeof(double));
  tempd1=(double*)mxCalloc(tempx*tempy,sizeof(double));
  tempd2=(double*)mxCalloc(tempx*tempy,sizeof(double));


  // sull'immagine ingrandita ricopio l'immagine originale 
       
  for (jj=0;jj<imy;jj++)
  { out_scand=tempd+tempx*(jj+padding)+padding;
    in_scand=imd+imx*jj;
    for (ii=0;ii<imx;ii++)
      {*out_scand=*in_scand;
        out_scand++;
        in_scand++;
      }
  }
  
 }
else 
 {
   

   tempx=imx;                // 2 soli ingressi: nessun padding
   tempy=imy;
   tempd=(double*)mxCalloc(tempx*tempy,sizeof(double));
   tempd1=(double*)mxCalloc(tempx*tempy,sizeof(double));
   tempd2=(double*)mxCalloc(tempx*tempy,sizeof(double));

   for (jj=0;jj<imy;jj++)
     { out_scand=tempd+tempx*jj;
       in_scand=imd+imx*jj;
       for (ii=0;ii<imx;ii++)
         {*out_scand=*in_scand;
           out_scand++;
           in_scand++;
         }
     }
  
 }


// ora tempd punta ad immagine 2D (paddata o meno con zeri)  
for (jj=0;jj<tempy;jj++)
 {
   bc=*(tempd+tempx*jj);
   *(tempd1+tempx*jj+0)=B**(tempd+jj*tempx+0)+b1*bc+b2*bc+b3*bc;
   *(tempd1+tempx*jj+1)=B**(tempd+jj*tempx+1)+b1*bc+b2*bc+b3*bc;
   *(tempd1+tempx*jj+2)=B**(tempd+jj*tempx+2)+b1*bc+b2*bc+b3*bc;
    
   for (ii=3;ii<tempx;ii++)
     {
       *(tempd1+tempx*jj+ii)=B**(tempd+tempx*jj+ii)+b1**(tempd1+tempx*jj+ii-1)+b2**(tempd1+tempx*jj+ii-2)+b3**(tempd1+tempx*jj+ii-3);
     }
   
   bc=*(tempd1+tempx*jj+tempx-1);
   *(tempd2+tempx*jj+tempx-1-0)=B**(tempd1+tempx*jj+tempx-1-0)+b1*bc+b2*bc+b3*bc;
   *(tempd2+tempx*jj+tempx-1-1)=B**(tempd1+tempx*jj+tempx-1-1)+b1*bc+b2*bc+b3*bc;
   *(tempd2+tempx*jj+tempx-1-2)=B**(tempd1+tempx*jj+tempx-1-2)+b1*bc+b2*bc+b3*bc;

   for (ii=tempx-4;ii>=0;ii--)
     {
       *(tempd2+tempx*jj+ii)=B**(tempd1+tempx*jj+ii)+b1**(tempd2+tempx*jj+ii+1)+b2**(tempd2+tempx*jj+ii+2)+b3**(tempd2+tempx*jj+ii+3);  
     }
  
 }


for (ii=0;ii<tempx;ii++)
 {
   bc=*(tempd2+ii);
   *(tempd1+tempx*0+ii)=B**(tempd2+tempx*0+ii)+b1*bc+b2*bc+b3*bc;
   *(tempd1+tempx*1+ii)=B**(tempd2+tempx*1+ii)+b1*bc+b2*bc+b3*bc;
   *(tempd1+tempx*2+ii)=B**(tempd2+tempx*2+ii)+b1*bc+b2*bc+b3*bc;

   for (jj=3;jj<tempy;jj++)
    {
      *(tempd1+tempx*jj+ii)=B**(tempd2+tempx*jj+ii)+b1**(tempd1+tempx*(jj-1)+ii)+b2**(tempd1+tempx*(jj-2)+ii)+b3**(tempd1+tempx*(jj-3)+ii);
    }
   
  
   bc=*(tempd1+tempx*(tempy-1)+ii);
   *(tempd2+tempx*(tempy-1-0)+ii)=B**(tempd1+tempx*(tempy-1-0)+ii)+b1*bc+b2*bc+b3*bc;
   *(tempd2+tempx*(tempy-1-1)+ii)=B**(tempd1+tempx*(tempy-1-1)+ii)+b1*bc+b2*bc+b3*bc;
   *(tempd2+tempx*(tempy-1-2)+ii)=B**(tempd1+tempx*(tempy-1-2)+ii)+b1*bc+b2*bc+b3*bc;
  

   for (jj=tempy-4;jj>=0;jj--)
    {
      *(tempd2+tempx*jj+ii)=B**(tempd1+tempx*jj+ii)+b1**(tempd2+tempx*(jj+1)+ii)+b2**(tempd2+tempx*(jj+2)+ii)+b3**(tempd2+tempx*(jj+3)+ii);
    }



 }


plhs[0] = mxCreateDoubleMatrix(imx,imy,mxREAL);
outd=mxGetPr(plhs[0]);

if (nrhs==3)
{
  for (jj=0;jj<imy;jj++)
    { in_scand=tempd2+tempx*(jj+padding)+padding;
      out_scand=outd+imx*jj;
      for (ii=0;ii<imx;ii++)
        {*out_scand=*in_scand;
          out_scand++;
          in_scand++;
        }
    }
}

else 
 { 
   for (jj=0;jj<imy;jj++)
    { in_scand=tempd2+tempx*jj;
      out_scand=outd+imx*jj;
      for (ii=0;ii<imx;ii++)
        {*out_scand=*in_scand;
          out_scand++;
          in_scand++;
        }
    }

 }

mxFree(tempd);
mxFree(tempd1);
mxFree(tempd2);
return;


}


if ((imx==1)||(imy==1))
{
  if (imx==1)
    {imx=imy;
     riga=0;
    }
   else
    { riga=1;
    }


  if (nrhs==3)
  {

    padding=mxGetScalar(prhs[2]);
   
    tempx=imx+2*padding;
    
    tempd=(double*)mxCalloc(tempx,sizeof(double));
    tempd1=(double*)mxCalloc(tempx,sizeof(double));
    tempd2=(double*)mxCalloc(tempx,sizeof(double));

    for (ii=0;ii<imx;ii++)
    {
      *(tempd+ii+padding)=*(imd+ii);
    }
  } 
  else
  {

    tempx=imx;                // 2 soli ingressi: nessun padding
    
    tempd=(double*)mxCalloc(tempx,sizeof(double));
    tempd1=(double*)mxCalloc(tempx,sizeof(double));
    tempd2=(double*)mxCalloc(tempx,sizeof(double));

    for (ii=0;ii<imx;ii++)
    {
      *(tempd+ii)=*(imd+ii);
    }


  }

   bc=*(tempd);
   *(tempd1+0)=B**(tempd+0)+b1*bc+b2*bc+b3*bc;
   *(tempd1+1)=B**(tempd+1)+b1*bc+b2*bc+b3*bc;
   *(tempd1+2)=B**(tempd+2)+b1*bc+b2*bc+b3*bc;
    
   for (ii=3;ii<tempx;ii++)
     {
       *(tempd1+ii)=B**(tempd+ii)+b1**(tempd1+ii-1)+b2**(tempd1+ii-2)+b3**(tempd1+ii-3);
     }
   
   bc=*(tempd1+tempx-1);
   *(tempd2+tempx-1-0)=B**(tempd1+tempx-1-0)+b1*bc+b2*bc+b3*bc;
   *(tempd2+tempx-1-1)=B**(tempd1+tempx-1-1)+b1*bc+b2*bc+b3*bc;
   *(tempd2+tempx-1-2)=B**(tempd1+tempx-1-2)+b1*bc+b2*bc+b3*bc;

   for (ii=tempx-4;ii>=0;ii--)
     {
       *(tempd2+ii)=B**(tempd1+ii)+b1**(tempd2+ii+1)+b2**(tempd2+ii+2)+b3**(tempd2+ii+3);  
     }

   

    if (riga==1)
     plhs[0] = mxCreateDoubleMatrix(imx,1,mxREAL);
    else
     plhs[0] = mxCreateDoubleMatrix(1,imx,mxREAL);


    outd=mxGetPr(plhs[0]);

    if (nrhs==3)
    {
       for (ii=0;ii<imx;ii++)
        {
          *(outd+ii)=*(tempd2+ii+padding);
        }
    }
    else
    {
       for (ii=0;ii<imx;ii++)
        {
          *(outd+ii)=*(tempd2+ii);
        }

    }
   
  

    mxFree(tempd);
    mxFree(tempd1);
    mxFree(tempd2);
    return;
 


  


}


} // fine programma	 
