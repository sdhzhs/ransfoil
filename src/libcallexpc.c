#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "aero2dcom.h"

int main(void)
{
  int i,prlv,iw1,iw2,iw3,status;
  char *mode;
  double max,min;

  FILE *fp;
  mode="C";
  maxs=2000;
  delta=1e-4;
  rau=7e-1;
  rap=3e-1;
  rae=7e-1;
  rat=3e-1;
  fd=1e-3;
  ifd=5;
  jp=75;
  lfar=10.0;
  c=0.5334;
  vfar=75.;
  aoa=4.;
  ta=263.15;
  tf=273.15;
  po=100000.;
  tvr=10.;
  ksi=0.;
  turmodflag=SA;
  proctrlflag=INCOM;
  walltreatflag=WF;
  fstypeflag=ALLFIX;
  solctrlflag=SIMPLE;
  discretflag=SUP;
  linsolflag=SORGS;
  tmptypeflag=FIXED;
  gtypeflag=CTYPE;
  initflag=INITBC;
  pntctrlflag=false;
  energyflag=true;
  visheatflag=false;
  stagflag=false;
  matairflag=true;
  wallfunutype=PARVEL;
  wallfunktype=LOGLAW;
  sstlowre=false;
  fp=fopen("../cases/NACA0012.xyz","rt");
  status=fscanf(fp,"%d",&iwd);
  cxwd=(double *) malloc(iwd*sizeof(double));
  cywd=(double *) malloc(iwd*sizeof(double));
  for(i=0;i<iwd;i++)
    status=fscanf(fp,"%le %le",&cxwd[i],&cywd[i]);
  status=fscanf(fp,"%d",&iwu);
  cxwu=(double *) malloc(iwu*sizeof(double));
  cywu=(double *) malloc(iwu*sizeof(double));
  for(i=0;i<iwu;i++)
    status=fscanf(fp,"%le %le",&cxwu[i],&cywu[i]);
  fclose(fp);
  printf("Read airfoil coordinates completed!\n");
  iw1=(iwd>iwu)?iwd:iwu;
  iw2=iw1+iwd-1;
  iw3=iw2+iwu-1;
  ip=iw3+iw1-1;
  ic=ip-1;
  jc=jp-1;
  ib1=iw1;
  ib2=iw3-1;
  allocarray_();
  prlv=1;
  aero2d_(mode,&prlv,"\0",strlen(mode),1);
  printf("%le,%le,%le,%le\n",cl,cd,cf,cm);
  printf("%le,%le\n",xpc,ypc);
  max=0;
  for(i=0;i<ic*jc;i++)
    if(cp[i]>max) max=cp[i];
  min=1e+5;
  for(i=0;i<ic*jc;i++)
    if(cp[i]<min) min=cp[i];
  printf("Maximum and minimum values of pressure are %le,%le.\n",max,min);
  deallocarray_();
  free(cxwd);
  free(cywd);
  free(cxwu);
  free(cywu);
  
  return 0;
}
