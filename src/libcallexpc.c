#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "aero2dcom.h"

int main(void)
{
  int i,prlv,iw1,iw2,iw3;
  char *mode;
  double max,min;

  FILE *fp;
  mode="C";
  strncpy(cpntctrl,"N",1);
  strncpy(cproctrl,"incom",5);
  strncpy(cenergy,"Y",1);
  strncpy(cvisheat,"N",1);
  strncpy(cturmod,"sa",2);
  strncpy(cwalltreat,"wf",2);
  strncpy(csolctrl,"SIMPLE",6);
  strncpy(cdiscret,"2upwind",7);
  strncpy(cinit,"N",1);
  strncpy(cstag,"N",1);
  maxs=2000;
  delta=1e-4;
  rau=7e-1;
  rap=3e-1;
  rae=7e-1;
  rat=3e-1;
  fd=1e-3;
  ifd=5;
  jp=75;
  c=0.5334;
  vfar=75.;
  aoa=4.;
  ta=263.15;
  tf=273.15;
  po=100000.;
  tvr=10.;
  ksi=0.;
  fp=fopen("cases/NACA0012.xyz","rt");
  fscanf(fp,"%d",&iwd);
  cxwd=(double *) malloc(iwd*sizeof(double));
  cywd=(double *) malloc(iwd*sizeof(double));
  for(i=0;i<iwd;i++)
    fscanf(fp,"%le %le",&cxwd[i],&cywd[i]);
  fscanf(fp,"%d",&iwu);
  cxwu=(double *) malloc(iwu*sizeof(double));
  cywu=(double *) malloc(iwu*sizeof(double));
  for(i=0;i<iwu;i++)
    fscanf(fp,"%le %le",&cxwu[i],&cywu[i]);
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
  allocarray_(mode,strlen(mode));
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
  printf("Maximum and minmum values of pressure are %le,%le.\n",max,min);
  deallocarray_(mode,strlen(mode));
  
  return 0;
}
