#include <stdbool.h>

extern int iwd,iwu,iw0,iw,ip,jp,ic,jc,ib1,ib2,ifd,maxs;
extern double fd,fb,eb,lfar,delta,rau,rap,rae,rat,c,vfar,aoa,ta,tf,po,itur,tvr,ksi,cl,cd,cf,cm,xpc,ypc;
extern char cfilename[64],cdir[64],cmatfile[64];
extern double *cxwd,*cywd,*cxwu,*cywu,*cxwp0,*cywp0;
extern double *cxw,*cyw,*csw,*cyplus,*cystar,*chcv,*cax,*cay;
extern double *cxg,*cyg,*cxc,*cyc,*crho,*cmu,*cp,*cvx,*cvy,*ct,*ctn,*ctk,*cte,*ctw,*cmut;
extern int turmodflag,proctrlflag,walltreatflag,fstypeflag,solctrlflag,discretflag,denfaceflag,linsolflag,tmptypeflag,gtypeflag,initflag,
wallfunutype,wallfunktype;
const int INV=0,LAM=1,SA=2,KE=3,SST=4,INCOM=0,COM=1,WF=0,LR=1,ALLFIX=0,VINPOUT=1,SIMPLE=0,SIMPLEC=1,FUP=0,SUP=1,QUICK=2,TVD=3,CDS=4,FROMM=5,
SORGS=0,PBICG=1,FIXED=0,FLUX=1,CTYPE=0,OTYPE=1,INITBC=0,INITFILE=1,INITMEM=2,ORIVEL=0,PARVEL=1,LOGLAW=0,GENLAW=1;
extern bool pntctrlflag,energyflag,visheatflag,stagflag,matairflag,sstlowre;
void aero2d_(char *mode,int *prlv,char *scptname,size_t nmod,size_t nscpt);
void allocarray_();
void deallocarray_();