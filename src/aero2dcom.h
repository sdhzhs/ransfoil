extern int iwd,iwu,iw0,iw,ip,jp,ic,jc,ib1,ib2,ifd,maxs;
extern double fd,fb,eb,lfar,delta,rau,rap,rae,rat,c,vfar,aoa,ta,tf,po,itur,tvr,ksi,cl,cd,cf,cm,xpc,ypc;
extern char cproctrl[8],cenergy[8],cvisheat[8],cturmod[8],cwalltreat[8],csolctrl[8],cdiscret[8],cinit[8],cstag[8],cpntctrl[8],clinsol[8],
ctmptype[8],cgtype[8],cmatair[8];
extern char cfilename[64],cdir[64],cmatfile[64];
extern double *cxwd,*cywd,*cxwu,*cywu,*cxwp0,*cywp0;
extern double *cxw,*cyw,*csw,*cyplus,*cystar,*chcv,*cax,*cay;
extern double *cxg,*cyg,*cxc,*cyc,*crho,*cmu,*cp,*cvx,*cvy,*ct,*ctn,*ctk,*cte,*ctw,*cmut;
void aero2d_(char *mode,int *prlv,char *scptname,size_t nmod,size_t nscpt);
void allocarray_(char *mode,size_t nmod);
void deallocarray_(char *mode,size_t nmod);