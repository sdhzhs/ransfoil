extern int iwd,iwu,ip,jp,ic,jc,ib1,ib2,maxs;
extern double fd,delta,rau,rap,rae,rar,rat,c,vfar,aoa,ta,tf,po,itur,tvr,ksi,cl,cd,cf,cm,xpc,ypc;
extern char cproctrl[8],cenergy[8],cvisheat[8],cturmod[8],cwalltreat[8],csolctrl[8],cdiscret[8],cinit[8],cstag[8];
extern char cfilename[64],cdir[64];
extern double *cxwd,*cywd,*cxwu,*cywu;
extern double *cxw,*cyw,*csw,*cyplus,*cystar,*chcv,*cax,*cay;
extern double *cxg,*cyg,*cxc,*cyc,*crou,*cmiu,*cp,*cvx,*cvy,*ct,*ctn,*ctk,*cte,*ctw,*cmiut;
void aero2d_(char *mode,int *prlv,char *scptname,size_t nmod,size_t nscpt);
void allocarray_(void);
void deallocarray_(void);
