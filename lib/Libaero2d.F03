Subroutine Aero2D(libmod,prtl,scptname)
#include <config.h>
use Aero2DCOM
implicit none
integer(C_INT)::prtl
character(*) libmod,scptname
integer s,ioerr
character(128) ioerrmsg
integer time0(8),time1(8)
real(8) rtime0,rtime1,cputime,walltime
logical(1) lconv,ldiv
#ifdef HYPRE
integer(8) matA,vecb,vecx,solver,precond
integer solid
#endif
#ifdef HYPRE_GPU
!type(c_ptr) :: p_values
!integer(8) matA,vecb,vecx,solver,precond
integer solid
#endif

if(libmod=='S'.or.libmod=='I'.or.libmod=='M') then
 Call Readpara(libmod,scptname)
else if(libmod=='C') then
 Call ConfigparaC2F
end if
Call Genmesh(libmod)
Call Jacobi
Call Initial
Call Densityinterp
if(libmod=='S'.or.libmod=='I'.or.libmod=='M') then
 filename(2)=trim(dir)//'/iteration.txt'
 open(unit=2,file=filename(2),status='replace',IOSTAT=ioerr,IOMSG=ioerrmsg,recl=128)
  if(ioerr>0) stop ioerrmsg
  write(2,*) 'Residual of RANS equations:'
  if(TurmodFlag==SA) then
   write(2,*) 'U equation   ','V equation   ','T equation   ','v equation   ','c equation   ','Iteration step'
  else if(TurmodFlag==KE) then
   write(2,*) 'U equation   ','V equation   ','T equation   ','k equation   ','e equation   ','c equation   ','Iteration step'
  else if(TurmodFlag==SST) then
   write(2,*) 'U equation   ','V equation   ','T equation   ','k equation   ','w equation   ','c equation   ','Iteration step'
  else if(TurmodFlag==LAM.or.TurmodFlag==INV) then
   write(2,*) 'U equation   ','V equation   ','T equation   ','c equation   ','Iteration step'
  end if
end if
#ifdef HYPRE
  solid=3
  Call hypreinit(matA,vecb,vecx,solver,precond,Ic,Jc,Ib1,Ib2,solid)
#endif
#ifdef HYPRE_GPU
  solid=3
  !Call hypreinit_gpu(matA,vecb,vecx,solver,precond,p_values,Ic,Jc,Ib1,Ib2,solid)
#endif
print *,'RANS equations iteration is beginning...'
Call DATE_AND_TIME(VALUES=time0)
Call CPU_TIME(rtime0)
DO s=1,maxs
  Call Derivatives(MUX)
  Call Derivatives(MUY)
  Call Derivatives(MVX)
  Call Derivatives(MVY)
  Call Wallfunc
  Call Condiff(VELX)
  Call ImplicitRelax(aM,b,U0,Rau,Ic,Jc,Ib1,Ib2,VELX)
#if defined(HYPRE)
  Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,U,Ic,Jc,Ib1,solid,VELX,FstypeFlag)
#elif defined(HYPRE_GPU)
  Call hypresolve_gpu(aM,b,U,Ic,Jc,Ib1,Ib2,solid,VELX,FstypeFlag)
  !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,U,Ic,Jc,Ib1,solid,VELX,FstypeFlag)
#else
  if(LinsolFlag==SORGS) then
   Call sor(aM,b,U,Ic,Jc,Ib1,Ib2,VELX,FstypeFlag)
  else if(LinsolFlag==PBICG) then
   Call CGSTAB(aM,b,U,Ic,Jc,Ib1,Ib2,VELX,FstypeFlag)
  end if
#endif
  Call Condiff(VELY)
  Call ImplicitRelax(aM,b,V0,Rau,Ic,Jc,Ib1,Ib2,VELY)
#if defined(HYPRE)
  Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,V,Ic,Jc,Ib1,solid,VELY,FstypeFlag)
#elif defined(HYPRE_GPU)
  Call hypresolve_gpu(aM,b,V,Ic,Jc,Ib1,Ib2,solid,VELY,FstypeFlag)
  !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,V,Ic,Jc,Ib1,solid,VELY,FstypeFlag)
#else
  if(LinsolFlag==SORGS) then
   Call sor(aM,b,V,Ic,Jc,Ib1,Ib2,VELY,FstypeFlag)
  else if(LinsolFlag==PBICG) then
   Call CGSTAB(aM,b,V,Ic,Jc,Ib1,Ib2,VELY,FstypeFlag)
  end if
#endif
  Call PCEcoe
#if defined(HYPRE)
  Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,dP,Ic,Jc,Ib1,solid,DPRES,FstypeFlag)
#elif defined(HYPRE_GPU)
  Call hypresolve_gpu(aM,b,dP,Ic,Jc,Ib1,Ib2,solid,DPRES,FstypeFlag)
  !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,dP,Ic,Jc,Ib1,solid,DPRES,FstypeFlag)
#else
  !Call sor(aM,b,dP,Ic,Jc,Ib1,Ib2,DPRES,FstypeFlag)
  Call CGSTAB(aM,b,dP,Ic,Jc,Ib1,Ib2,DPRES,FstypeFlag)
#endif
  if(ProctrlFlag==COM) then
   Call StabSol(DPRES)
  end if
  P=P+Rap*dP
  Call Derivatives(DPRES)
  Call Velcorrect
  if(ProctrlFlag==COM) then
   if(FstypeFlag==VINPOUT) then
    rho(:,1:Jc-1)=rho(:,1:Jc-1)+Rap*dP(:,1:Jc-1)*Ma/(R*T(:,1:Jc-1))
    where(Vna(:,Jc)>0.0) rho(:,Jc)=rho(:,Jc)+Rap*dP(:,Jc)*Ma/(R*T(:,Jc))
   else 
    rho=rho+Rap*dP*Ma/(R*T)
   end if
   Call Derivatives(DENS)
   Call Densityinterp
  end if
  Call Derivatives(VELX)
  Call Derivatives(VELY)
  Call Derivatives(PRES)
  if(EnergyFlag) then
   Call Wallfunc
   Call Condiff(TEMP)
   Call ImplicitRelax(aM,b,T0,Rae,Ic,Jc,Ib1,Ib2,TEMP)
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,T,Ic,Jc,Ib1,solid,TEMP,FstypeFlag)
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,T,Ic,Jc,Ib1,Ib2,solid,TEMP,FstypeFlag)
   !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,T,Ic,Jc,Ib1,solid,TEMP,FstypeFlag)
#else
   if(LinsolFlag==SORGS) then
    Call sor(aM,b,T,Ic,Jc,Ib1,Ib2,TEMP,FstypeFlag)
   else if(LinsolFlag==PBICG) then
    Call CGSTAB(aM,b,T,Ic,Jc,Ib1,Ib2,TEMP,FstypeFlag)
   end if
#endif
   if(ProctrlFlag==COM) then
    Call StabSol(TEMP)
    if(FstypeFlag==VINPOUT) then
     rho(:,1:Jc-1)=(Po+P(:,1:Jc-1))*Ma/(R*T(:,1:Jc-1))
     where(Vna(:,Jc)>0.0) rho(:,Jc)=(Po+P(:,Jc))*Ma/(R*T(:,Jc))
    else
     rho=(Po+P)*Ma/(R*T)
    end if
    Call Derivatives(DENS)
    Call Densityinterp
   end if
   Call Derivatives(TEMP)
   if(TurmodFlag/=INV.and.MatairFlag) then
    mu=mu0*(T/Ti)**1.5*(Ti+Si)/(T+Si)
   end if
   Pr=mu*ca/ka
   Pc=9.24*((Pr/Prt)**0.75-1)*(1+0.28*exp(-0.007*Pr/Prt))
  end if
  if(TurmodFlag==SA) then
   Call Wallfunc
   Call Condiff(TURBNU)
   Call ImplicitRelax(aM,b,Tn0,Rat,Ic,Jc,Ib1,Ib2,TURBNU)
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Tn,Ic,Jc,Ib1,solid,TURBNU,FstypeFlag)
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Tn,Ic,Jc,Ib1,Ib2,solid,TURBNU,FstypeFlag)
   !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,Tn,Ic,Jc,Ib1,solid,TURBNU,FstypeFlag) 
#else
   if(LinsolFlag==SORGS) then
    Call sor(aM,b,Tn,Ic,Jc,Ib1,Ib2,TURBNU,FstypeFlag)
   else if(LinsolFlag==PBICG) then
    Call CGSTAB(aM,b,Tn,Ic,Jc,Ib1,Ib2,TURBNU,FstypeFlag)
   end if
#endif
   Call StabSol(TURBNU)
   Call Derivatives(TURBNU)
  else if(TurmodFlag==KE) then
   Call Wallfunc
   Call Condiff(TURBK)
   Call ImplicitRelax(aM,b,Tk0,Rat,Ic,Jc,Ib1,Ib2,TURBK)
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Tk,Ic,Jc,Ib1,solid,TURBK,FstypeFlag)
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Tk,Ic,Jc,Ib1,Ib2,solid,TURBK,FstypeFlag)
#else
   if(LinsolFlag==SORGS) then
    Call sor(aM,b,Tk,Ic,Jc,Ib1,Ib2,TURBK,FstypeFlag)
   else if(LinsolFlag==PBICG) then
    Call CGSTAB(aM,b,Tk,Ic,Jc,Ib1,Ib2,TURBK,FstypeFlag)
   end if
#endif
   Call StabSol(TURBK)
   Call Condiff(TURBE)
   Call ImplicitRelax(aM,b,Te0,Rat,Ic,Jc,Ib1,Ib2,TURBE)
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Te,Ic,Jc,Ib1,solid,TURBE,FstypeFlag)
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Te,Ic,Jc,Ib1,Ib2,solid,TURBE,FstypeFlag)
#else
   if(LinsolFlag==SORGS) then
    Call sor(aM,b,Te,Ic,Jc,Ib1,Ib2,TURBE,FstypeFlag)
   else if(LinsolFlag==PBICG) then
    Call CGSTAB(aM,b,Te,Ic,Jc,Ib1,Ib2,TURBE,FstypeFlag)
   end if
#endif
   Call StabSol(TURBE)
   Call Derivatives(TURBK)
   Call Derivatives(TURBE)
  else if(TurmodFlag==SST) then
   Call Wallfunc
   Call Condiff(TURBK)
   Call ImplicitRelax(aM,b,Tk0,Rat,Ic,Jc,Ib1,Ib2,TURBK)
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Tk,Ic,Jc,Ib1,solid,TURBK,FstypeFlag)
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Tk,Ic,Jc,Ib1,Ib2,solid,TURBK,FstypeFlag)
#else
   if(LinsolFlag==SORGS) then
    Call sor(aM,b,Tk,Ic,Jc,Ib1,Ib2,TURBK,FstypeFlag)
   else if(LinsolFlag==PBICG) then
    Call CGSTAB(aM,b,Tk,Ic,Jc,Ib1,Ib2,TURBK,FstypeFlag)
   end if
#endif
   Call StabSol(TURBK)
   Call Condiff(TURBW)
   Call ImplicitRelax(aM,b,Tw0,Rat,Ic,Jc,Ib1,Ib2,TURBW)
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Tw,Ic,Jc,Ib1,solid,TURBW,FstypeFlag)
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Tw,Ic,Jc,Ib1,Ib2,solid,TURBW,FstypeFlag)
#else
   if(LinsolFlag==SORGS) then
    Call sor(aM,b,Tw,Ic,Jc,Ib1,Ib2,TURBW,FstypeFlag)
   else if(LinsolFlag==PBICG) then
    Call CGSTAB(aM,b,Tw,Ic,Jc,Ib1,Ib2,TURBW,FstypeFlag)
   end if
#endif
   Call StabSol(TURBW)
   Call Derivatives(TURBK)
   Call Derivatives(TURBW)
  end if
  Call Turvis
  Call Postprocess(s)
  if((libmod=='S'.or.libmod=='I'.or.libmod=='M').and.(s==1.or.mod(s,100)==0)) then
   Call Saveresults
   print *,'Output results completed!'
  end if
  Call RMS(U,U0,Ic,Jc,rmsu)
  Call RMS(V,V0,Ic,Jc,rmsv)
  Call RMS(T,T0,Ic,Jc,rmst)
  if(TurmodFlag==SA) then
   Call RMS(Tn,Tn0,Ic,Jc,rmsn)
  else if(TurmodFlag==KE) then
   Call RMS(Tk,Tk0,Ic,Jc,rmsk)
   Call RMS(Te,Te0,Ic,Jc,rmse)
  else if(TurmodFlag==SST) then
   Call RMS(Tk,Tk0,Ic,Jc,rmsk)
   Call RMS(Tw,Tw0,Ic,Jc,rmsw)
  end if
  if(prtl==1) then
   print *,'mass residual:'
   print *,rmsm,s
  end if
  if(libmod=='S'.or.libmod=='I'.or.libmod=='M') then
   if(TurmodFlag==SA) then
    write(2,'(5(ES11.4,2X),I5)') rmsu,rmsv,rmst,rmsn,rmsm,s
   else if(TurmodFlag==KE) then
    write(2,'(6(ES11.4,2X),I5)') rmsu,rmsv,rmst,rmsk,rmse,rmsm,s
   else if(TurmodFlag==SST) then
    write(2,'(6(ES11.4,2X),I5)') rmsu,rmsv,rmst,rmsk,rmsw,rmsm,s
   else if(TurmodFlag==LAM.or.TurmodFlag==INV) then
    write(2,'(4(ES11.4,2X),I5)') rmsu,rmsv,rmst,rmsm,s
   end if
  end if
  Call Convergence(lconv,ldiv)
  if(lconv) then
   print *,'RANS Equations converge!'
   exit
  end if
  if(ldiv) then
   stop 'RANS Equations diverge!'
  end if
  U0=U
  V0=V
  T0=T
  if(TurmodFlag==SA) then
   Tn0=Tn
  else if(TurmodFlag==KE) then
   Tk0=Tk
   Te0=Te
  else if(TurmodFlag==SST) then
   Tk0=Tk
   Tw0=Tw
  end if
end DO
Call DATE_AND_TIME(VALUES=time1)
Call CPU_TIME(rtime1)
cputime=rtime1-rtime0
if(time1(6)-time0(6)>=0) then
 if(time1(7)-time0(7)>=0) then
  print *,'Iteration elapsed time:',time1(6)-time0(6),'m',time1(7)-time0(7),'s'
 else
  print *,'Iteration elapsed time:',time1(6)-time0(6)-1,'m',60+time1(7)-time0(7),'s'
 end if
 walltime=60.*(time1(6)-time0(6))+time1(7)-time0(7)
else
 if(time1(7)-time0(7)>=0) then
  print *,'Iteration elapsed time:',60+time1(6)-time0(6),'m',time1(7)-time0(7),'s'
 else
  print *,'Iteration elapsed time:',60+time1(6)-time0(6)-1,'m',60+time1(7)-time0(7),'s'
 end if
 walltime=60.*(60+time1(6)-time0(6))+time1(7)-time0(7)
end if
write(*,'(1X,A19,14X,F7.2,A1)') 'Iteration CPU time:',cputime,'s'
write(*,'(1X,A8,25X,F7.2)') 'Scaling:',cputime/walltime
if(s==maxs+1) then
 print *,'Reach maximum iteration steps!'
 if(s>1500) then
  Call AverageSol
 end if
end if
if(libmod=='S'.or.libmod=='I'.or.libmod=='M') then
 close(2)
 Call Saveresults
 print *,'Save results successfully!'
end if
if(libmod=='S'.or.libmod=='I') then
 Call Deallocarray
 deallocate(Xwd,Ywd,Xwu,Ywu)
 if(PntctrlFlag) then
  deallocate(Xwp,Ywp,Xwp0,Ywp0)
 end if
else if(libmod=='C') then
 Call ArraypointerF2C
 if(PntctrlFlag) then
  deallocate(Xwp,Ywp,Xwp0,Ywp0)
 else
  deallocate(Xwd,Ywd,Xwu,Ywu)
 end if
end if
#ifdef HYPRE
  Call hyprerelease(matA,vecb,vecx,solver,precond,solid)
#endif
#ifdef HYPRE_GPU
  !Call hyprerelease_gpu(matA,vecb,vecx,solver,precond,p_values,solid)
#endif

end Subroutine Aero2D

Subroutine ConfigparaC2F
use Aero2DCOM
implicit none
integer i

DO i=1,64
 if(cfilename(i)(1:1)/=C_NULL_CHAR) then
  filename(9)(i:i)=cfilename(i)(1:1)
 else
  filename(9)(i:i)=' '
 end if
 if(cdir(i)(1:1)/=C_NULL_CHAR) then
  dir(i:i)=cdir(i)(1:1)
 else
  dir(i:i)=' '
 end if
 if(cmatfile(i)(1:1)/=C_NULL_CHAR) then
  matfile(i:i)=cmatfile(i)(1:1)
 else
  matfile(i:i)=' '
 end if
end DO

end Subroutine ConfigparaC2F

Subroutine ArraypointerF2C
use Aero2DCOM
implicit none

cXw=C_LOC(Xw)
cYw=C_LOC(Yw)
cSw=C_LOC(Sw)
cYplus=C_LOC(Yplus)
cYstar=C_LOC(Ystar)
chcv=C_LOC(hcv)
cAx=C_LOC(Ax)
cAy=C_LOC(Ay)
cXg=C_LOC(Xg)
cYg=C_LOC(Yg)
cXc=C_LOC(Xc)
cYc=C_LOC(Yc)
crho=C_LOC(rho)
cmu=C_LOC(mu)
cP=C_LOC(P)
cVx=C_LOC(U)
cVy=C_LOC(V)
cT=C_LOC(T)
cTn=C_LOC(Tn)
cTk=C_LOC(Tk)
cTe=C_LOC(Te)
cTw=C_LOC(Tw)
cmut=C_LOC(mut)
if(PntctrlFlag) then
 cXwd=C_LOC(Xwd)
 cYwd=C_LOC(Ywd)
 cXwu=C_LOC(Xwu)
 cYwu=C_LOC(Ywu)
end if

end Subroutine ArraypointerF2C

Subroutine AverageSol
use Aero2DCOM
implicit none
integer i,n

n=3
DO i=1,n
 Cl=Cl+Clrec(i)
 Cd=Cd+Cdrec(i)
 Cf=Cf+Cfrec(i)
 Cm=Cm+Cmrec(i)
 Xpc=Xpc+Xpcrec(i)
 Ypc=Ypc+Ypcrec(i)
 hcv=hcv+hcvnw(:,i)
 Ax=Ax+Axnw(:,i)
 Ay=Ay+Aynw(:,i)
 if(TurmodFlag==SA.or.TurmodFlag==SST) then
  Yplus=Yplus+Ypnw(:,i)
 else if(TurmodFlag==KE) then
  Ystar=Ystar+Ypnw(:,i)
 end if
 P(Ib1:Ib2,1)=P(Ib1:Ib2,1)+Pnw(:,i)
 U(Ib1:Ib2,1)=U(Ib1:Ib2,1)+Unw(:,i)
 V(Ib1:Ib2,1)=V(Ib1:Ib2,1)+Vnw(:,i)
 T(Ib1:Ib2,1)=T(Ib1:Ib2,1)+Tnw(:,i)
 mut(Ib1:Ib2,1)=mut(Ib1:Ib2,1)+mutnw(:,i)
end DO
Cl=Cl/(n+1)
Cd=Cd/(n+1)
Cf=Cf/(n+1)
Cm=Cm/(n+1)
Xpc=Xpc/(n+1)
Ypc=Ypc/(n+1)
hcv=hcv/(n+1)
Ax=Ax/(n+1)
Ay=Ay/(n+1)
if(TurmodFlag==SA.or.TurmodFlag==SST) then
 Yplus=Yplus/(n+1)
else if(TurmodFlag==KE) then
 Ystar=Ystar/(n+1)
end if
P(Ib1:Ib2,1)=P(Ib1:Ib2,1)/(n+1)
U(Ib1:Ib2,1)=U(Ib1:Ib2,1)/(n+1)
V(Ib1:Ib2,1)=V(Ib1:Ib2,1)/(n+1)
T(Ib1:Ib2,1)=T(Ib1:Ib2,1)/(n+1)
mut(Ib1:Ib2,1)=mut(Ib1:Ib2,1)/(n+1)

end Subroutine AverageSol

Subroutine Convergence(lconv,ldiv)
use ieee_arithmetic
use Aero2DCOM
implicit none
logical(1) lconv,ldiv
real(8) divnum

lconv=.false.
ldiv=.false.
divnum=1e+20

if(TurmodFlag==SA) then
 if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsn<delta) then
  lconv=.true.
 end if
 if(rmsm>divnum.or.rmsu>divnum.or.rmsv>divnum.or.rmst>divnum.or.rmsn>divnum.or.(.not.ieee_is_normal(rmsm)).or.(.not.ieee_is_normal(rmsu)).or.&
 (.not.ieee_is_normal(rmsv)).or.(.not.ieee_is_normal(rmst)).or.(.not.ieee_is_normal(rmsn))) then
  ldiv=.true.
 end if
else if(TurmodFlag==KE) then
 if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsk<delta.and.rmse<delta) then
  lconv=.true.
 end if
 if(rmsm>divnum.or.rmsu>divnum.or.rmsv>divnum.or.rmst>divnum.or.rmsk>divnum.or.rmse>divnum.or.(.not.ieee_is_normal(rmsm)).or.(.not.ieee_is_normal(rmsu)).or.&
 (.not.ieee_is_normal(rmsv)).or.(.not.ieee_is_normal(rmst)).or.(.not.ieee_is_normal(rmsk)).or.(.not.ieee_is_normal(rmse))) then
  ldiv=.true.
 end if
else if(TurmodFlag==SST) then
 if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsk<delta.and.rmsw<delta) then
  lconv=.true.
 end if
 if(rmsm>divnum.or.rmsu>divnum.or.rmsv>divnum.or.rmst>divnum.or.rmsk>divnum.or.rmsw>divnum.or.(.not.ieee_is_normal(rmsm)).or.(.not.ieee_is_normal(rmsu)).or.&
 (.not.ieee_is_normal(rmsv)).or.(.not.ieee_is_normal(rmst)).or.(.not.ieee_is_normal(rmsk)).or.(.not.ieee_is_normal(rmsw))) then
  ldiv=.true.
 end if
else if(TurmodFlag==LAM.or.TurmodFlag==INV) then
 if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta) then
  lconv=.true.
 end if
 if(rmsm>divnum.or.rmsu>divnum.or.rmsv>divnum.or.rmst>divnum.or.(.not.ieee_is_normal(rmsm)).or.(.not.ieee_is_normal(rmsu)).or.(.not.ieee_is_normal(rmsv)).or.&
 (.not.ieee_is_normal(rmst))) then
  ldiv=.true.
 end if
end if

end Subroutine Convergence

Subroutine StabSol(scalar)
use Aero2DCOM
implicit none
integer scalar

if(scalar==DPRES) then
 dP=min(max(dP,(1-Po-P)/Rap),(5e+10-Po-P)/Rap)
else if(scalar==TEMP) then
 T=max(min(5e+3,T),1.0)
else if(scalar==TURBNU) then
 Tn=min(max(Tn,1e-15),1e+5*mu/rho)
else if(scalar==TURBK) then
 Tk=max(Tk,1e-15)
else if(scalar==TURBE) then
 Te=max(Te,(rho*Cu*Tk**2)/(1e+5*mu))
else if(scalar==TURBW) then
 Tw=max(Tw,(rho*Tk)/(1e+5*mu))
end if

end Subroutine StabSol

Subroutine Allocarray
use Aero2DCOM
implicit none

allocate(U0(Ic,Jc),V0(Ic,Jc),T0(Ic,Jc),U(Ic,Jc),V(Ic,Jc),T(Ic,Jc),rho(Ic,Jc),mu(Ic,Jc),P(Ic,Jc),dP(Ic,Jc),mut(Ic,Jc),&
Pr(Ic,Jc),Pc(Ic,Jc),auP(Ic,Jc),auNB(Ic,Jc),aM(5,Ic,Jc),b(Ic,Jc),Xg(Ip,Jp),Yg(Ip,Jp),Xc(Ic,Jc),Yc(Ic,Jc),Xfa(Ic,Jp),Yfa(Ic,Jp),&
Xfk(Ip,Jc),Yfk(Ip,Jc),dkw(Ip,Jc),daw(Ic,Jp),dkd(Ip,Jc),dad(Ic,Jp),dk(Ic,Jc),da(Ic,Jc),Vol(Ic,Jc),Unk(Ip,Jc),Vna(Ic,Jp),duk(Ip,Jc),dva(Ic,Jp),&
Ux(Ic,Jc),Uy(Ic,Jc),Vx(Ic,Jc),Vy(Ic,Jc),Px(Ic,Jc),Py(Ic,Jc),dPx(Ic,Jc),dPy(Ic,Jc),Tx(Ic,Jc),Ty(Ic,Jc),rhox(Ic,Jc),rhoy(Ic,Jc),muxx(Ic,Jc),&
muxy(Ic,Jc),muyx(Ic,Jc),mvxy(Ic,Jc),mvyx(Ic,Jc),mvyy(Ic,Jc),rhok(Ip,Jc),rhoa(Ic,Jp),d(Ic,Jc))
if(TurmodFlag==SA) then
 allocate(Tn0(Ic,Jc),Tn(Ic,Jc),Tnx(Ic,Jc),Tny(Ic,Jc))
else if(TurmodFlag==KE) then
 allocate(Tk0(Ic,Jc),Te0(Ic,Jc),Tk(Ic,Jc),Te(Ic,Jc),Tkx(Ic,Jc),Tky(Ic,Jc),Tex(Ic,Jc),Tey(Ic,Jc))
else if(TurmodFlag==SST) then
 allocate(Tk0(Ic,Jc),Tw0(Ic,Jc),Tk(Ic,Jc),Tw(Ic,Jc),Tkx(Ic,Jc),Tky(Ic,Jc),Twx(Ic,Jc),Twy(Ic,Jc),sigmatk(Ic,Jc),sigmatw(Ic,Jc))
end if
allocate(Xw(Ib1:Ib2),Yw(Ib1:Ib2),Yp(Ib1:Ib2),DR(Ib1:Ib2),Sw(Ib1:Ib2),ks(Ib1:Ib2),Q(Ib1:Ib2),Yplus(Ib1:Ib2),Ystar(Ib1:Ib2),&
ustar(Ib1:Ib2),Uplus(Ib1:Ib2),Tplus(Ib1:Ib2),hcv(Ib1:Ib2),Ax(Ib1:Ib2),Ay(Ib1:Ib2))
allocate(Pnw(Ib1:Ib2,3),Unw(Ib1:Ib2,3),Vnw(Ib1:Ib2,3),Tnw(Ib1:Ib2,3),mutnw(Ib1:Ib2,3),hcvnw(Ib1:Ib2,3),Axnw(Ib1:Ib2,3),Aynw(Ib1:Ib2,3),Ypnw(Ib1:Ib2,3))
allocate(Clrec(3),Cdrec(3),Cfrec(3),Cmrec(3),Xpcrec(3),Ypcrec(3))

end Subroutine Allocarray

Subroutine Deallocarray
use Aero2DCOM
implicit none

deallocate(U0,V0,T0,U,V,T,rho,mu,mut,P,dP,Pr,Pc,auP,auNB,aM,b,Xg,Yg,Xc,Yc,Xfa,Yfa,Xfk,Yfk,dkw,daw,Vol,dkd,dad,dk,da,&
Unk,Vna,duk,dva,Ux,Uy,Vx,Vy,Px,Py,dPx,dPy,Tx,Ty,rhox,rhoy,muxx,muxy,muyx,mvxy,mvyx,mvyy,rhok,rhoa,d)
if(TurmodFlag==SA) then
 deallocate(Tn0,Tn,Tnx,Tny)
else if(TurmodFlag==KE) then
 deallocate(Tk0,Te0,Tk,Te,Tkx,Tky,Tex,Tey)
else if(TurmodFlag==SST) then
 deallocate(Tk0,Tw0,Tk,Tw,Tkx,Tky,Twx,Twy,sigmatk,sigmatw)
end if
deallocate(Xw,Yw,Yp,DR,Sw,ks,Q,Yplus,Ystar,ustar,Uplus,Tplus,hcv,Ax,Ay)

end Subroutine Deallocarray