Subroutine Aero2D(libmod,prtl,scptname)
#include <config.h>
use Aero2DCOM
implicit none
integer s
integer(C_INT)::prtl
character(*) libmod,scptname
integer time0(8),time1(8)
real(8) rtime0,rtime1,cputime,walltime
logical(1) lconv,ldiv
#ifdef HYPRE
integer(8) matA,vecb,vecx,solver,precond
integer(1) solid
#endif
#ifdef HYPRE_GPU
!type(c_ptr) :: p_values
!integer(8) matA,vecb,vecx,solver,precond
integer(1) solid
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
 open(unit=2,file=filename(2),status='replace')
  write(2,*) 'Residual of RANS equations:'
  if(Turmod=='sa') then
   write(2,*) 'U equation   ','V equation   ','T equation   ','v equation   ','c equation   ','Iteration step'
  else if(Turmod=='ke') then
   write(2,*) 'U equation   ','V equation   ','T equation   ','k equation   ','e equation   ','c equation   ','Iteration step'
  else if(Turmod=='sst') then
   write(2,*) 'U equation   ','V equation   ','T equation   ','k equation   ','w equation   ','c equation   ','Iteration step'
  else if(Turmod=='lam'.or.Turmod=='inv') then
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
  Call Derivatives('mux')
  Call Derivatives('muy')
  Call Derivatives('mvx')
  Call Derivatives('mvy')
  Call Wallfunc
  Call Condiff('U')
#if defined(HYPRE)
  Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,U,U0,Rau,Ic,Jc,Ib1,Ib2,solid,'U')
#elif defined(HYPRE_GPU)
  Call hypresolve_gpu(aM,b,U,U0,Rau,Ic,Jc,Ib1,Ib2,solid,'U')
  !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,U,U0,Rau,Ic,Jc,Ib1,Ib2,solid,'U')
#else
  if(Linsol=='sor') then
   Call sor(aM,b,U,U0,Rau,Ic,Jc,Ib1,Ib2,'U')
  else if(Linsol=='pbicg') then
   Call CGSTAB(aM,b,U,U0,Rau,Ic,Jc,Ib1,Ib2,'U')
  end if
#endif
  Call Condiff('V')
#if defined(HYPRE)
  Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,V,V0,Rau,Ic,Jc,Ib1,Ib2,solid,'V')
#elif defined(HYPRE_GPU)
  Call hypresolve_gpu(aM,b,V,V0,Rau,Ic,Jc,Ib1,Ib2,solid,'V')
  !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,V,V0,Rau,Ic,Jc,Ib1,Ib2,solid,'V')
#else
  if(Linsol=='sor') then
   Call sor(aM,b,V,V0,Rau,Ic,Jc,Ib1,Ib2,'V')
  else if(Linsol=='pbicg') then
   Call CGSTAB(aM,b,V,V0,Rau,Ic,Jc,Ib1,Ib2,'V')
  end if
#endif
  Call PCEcoe
  rmsm=sum(abs(b))/(Ic*Jc)
#if defined(HYPRE)
  Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,dP,dP,1.0,Ic,Jc,Ib1,Ib2,solid,'dP')
#elif defined(HYPRE_GPU)
  Call hypresolve_gpu(aM,b,dP,dP,1.0,Ic,Jc,Ib1,Ib2,solid,'dP')
  !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,dP,dP,1.0,Ic,Jc,Ib1,Ib2,solid,'dP')
#else
  !Call sor(aM,b,dP,dP,1.0,Ic,Jc,Ib1,Ib2,'dP')
  Call CGSTAB(aM,b,dP,dP,1.0,Ic,Jc,Ib1,Ib2,'dP')
#endif
  if(Proctrl=='com') then
   Call StabSol('dP')
  end if
  P=P+Rap*dP
  Call Derivatives('dP')
  Call Velcorrect
  if(Proctrl=='com') then
   rho=rho+Rap*dP*Ma/(R*T)
   Call Densityinterp
  end if
  Call Derivatives('U')
  Call Derivatives('V')
  Call Derivatives('P')
  if(Energy=='Y') then
   Call Wallfunc
   Call Condiff('T')
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,T,T0,Rae,Ic,Jc,Ib1,Ib2,solid,'T')
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,T,T0,Rae,Ic,Jc,Ib1,Ib2,solid,'T')
   !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,T,T0,Rae,Ic,Jc,Ib1,Ib2,solid,'T')
#else
   if(Linsol=='sor') then
    Call sor(aM,b,T,T0,Rae,Ic,Jc,Ib1,Ib2,'T')
   else if(Linsol=='pbicg') then
    Call CGSTAB(aM,b,T,T0,Rae,Ic,Jc,Ib1,Ib2,'T')
   end if
#endif
   if(Proctrl=='com') then
    Call StabSol('T')
    rho=(Po+P)*Ma/(R*T)
    Call Densityinterp
   end if
   if(Turmod/='inv'.and.Matair=='Y') then
    mu=mu0*(T/Ti)**1.5*(Ti+Si)/(T+Si)
   end if
   Pr=mu*ca/ka
   Pc=9.24*((Pr/Prt)**0.75-1)*(1+0.28*exp(-0.007*Pr/Prt))
  end if
  if(Turmod=='sa') then
   Call Derivatives('Tn')
   Call Wallfunc
   Call Condiff('Tn')
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tn')
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tn')
   !Call hyprecompute_gpu(matA,vecb,vecx,solver,precond,p_values,aM,b,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tn') 
#else
   if(Linsol=='sor') then
    Call sor(aM,b,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,'Tn')
   else if(Linsol=='pbicg') then
    Call CGSTAB(aM,b,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,'Tn')
   end if
#endif
   Call StabSol('Tn')
  else if(Turmod=='ke') then
   Call Wallfunc
   Call Condiff('Tk')
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tk')
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tk')
#else
   if(Linsol=='sor') then
    Call sor(aM,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
   else if(Linsol=='pbicg') then
    Call CGSTAB(aM,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
   end if
#endif
   Call StabSol('Tk')
   Call Condiff('Te')
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,solid,'Te')
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,solid,'Te')
#else
   if(Linsol=='sor') then
    Call sor(aM,b,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,'Te')
   else if(Linsol=='pbicg') then
    Call CGSTAB(aM,b,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,'Te')
   end if
#endif
   Call StabSol('Te')
  else if(Turmod=='sst') then
   Call Derivatives('Tk')
   Call Derivatives('Tw')
   Call Wallfunc
   Call Condiff('Tk')
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tk')
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tk')
#else
   if(Linsol=='sor') then
    Call sor(aM,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
   else if(Linsol=='pbicg') then
    Call CGSTAB(aM,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
   end if
#endif
   Call StabSol('Tk')
   Call Condiff('Tw')
#if defined(HYPRE)
   Call hypresolve(matA,vecb,vecx,solver,precond,aM,b,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tw')
#elif defined(HYPRE_GPU)
   Call hypresolve_gpu(aM,b,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,solid,'Tw')
#else
   if(Linsol=='sor') then
    Call sor(aM,b,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,'Tw')
   else if(Linsol=='pbicg') then
    Call CGSTAB(aM,b,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,'Tw')
   end if
#endif
   Call StabSol('Tw')
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
  if(Turmod=='sa') then
   Call RMS(Tn,Tn0,Ic,Jc,rmsn)
  else if(Turmod=='ke') then
   Call RMS(Tk,Tk0,Ic,Jc,rmsk)
   Call RMS(Te,Te0,Ic,Jc,rmse)
  else if(Turmod=='sst') then
   Call RMS(Tk,Tk0,Ic,Jc,rmsk)
   Call RMS(Tw,Tw0,Ic,Jc,rmsw)
  end if
  if(prtl==1) then
   print *,'mass residual:'
   print *,rmsm,s
  end if
  if(libmod=='S'.or.libmod=='I'.or.libmod=='M') then
   if(Turmod=='sa') then
    write(2,'(5(ES11.4,2X),I5)') rmsu,rmsv,rmst,rmsn,rmsm,s
   else if(Turmod=='ke') then
    write(2,'(6(ES11.4,2X),I5)') rmsu,rmsv,rmst,rmsk,rmse,rmsm,s
   else if(Turmod=='sst') then
    write(2,'(6(ES11.4,2X),I5)') rmsu,rmsv,rmst,rmsk,rmsw,rmsm,s
   else if(Turmod=='lam'.or.Turmod=='inv') then
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
  if(Turmod=='sa') then
   Tn0=Tn
  else if(Turmod=='ke') then
   Tk0=Tk
   Te0=Te
  else if(Turmod=='sst') then
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
 Call Deallocarray(libmod)
 deallocate(Xwd,Ywd,Xwu,Ywu)
 if(Pntctrl=='Y') then
  deallocate(Xwp,Ywp,Xwp0,Ywp0)
 end if
else if(libmod=='C') then
 Call ArraypointerF2C
 if(Pntctrl=='Y') then
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

DO i=1,8
 if(cPntctrl(i)(1:1)/=C_NULL_CHAR) then
  Pntctrl(i:i)=cPntctrl(i)(1:1)
 else
  Pntctrl(i:i)=' '
 end if
 if(cProctrl(i)(1:1)/=C_NULL_CHAR) then
  Proctrl(i:i)=cProctrl(i)(1:1)
 else
  Proctrl(i:i)=' '
 end if
 if(cEnergy(i)(1:1)/=C_NULL_CHAR) then
  Energy(i:i)=cEnergy(i)(1:1)
 else
  Energy(i:i)=' '
 end if
 if(cvisheat(i)(1:1)/=C_NULL_CHAR) then
  visheat(i:i)=cvisheat(i)(1:1)
 else
  visheat(i:i)=' '
 end if
 if(cTurmod(i)(1:1)/=C_NULL_CHAR) then
  Turmod(i:i)=cTurmod(i)(1:1)
 else
  Turmod(i:i)=' '
 end if
 if(cWalltreat(i)(1:1)/=C_NULL_CHAR) then
  Walltreat(i:i)=cWalltreat(i)(1:1)
 else
  Walltreat(i:i)=' '
 end if
 if(csolctrl(i)(1:1)/=C_NULL_CHAR) then
  solctrl(i:i)=csolctrl(i)(1:1)
 else
  solctrl(i:i)=' '
 end if
 if(cDiscret(i)(1:1)/=C_NULL_CHAR) then
  Discret(i:i)=cDiscret(i)(1:1)
 else
  Discret(i:i)=' '
 end if
 if(cdenface(i)(1:1)/=C_NULL_CHAR) then
  denface(i:i)=cdenface(i)(1:1)
 else
  denface(i:i)=' '
 end if
 if(cLinsol(i)(1:1)/=C_NULL_CHAR) then
  Linsol(i:i)=cLinsol(i)(1:1)
 else
  Linsol(i:i)=' '
 end if
 if(cTmptype(i)(1:1)/=C_NULL_CHAR) then
  Tmptype(i:i)=cTmptype(i)(1:1)
 else
  Tmptype(i:i)=' '
 end if
 if(cInit(i)(1:1)/=C_NULL_CHAR) then
  Init(i:i)=cInit(i)(1:1)
 else
  Init(i:i)=' '
 end if
 if(cStag(i)(1:1)/=C_NULL_CHAR) then
  Stag(i:i)=cStag(i)(1:1)
 else
  Stag(i:i)=' '
 end if
 if(cgtype(i)(1:1)/=C_NULL_CHAR) then
  gtype(i:i)=cgtype(i)(1:1)
 else
  gtype(i:i)=' '
 end if
 if(cMatair(i)(1:1)/=C_NULL_CHAR) then
  Matair(i:i)=cMatair(i)(1:1)
 else
  Matair(i:i)=' '
 end if
end DO
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
if(Pntctrl=='Y') then
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
 if(Turmod=='sa'.or.Turmod=='sst') then
  Yplus=Yplus+Ypnw(:,i)
 else if(Turmod=='ke') then
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
if(Turmod=='sa'.or.Turmod=='sst') then
 Yplus=Yplus/(n+1)
else if(Turmod=='ke') then
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
divnum=1e+10

if(Turmod=='sa') then
 if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsn<delta) then
  lconv=.true.
 end if
 if(rmsm>divnum.or.rmsu>divnum.or.rmsv>divnum.or.rmst>divnum.or.rmsn>divnum.or.(.not.ieee_is_normal(rmsm)).or.(.not.ieee_is_normal(rmsu)).or.&
 (.not.ieee_is_normal(rmsv)).or.(.not.ieee_is_normal(rmst)).or.(.not.ieee_is_normal(rmsn))) then
  ldiv=.true.
 end if
else if(Turmod=='ke') then
 if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsk<delta.and.rmse<delta) then
  lconv=.true.
 end if
 if(rmsm>divnum.or.rmsu>divnum.or.rmsv>divnum.or.rmst>divnum.or.rmsk>divnum.or.rmse>divnum.or.(.not.ieee_is_normal(rmsm)).or.(.not.ieee_is_normal(rmsu)).or.&
 (.not.ieee_is_normal(rmsv)).or.(.not.ieee_is_normal(rmst)).or.(.not.ieee_is_normal(rmsk)).or.(.not.ieee_is_normal(rmse))) then
  ldiv=.true.
 end if
else if(Turmod=='sst') then
 if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsk<delta.and.rmsw<delta) then
  lconv=.true.
 end if
 if(rmsm>divnum.or.rmsu>divnum.or.rmsv>divnum.or.rmst>divnum.or.rmsk>divnum.or.rmsw>divnum.or.(.not.ieee_is_normal(rmsm)).or.(.not.ieee_is_normal(rmsu)).or.&
 (.not.ieee_is_normal(rmsv)).or.(.not.ieee_is_normal(rmst)).or.(.not.ieee_is_normal(rmsk)).or.(.not.ieee_is_normal(rmsw))) then
  ldiv=.true.
 end if
else if(Turmod=='lam'.or.Turmod=='inv') then
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
character(*) scalar

if(scalar=='dP') then
 dP=min(max(dP,(1-Po-P)/Rap),(5e+10-Po-P)/Rap)
else if(scalar=='T') then
 T=max(min(5e+3,T),1.0)
else if(scalar=='Tn') then
 Tn=min(max(Tn,0e+0),1e+5*mu/rho)
else if(scalar=='Tk') then
 Tk=max(Tk,1e-15)
else if(scalar=='Te') then
 Te=max(Te,(rho*Cu*Tk**2)/(1e+5*mu))
else if(scalar=='Tw') then
 Tw=max(Tw,(rho*Tk)/(1e+5*mu))
end if

end Subroutine StabSol

Subroutine Allocarray(libmod)
use Aero2DCOM
implicit none
integer i
character(*) libmod

if(libmod=='C') then
 DO i=1,8
  if(cTurmod(i)(1:1)/=C_NULL_CHAR) then
   Turmod(i:i)=cTurmod(i)(1:1)
  else
   Turmod(i:i)=' '
  end if
 end DO
end if
allocate(U0(Ic,Jc),V0(Ic,Jc),T0(Ic,Jc),U(Ic,Jc),V(Ic,Jc),T(Ic,Jc),rho(Ic,Jc),mu(Ic,Jc),P(Ic,Jc),dP(Ic,Jc),mut(Ic,Jc),&
Pr(Ic,Jc),Pc(Ic,Jc),auP(Ic,Jc),auNB(Ic,Jc),aM(5,Ic,Jc),b(Ic,Jc),Xg(Ip,Jp),Yg(Ip,Jp),Xc(Ic,Jc),Yc(Ic,Jc),Xga(Ic,Jc),Xgk(Ic,Jc),&
Yga(Ic,Jc),Ygk(Ic,Jc),dk(Ic,Jc),da(Ic,Jc),Jg(Ic,Jc),a1(Ic,Jc),y1(Ic,Jc),b1(Ic,Jc),Un(Ic,Jc),Vn(Ic,Jc),Unk(Ip,Jc),Vna(Ic,Jp),&
duk(Ip,Jc),dva(Ic,Jp),Ux(Ic,Jc),Uy(Ic,Jc),Vx(Ic,Jc),Vy(Ic,Jc),Px(Ic,Jc),Py(Ic,Jc),dPx(Ic,Jc),dPy(Ic,Jc),muxx(Ic,Jc),muxy(Ic,Jc),&
muyx(Ic,Jc),mvxy(Ic,Jc),mvyx(Ic,Jc),mvyy(Ic,Jc),rhok(Ip,Jc),rhoa(Ic,Jp),d(Ic,Jc))
if(Turmod=='sa') then
 allocate(Tn0(Ic,Jc),Tn(Ic,Jc),Tnx(Ic,Jc),Tny(Ic,Jc))
else if(Turmod=='ke') then
 allocate(Tk0(Ic,Jc),Te0(Ic,Jc),Tk(Ic,Jc),Te(Ic,Jc),rhox(Ic,Jc),rhoy(Ic,Jc))
else if(Turmod=='sst') then
 allocate(Tk0(Ic,Jc),Tw0(Ic,Jc),Tk(Ic,Jc),Tw(Ic,Jc),Tkx(Ic,Jc),Tky(Ic,Jc),Twx(Ic,Jc),Twy(Ic,Jc),sigmatk(Ic,Jc),sigmatw(Ic,Jc))
end if
allocate(Xw(Ib1:Ib2),Yw(Ib1:Ib2),Yp(Ib1:Ib2),DR(Ib1:Ib2),Sw(Ib1:Ib2),ks(Ib1:Ib2),Q(Ib1:Ib2),Yplus(Ib1:Ib2),Ystar(Ib1:Ib2),&
ustar(Ib1:Ib2),Uplus(Ib1:Ib2),Tplus(Ib1:Ib2),hcv(Ib1:Ib2),Ax(Ib1:Ib2),Ay(Ib1:Ib2))
allocate(Pnw(Ib1:Ib2,3),Unw(Ib1:Ib2,3),Vnw(Ib1:Ib2,3),Tnw(Ib1:Ib2,3),mutnw(Ib1:Ib2,3),hcvnw(Ib1:Ib2,3),Axnw(Ib1:Ib2,3),Aynw(Ib1:Ib2,3),Ypnw(Ib1:Ib2,3))
allocate(Clrec(3),Cdrec(3),Cfrec(3),Cmrec(3),Xpcrec(3),Ypcrec(3))

end Subroutine Allocarray

Subroutine Deallocarray(libmod)
use Aero2DCOM
implicit none
integer i
character(*) libmod

if(libmod=='C') then
 DO i=1,8
 if(cTurmod(i)(1:1)/=C_NULL_CHAR) then
  Turmod(i:i)=cTurmod(i)(1:1)
 else
  Turmod(i:i)=' '
 end if
 end DO
end if
deallocate(U0,V0,T0,U,V,T,rho,mu,mut,P,dP,Pr,Pc,auP,auNB,aM,b,Xg,Yg,Xc,Yc,Xga,Xgk,Yga,Ygk,dk,da,Jg,a1,y1,b1,&
Un,Vn,Unk,Vna,duk,dva,Ux,Uy,Vx,Vy,Px,Py,dPx,dPy,muxx,muxy,muyx,mvxy,mvyx,mvyy,rhok,rhoa,d)
if(Turmod=='sa') then
 deallocate(Tn0,Tn,Tnx,Tny)
else if(Turmod=='ke') then
 deallocate(Tk0,Te0,Tk,Te,rhox,rhoy)
else if(Turmod=='sst') then
 deallocate(Tk0,Tw0,Tk,Tw,Tkx,Tky,Twx,Twy,sigmatk,sigmatw)
end if
deallocate(Xw,Yw,Yp,DR,Sw,ks,Q,Yplus,Ystar,ustar,Uplus,Tplus,hcv,Ax,Ay)

end Subroutine Deallocarray
