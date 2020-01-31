Subroutine Aero2D(libmod,prtl,scptname)
#include <config.h>
use Aero2DCOM
implicit none
integer i,s
integer(C_INT),bind(C)::prtl
character(*) libmod,scptname
integer time0(8),time1(8)

if(libmod=='S'.or.libmod=='I') then
 Call Readpara(libmod,scptname)
else if(libmod=='C') then
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
 end DO
end if
Call Genmesh(libmod)
Call Jacobi
Call Initial
Call Densityinterp
if(libmod=='S'.or.libmod=='I') then
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
print *,'RANS equations iteration is beginning...'
Call DATE_AND_TIME(VALUES=time0)
DO s=1,maxs
  Call Derivatives('mux')
  Call Derivatives('muy')
  Call Derivatives('mvx')
  Call Derivatives('mvy')
  Call Wallfunc
  Call Condiff('U')
#ifdef HYPRE
  Call hypresolve(aP,aW,aE,aS,aN,b,U,U0,Rau,Ic,Jc,Ib1,Ib2,3,'U')
#else
  Call sor(aP,aW,aE,aS,aN,b,U,U0,Rau,Ic,Jc,Ib1,Ib2,'U')
  !Call CGSTAB(aP,aW,aE,aS,aN,b,U,U0,Rau,Ic,Jc,Ib1,Ib2,'U')
#endif
  Call Condiff('V')
#ifdef HYPRE
  Call hypresolve(aP,aW,aE,aS,aN,b,V,V0,Rau,Ic,Jc,Ib1,Ib2,3,'V')
#else
  Call sor(aP,aW,aE,aS,aN,b,V,V0,Rau,Ic,Jc,Ib1,Ib2,'V')
  !Call CGSTAB(aP,aW,aE,aS,aN,b,V,V0,Rau,Ic,Jc,Ib1,Ib2,'V')
#endif
  Call PCEcoe
  rmsm=sum(abs(b))/(Ic*Jc)
#ifdef HYPRE
  Call hypresolve(aP,aW,aE,aS,aN,b,dP,dP,1.0,Ic,Jc,Ib1,Ib2,3,'dP')
#else
  !Call sor(aP,aW,aE,aS,aN,b,dP,dP,1.0,Ic,Jc,Ib1,Ib2,'dP')
  Call CGSTAB(aP,aW,aE,aS,aN,b,dP,dP,1.0,Ic,Jc,Ib1,Ib2,'dP')
#endif
  if(Proctrl=='com') then
   dP=min(max(dP,(1-Po-P)/Rap),(5e+10-Po-P)/Rap)
  end if
  P=P+Rap*dP
  if(Proctrl=='com') then
   rho=rho+Rap*dP*Ma/(R*T)
   Call Densityinterp
  end if
  Call Derivatives('dP')
  Call Velcorrect
  Call Derivatives('U')
  Call Derivatives('V')
  Call Derivatives('P')
  if(Energy=='Y') then
   Call Wallfunc
   Call Condiff('T')
#ifdef HYPRE
   Call hypresolve(aP,aW,aE,aS,aN,b,T,T0,Rae,Ic,Jc,Ib1,Ib2,3,'T')
#else
   Call sor(aP,aW,aE,aS,aN,b,T,T0,Rae,Ic,Jc,Ib1,Ib2,'T')
   !Call CGSTAB(aP,aW,aE,aS,aN,b,T,T0,Rae,Ic,Jc,Ib1,Ib2,'T')
#endif
   if(Proctrl=='com') then
    T=max(min(5e+3,T),1.0)
    rho=(Po+P)*Ma/(R*T)
    Call Densityinterp
   end if
   if(Turmod/='inv') then
    mu=mu0*(T/Ti)**1.5*(Ti+Si)/(T+Si)
   end if
  end if
  if(Turmod=='sa') then
   Call Derivatives('Tn')
   Call Wallfunc
   Call Condiff('Tn')
#ifdef HYPRE
   Call hypresolve(aP,aW,aE,aS,aN,b,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,3,'Tn')
#else
   Call sor(aP,aW,aE,aS,aN,b,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,'Tn')
   !Call CGSTAB(aP,aW,aE,aS,aN,b,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,'Tn')
#endif
   Tn=min(max(Tn,1e-15),1e+5*mu/rho)
  else if(Turmod=='ke') then
   Call Wallfunc
   Call Condiff('Tk')
#ifdef HYPRE
   Call hypresolve(aP,aW,aE,aS,aN,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,3,'Tk')
#else
   Call sor(aP,aW,aE,aS,aN,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
   !Call CGSTAB(aP,aW,aE,aS,aN,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
#endif
   Tk=max(Tk,1e-15)
   Call Condiff('Te')
#ifdef HYPRE
   Call hypresolve(aP,aW,aE,aS,aN,b,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,3,'Te')
#else
   Call sor(aP,aW,aE,aS,aN,b,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,'Te')
   !Call CGSTAB(aP,aW,aE,aS,aN,b,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,'Te')
#endif
   Te=max(Te,(rho*Cu*Tk**2)/(1e+5*mu))
  else if(Turmod=='sst') then
   Call Derivatives('Tk')
   Call Derivatives('Tw')
   Call Wallfunc
   Call Condiff('Tk')
#ifdef HYPRE
   Call hypresolve(aP,aW,aE,aS,aN,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,3,'Tk')
#else
   Call sor(aP,aW,aE,aS,aN,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
   !Call CGSTAB(aP,aW,aE,aS,aN,b,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
#endif
   Tk=max(Tk,1e-15)
   Call Condiff('Tw')
#ifdef HYPRE
   Call hypresolve(aP,aW,aE,aS,aN,b,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,3,'Tw')
#else
   Call sor(aP,aW,aE,aS,aN,b,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,'Tw')
   !Call CGSTAB(aP,aW,aE,aS,aN,b,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,'Tw')
#endif
   Tw=max(Tw,(rho*Tk)/(1e+5*mu))
  end if
  Call Turvis
  Call Postprocess
  if((libmod=='S'.or.libmod=='I').and.mod(s,100)==0) then
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
  if(libmod=='S'.or.libmod=='I') then
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
  if(Turmod=='sa') then
   if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsn<delta) then
    print *,'RANS Equations converge!'
    exit
   end if
   if(rmsm>1e+5.or.rmsu>1e+5.or.rmsv>1e+5.or.rmst>1e+5.or.rmsn>1e+5) then
    print *,'RANS Equations divergence!'
    stop
   end if
  else if(Turmod=='ke') then
   if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsk<delta.and.rmse<delta) then
    print *,'RANS Equations converge!'
    exit
   end if
   if(rmsm>1e+5.or.rmsu>1e+5.or.rmsv>1e+5.or.rmst>1e+5.or.rmsk>1e+5.or.rmse>1e+5) then
    print *,'RANS Equations divergence!'
    stop
   end if
  else if(Turmod=='sst') then
   if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsk<delta.and.rmsw<delta) then
    print *,'RANS Equations converge!'
    exit
   end if
   if(rmsm>1e+5.or.rmsu>1e+5.or.rmsv>1e+5.or.rmst>1e+5.or.rmsk>1e+5.or.rmsw>1e+5) then
    print *,'RANS Equations divergence!'
    stop
   end if
  else if(Turmod=='lam'.or.Turmod=='inv') then
   if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta) then
    print *,'RANS Equations converge!'
    exit
   end if
   if(rmsm>1e+5.or.rmsu>1e+5.or.rmsv>1e+5.or.rmst>1e+5) then
    print *,'RANS Equations divergence!'
    stop
   end if
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
if(time1(6)-time0(6)>=0) then
 if(time1(7)-time0(7)>=0) then
  print *,'Iteration time:',time1(6)-time0(6),'m',time1(7)-time0(7),'s'
 else
  print *,'Iteration time:',time1(6)-time0(6)-1,'m',60+time1(7)-time0(7),'s'
 end if
else
 if(time1(7)-time0(7)>=0) then
  print *,'Iteration time:',60+time1(6)-time0(6),'m',time1(7)-time0(7),'s'
 else
  print *,'Iteration time:',60+time1(6)-time0(6)-1,'m',60+time1(7)-time0(7),'s'
 end if
end if
if(s==maxs+1) print *,'Reach maximum iteration steps!'
if(libmod=='S'.or.libmod=='I') then
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
  deallocate(Xwp,Ywp,Xwp0,Ywp0)
 else
  deallocate(Xwd,Ywd,Xwu,Ywu)
 end if
end if

end Subroutine Aero2D

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
deallocate(U0,V0,T0,U,V,T,rho,mu,mut,P,dP,Pr,Pc,auP,auNB,aP,aW,aE,aS,aN,b,Xg,Yg,Xc,Yc,Xga,Xgk,Yga,Ygk,dk,da,Jg,a1,y1,b1,&
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
