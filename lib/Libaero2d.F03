Subroutine Aero2D(libmod,prtl,scptname)
#include <config.h>
use Aero2DCOM
implicit none
integer i,s
integer(C_INT),bind(C)::prtl
character(*) libmod,scptname
real(8) time0,time1
if(libmod=='S'.or.libmod=='I') then
Call Readpara(libmod,scptname)
else if(libmod=='C') then
DO i=1,8
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
Ram=1
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
Call CPU_TIME(time0)
DO s=1,maxs
  Call Derivatives('mux')
  Call Derivatives('muy')
  Call Derivatives('mvx')
  Call Derivatives('mvy')
  Call Wallfunc
  Call Condiff('U')
  Call Condiff('V')
#ifdef HYPRE
  Call hypresolve(auP,auW,auE,auS,auN,bu,U,U0,Rau,Ic,Jc,Ib1,Ib2,3,'U')
  Call hypresolve(auP,auW,auE,auS,auN,bv,V,V0,Rau,Ic,Jc,Ib1,Ib2,3,'V')
#else
  Call sor(auP,auW,auE,auS,auN,bu,U,U0,Rau,Ic,Jc,Ib1,Ib2,'U')
  Call sor(auP,auW,auE,auS,auN,bv,V,V0,Rau,Ic,Jc,Ib1,Ib2,'V')
  !Call CGSTAB(auP,auW,auE,auS,auN,bu,U,U0,Rau,Ic,Jc,Ib1,Ib2,'U')
  !Call CGSTAB(auP,auW,auE,auS,auN,bv,V,V0,Rau,Ic,Jc,Ib1,Ib2,'V')
#endif
  Call PCEcoe
  rmsm=sum(abs(bp))/(Ic*Jc)
#ifdef HYPRE
  Call hypresolve(apP,apW,apE,apS,apN,bp,dP,dP,1d+0,Ic,Jc,Ib1,Ib2,3,'dP')
#else
  !Call sor(apP,apW,apE,apS,apN,bp,dP,dP,1d+0,Ic,Jc,Ib1,Ib2,'dP')
  Call CGSTAB(apP,apW,apE,apS,apN,bp,dP,dP,1d+0,Ic,Jc,Ib1,Ib2,'dP')
#endif
  P=P+Rap*dP
  if(Proctrl=='com') then
  rou=rou+Rar*dP*Ma/(R*T)
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
  Call hypresolve(atP,atW,atE,atS,atN,bt,T,T0,Rae,Ic,Jc,Ib1,Ib2,3,'T')
#else
  Call sor(atP,atW,atE,atS,atN,bt,T,T0,Rae,Ic,Jc,Ib1,Ib2,'T')
  !Call CGSTAB(atP,atW,atE,atS,atN,bt,T,T0,Rae,Ic,Jc,Ib1,Ib2,'T')
#endif
  if(Turmod/='inv') then
  miu=miu0*(T/Ti)**1.5*(Ti+Si)/(T+Si)
  end if
  if(Proctrl=='com') then
  rou=(Po+P)*Ma/(R*T)
  Call Densityinterp
  end if
  end if
  if(Turmod=='sa') then
  Call Derivatives('Tn')
  Call Wallfunc
  Call Condiff('Tn')
#ifdef HYPRE
  Call hypresolve(anP,anW,anE,anS,anN,bn,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,3,'Tn')
#else
  Call sor(anP,anW,anE,anS,anN,bn,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,'Tn')
  !Call CGSTAB(anP,anW,anE,anS,anN,bn,Tn,Tn0,Rat,Ic,Jc,Ib1,Ib2,'Tn')
#endif
  if(Walltreat=='lr') then
  Xi=rou*Tn/miu+Cks*ksi/d
  else if(Walltreat=='wf') then
  Xi=rou*Tn/miu
  end if
  fniu1=Xi**3/(Xi**3+Cniu1**3)
  miut=(1-Ram)*miut+Ram*rou*Tn*fniu1
  else if(Turmod=='ke') then
  Call Wallfunc
  Call Condiff('Tk')
  Call Condiff('Te')
#ifdef HYPRE
  Call hypresolve(akP,akW,akE,akS,akN,bk,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,3,'Tk')
  Call hypresolve(aeP,aeW,aeE,aeS,aeN,be,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,3,'Te')
#else
  Call sor(akP,akW,akE,akS,akN,bk,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
  Call sor(aeP,aeW,aeE,aeS,aeN,be,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,'Te')
  !Call CGSTAB(akP,akW,akE,akS,akN,bk,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
  !Call CGSTAB(aeP,aeW,aeE,aeS,aeN,be,Te,Te0,Rat,Ic,Jc,Ib1,Ib2,'Te')
#endif
  miut=(1-Ram)*miut+Ram*(rou*Cu*Tk**2)/Te
  else if(Turmod=='sst') then
  Call Derivatives('Tk')
  Call Derivatives('Tw')
  Call Wallfunc
  Call Condiff('Tk')
  Call Condiff('Tw')
#ifdef HYPRE
  Call hypresolve(akP,akW,akE,akS,akN,bk,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,3,'Tk')
  Call hypresolve(awP,awW,awE,awS,awN,bw,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,3,'Tw')
#else
  Call sor(akP,akW,akE,akS,akN,bk,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
  Call sor(awP,awW,awE,awS,awN,bw,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,'Tw')
  !Call CGSTAB(akP,akW,akE,akS,akN,bk,Tk,Tk0,Rat,Ic,Jc,Ib1,Ib2,'Tk')
  !Call CGSTAB(awP,awW,awE,awS,awN,bw,Tw,Tw0,Rat,Ic,Jc,Ib1,Ib2,'Tw')
#endif
  fai2=max(2*sqrt(Tk)/(0.09*Tw*d),500*miu/(rou*d**2*Tw))
  F2=tanh(fai2**2)
  Ret=rou*Tk/(miu*Tw)
  alphastar=alphastarf*(0.024+Ret/Rk)/(1+Ret/Rk)
  miut=(1-Ram)*miut+Ram*(rou*Tk)/(Tw*max(1./alphastar,St*F2/(alpha1*Tw)))
  end if
  Call Postprocess
  if((libmod=='S'.or.libmod=='I').and.mod(s,200)==0) then
  Call Results
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
  else if(Turmod=='ke') then
  if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsk<delta.and.rmse<delta) then
  print *,'RANS Equations converge!'
  exit
  end if
  else if(Turmod=='sst') then
  if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta.and.rmsk<delta.and.rmsw<delta) then
  print *,'RANS Equations converge!'
  exit
  end if
  else if(Turmod=='lam'.or.Turmod=='inv') then
  if(rmsm<delta.and.rmsu<delta.and.rmsv<delta.and.rmst<delta) then
  print *,'RANS Equations converge!'
  exit
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
Call CPU_TIME(time1)
print *,'Iteration time:',time1-time0,'s'
if(s==maxs+1) print *,'Reach maximum iteration steps!'
if(libmod=='S'.or.libmod=='I') then
close(2)
Call Results
print *,'Save results successfully!'
end if
if(libmod=='S'.or.libmod=='I') then
Call Deallocarray
deallocate(Xwd,Ywd,Xwu,Ywu)
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
crou=C_LOC(rou)
cmiu=C_LOC(miu)
cP=C_LOC(P)
cVx=C_LOC(U)
cVy=C_LOC(V)
cT=C_LOC(T)
cTn=C_LOC(Tn)
cTk=C_LOC(Tk)
cTe=C_LOC(Te)
cTw=C_LOC(Tw)
cmiut=C_LOC(miut)
deallocate(Xwd,Ywd,Xwu,Ywu)
end if

end Subroutine Aero2D

Subroutine Deallocarray
use Aero2DCOM
implicit none
deallocate(U0,V0,T0,U,V,T,rou,miu,miut,P,dP,Pr,Pc,auP,auW,auE,auS,auN,bu,bv,apP,apW,apE,apS,apN,bp,&
Xg,Yg,Xc,Yc,Xga,Xgk,Yga,Ygk,Jg,a1,y1,b1,Un,Vn,Unw,Une,Vns,Vnn,wdu,edu,sdv,ndv,Ux,Uy,Vx,Vy,Px,Py,&
dPx,dPy,muxx,muxy,muyx,mvxy,mvyx,mvyy,ww,we,ws,wn,rouw,roue,rous,roun,d,St)
if(Energy=='Y') then
 deallocate(atP,atW,atE,atS,atN,bt)
end if
if(Turmod=='sa') then
 deallocate(Tn0,Tn,anP,anW,anE,anS,anN,bn,Tnx,Tny,Xi,fniu1)
else if(Turmod=='ke') then
 deallocate(Tk0,Te0,Tk,Te,akP,akW,akE,akS,akN,bk,aeP,aeW,aeE,aeS,aeN,be,roux,rouy)
else if(Turmod=='sst') then
 deallocate(Tk0,Tw0,Tk,Tw,akP,akW,akE,akS,akN,bk,awP,awW,awE,awS,awN,bw,Tkx,Tky,Twx,Twy,fai2,F2,Ret,&
 alphastar,sigmatk,sigmatw)
end if
deallocate(Xw,Yw,Yp,DR,Sw,ks,Q,Yplus,Ystar,ustar,Uplus,Tplus,hcv,Ax,Ay)
end Subroutine Deallocarray
