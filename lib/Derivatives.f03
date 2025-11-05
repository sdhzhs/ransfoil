Subroutine Derivatives(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Fw,Fe,Fs,Fn,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,Jgc
real(8),external:: interpl
real(8) F(Ic,Jc),Fwall(Ib1:Ib2),Fx(Ic,Jc),Fy(Ic,Jc)
character(*) scalar

!$OMP PARALLEL
if(scalar=='U') then
!$OMP WORKSHARE
 F=U
 Fwall=0
!$OMP END WORKSHARE
else if(scalar=='V') then
!$OMP WORKSHARE
 F=V
 Fwall=0
!$OMP END WORKSHARE
else if(scalar=='P') then
!$OMP WORKSHARE
 F=P
 Fwall=P(Ib1:Ib2,1)
!$OMP END WORKSHARE
else if(scalar=='dP') then
!$OMP WORKSHARE
 F=dP
 Fwall=dP(Ib1:Ib2,1)
!$OMP END WORKSHARE
else if(scalar=='rho') then
!$OMP WORKSHARE
 F=rho
 Fwall=1.5*rho(Ib1:Ib2,1)-0.5*rho(Ib1:Ib2,2)
!$OMP END WORKSHARE
else if(scalar=='Tn') then
!$OMP WORKSHARE
 F=Tn
 Fwall=1.5*Tn(Ib1:Ib2,1)-0.5*Tn(Ib1:Ib2,2)
!$OMP END WORKSHARE
else if(scalar=='Tk') then
!$OMP WORKSHARE
 F=Tk
 Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
 !Fwall=Tk(Ib1:Ib2,1)
 !Fwall=0
!$OMP END WORKSHARE
else if(scalar=='Tw') then
!$OMP WORKSHARE
 F=Tw
 Fwall=1.5*Tw(Ib1:Ib2,1)-0.5*Tw(Ib1:Ib2,2)
 !Fwall=Tw(Ib1:Ib2,1)
!$OMP END WORKSHARE
else if(scalar=='mux') then
!$OMP WORKSHARE
 F=(mu+mut)*Ux
 Fwall=1.5*(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Ux(Ib1:Ib2,1)-0.5*(mu(Ib1:Ib2,2)+mut(Ib1:Ib2,2))*Ux(Ib1:Ib2,2)
 !Fwall=(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Ux(Ib1:Ib2,1)
!$OMP END WORKSHARE
else if(scalar=='muy') then
!$OMP WORKSHARE
 F=(mu+mut)*Uy
 Fwall=1.5*(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Uy(Ib1:Ib2,1)-0.5*(mu(Ib1:Ib2,2)+mut(Ib1:Ib2,2))*Uy(Ib1:Ib2,2)
 !Fwall=(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Uy(Ib1:Ib2,1)
!$OMP END WORKSHARE
else if(scalar=='mvx') then
!$OMP WORKSHARE
 F=(mu+mut)*Vx
 Fwall=1.5*(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Vx(Ib1:Ib2,1)-0.5*(mu(Ib1:Ib2,2)+mut(Ib1:Ib2,2))*Vx(Ib1:Ib2,2)
 !Fwall=(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Vx(Ib1:Ib2,1)
!$OMP END WORKSHARE
else if(scalar=='mvy') then
!$OMP WORKSHARE
 F=(mu+mut)*Vy
 Fwall=1.5*(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Vy(Ib1:Ib2,1)-0.5*(mu(Ib1:Ib2,2)+mut(Ib1:Ib2,2))*Vy(Ib1:Ib2,2)
 !Fwall=(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Vy(Ib1:Ib2,1)
!$OMP END WORKSHARE
end if
!$OMP WORKSHARE
Fx=0
Fy=0
!$OMP END WORKSHARE
!$OMP DO PRIVATE(Fw,Fe,Fs,Fn,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,Jgc,i)
DO j=1,Jc-1
 DO i=Is,Ie
  Xgaw=-Yfk(i,j)
  Ygaw=Xfk(i,j)
  Xgae=-Yfk(i+1,j)
  Ygae=Xfk(i+1,j)
  Xgks=Yfa(i,j)
  Ygks=-Xfa(i,j)
  Xgkn=Yfa(i,j+1)
  Ygkn=-Xfa(i,j+1)
  Jgc=Vol(i,j)
  if(i==1) then
   Fw=interpl(F(Ic,j),F(i,j),dkw(i,j))
  else
   Fw=interpl(F(i-1,j),F(i,j),dkw(i,j))
  end if
  if(i==Ic) then
   Fe=interpl(F(i,j),F(1,j),dkw(i+1,j))
  else
   Fe=interpl(F(i,j),F(i+1,j),dkw(i+1,j))
  end if
  Fn=interpl(F(i,j),F(i,j+1),daw(i,j+1))
  if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
   Fs=Fwall(i)
  else if(j==1) then
   Fs=interpl(F(Ic+1-i,j),F(i,j),daw(i,j))
  else
   Fs=interpl(F(i,j-1),F(i,j),daw(i,j))
  end if
  Fx(i,j)=((Fe*Ygae-Fw*Ygaw)-(Fn*Ygkn-Fs*Ygks))/Jgc
  Fy(i,j)=(-(Fe*Xgae-Fw*Xgaw)+(Fn*Xgkn-Fs*Xgks))/Jgc
 end DO
end DO
!$OMP END DO
if(scalar=='U') then
!$OMP WORKSHARE
 Ux=Fx
 Uy=Fy
!$OMP END WORKSHARE
else if(scalar=='V') then
!$OMP WORKSHARE
 Vx=Fx
 Vy=Fy
!$OMP END WORKSHARE
else if(scalar=='P') then
!$OMP WORKSHARE
 Px=Fx
 Py=Fy
!$OMP END WORKSHARE
else if(scalar=='dP') then
!$OMP WORKSHARE
 dPx=Fx
 dPy=Fy
!$OMP END WORKSHARE
else if(scalar=='rho') then
!$OMP WORKSHARE
 rhox=Fx
 rhoy=Fy
!$OMP END WORKSHARE
else if(scalar=='Tn') then
!$OMP WORKSHARE
 Tnx=Fx
 Tny=Fy
!$OMP END WORKSHARE
else if(scalar=='Tk') then
!$OMP WORKSHARE
 Tkx=Fx
 Tky=Fy
!$OMP END WORKSHARE
else if(scalar=='Tw') then
!$OMP WORKSHARE
 Twx=Fx
 Twy=Fy
!$OMP END WORKSHARE
else if(scalar=='mux') then
!$OMP WORKSHARE
 muxx=Fx
 muxy=Fy
!$OMP END WORKSHARE
else if(scalar=='muy') then
!$OMP WORKSHARE
 muyx=Fx
!$OMP END WORKSHARE
else if(scalar=='mvx') then
!$OMP WORKSHARE
 mvxy=Fy
!$OMP END WORKSHARE
else if(scalar=='mvy') then
!$OMP WORKSHARE
 mvyx=Fx
 mvyy=Fy
!$OMP END WORKSHARE
end if
!$OMP END PARALLEL
end Subroutine Derivatives
