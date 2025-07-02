Subroutine Derivatives(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Fw,Fe,Fs,Fn,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,Jgc
real(8),external:: interpl
real(8) F(Ic,Jc),Fwall(Ib1:Ib2),Fx(Ic,Jc),Fy(Ic,Jc)
character(*) scalar
if(scalar=='U') then
 F=U
 Fwall=0
else if(scalar=='V') then
 F=V
 Fwall=0
else if(scalar=='P') then
 F=P
 Fwall=P(Ib1:Ib2,1)
else if(scalar=='dP') then
 F=dP
 Fwall=dP(Ib1:Ib2,1)
else if(scalar=='rho') then
 F=rho
 Fwall=1.5*rho(Ib1:Ib2,1)-0.5*rho(Ib1:Ib2,2)
else if(scalar=='Tn') then
 F=Tn
 Fwall=1.5*Tn(Ib1:Ib2,1)-0.5*Tn(Ib1:Ib2,2)
else if(scalar=='Tk') then
 F=Tk
 Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
 !Fwall=Tk(Ib1:Ib2,1)
 !Fwall=0
else if(scalar=='Tw') then
 F=Tw
 Fwall=1.5*Tw(Ib1:Ib2,1)-0.5*Tw(Ib1:Ib2,2)
 !Fwall=Tw(Ib1:Ib2,1)
else if(scalar=='mux') then
 F=(mu+mut)*Ux
 Fwall=1.5*(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Ux(Ib1:Ib2,1)-0.5*(mu(Ib1:Ib2,2)+mut(Ib1:Ib2,2))*Ux(Ib1:Ib2,2)
 !Fwall=(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Ux(Ib1:Ib2,1)
else if(scalar=='muy') then
 F=(mu+mut)*Uy
 Fwall=1.5*(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Uy(Ib1:Ib2,1)-0.5*(mu(Ib1:Ib2,2)+mut(Ib1:Ib2,2))*Uy(Ib1:Ib2,2)
 !Fwall=(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Uy(Ib1:Ib2,1)
else if(scalar=='mvx') then
 F=(mu+mut)*Vx
 Fwall=1.5*(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Vx(Ib1:Ib2,1)-0.5*(mu(Ib1:Ib2,2)+mut(Ib1:Ib2,2))*Vx(Ib1:Ib2,2)
 !Fwall=(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Vx(Ib1:Ib2,1)
else if(scalar=='mvy') then
 F=(mu+mut)*Vy
 Fwall=1.5*(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Vy(Ib1:Ib2,1)-0.5*(mu(Ib1:Ib2,2)+mut(Ib1:Ib2,2))*Vy(Ib1:Ib2,2)
 !Fwall=(mu(Ib1:Ib2,1)+mut(Ib1:Ib2,1))*Vy(Ib1:Ib2,1)
end if
Fx=0
Fy=0
DO j=1,Jc-1
 DO i=Is,Ie
  Xgaw=(Xg(i,j+1)-Xg(i,j))/dy
  Ygaw=(Yg(i,j+1)-Yg(i,j))/dy
  Xgae=(Xg(i+1,j+1)-Xg(i+1,j))/dy
  Ygae=(Yg(i+1,j+1)-Yg(i+1,j))/dy
  Xgks=(Xg(i+1,j)-Xg(i,j))/dx
  Ygks=(Yg(i+1,j)-Yg(i,j))/dx
  Xgkn=(Xg(i+1,j+1)-Xg(i,j+1))/dx
  Ygkn=(Yg(i+1,j+1)-Yg(i,j+1))/dx
  Jgc=Jg(i,j)
  if(i==1) then
   Fw=interpl(F(i,j),F(Ic,j),dk(i,j),dk(Ic,j))
  else
   Fw=interpl(F(i,j),F(i-1,j),dk(i,j),dk(i-1,j))
  end if
  if(i==Ic) then
   Fe=interpl(F(i,j),F(1,j),dk(i,j),dk(1,j))
  else
   Fe=interpl(F(i,j),F(i+1,j),dk(i,j),dk(i+1,j))
  end if
  Fn=interpl(F(i,j),F(i,j+1),da(i,j),da(i,j+1))
  if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
   Fs=Fwall(i)
  else if(j==1) then
   Fs=interpl(F(i,j),F(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
  else
   Fs=interpl(F(i,j),F(i,j-1),da(i,j),da(i,j-1))
  end if
  Fx(i,j)=((Fe*Ygae-Fw*Ygaw)/dx-(Fn*Ygkn-Fs*Ygks)/dy)/Jgc
  Fy(i,j)=(-(Fe*Xgae-Fw*Xgaw)/dx+(Fn*Xgkn-Fs*Xgks)/dy)/Jgc
 end DO
end DO
if(scalar=='U') then
 Ux=Fx
 Uy=Fy
else if(scalar=='V') then
 Vx=Fx
 Vy=Fy
else if(scalar=='P') then
 Px=Fx
 Py=Fy
else if(scalar=='dP') then
 dPx=Fx
 dPy=Fy
else if(scalar=='rho') then
 rhox=Fx
 rhoy=Fy
else if(scalar=='Tn') then
 Tnx=Fx
 Tny=Fy
else if(scalar=='Tk') then
 Tkx=Fx
 Tky=Fy
else if(scalar=='Tw') then
 Twx=Fx
 Twy=Fy
else if(scalar=='mux') then
 muxx=Fx
 muxy=Fy
else if(scalar=='muy') then
 muyx=Fx
else if(scalar=='mvx') then
 mvxy=Fy
else if(scalar=='mvy') then
 mvyx=Fx
 mvyy=Fy
end if
end Subroutine Derivatives
