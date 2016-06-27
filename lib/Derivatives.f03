Subroutine Derivatives(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Fw,Fe,Fs,Fn,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn
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
else if(scalar=='rou') then
F=rou
Fwall=1.5*rou(Ib1:Ib2,1)-0.5*rou(Ib1:Ib2,2)
else if(scalar=='Tn') then
F=Tn
Fwall=1.5*Tn(Ib1:Ib2,1)-0.5*Tn(Ib1:Ib2,2)
else if(scalar=='Tk') then
F=Tk
Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
else if(scalar=='Tw') then
F=Tw
Fwall=1.5*Tw(Ib1:Ib2,1)-0.5*Tw(Ib1:Ib2,2)
else if(scalar=='mux') then
F=(miu+miut)*Ux
Fwall=1.5*(miu(Ib1:Ib2,1)+miut(Ib1:Ib2,1))*Ux(Ib1:Ib2,1)-0.5*(miu(Ib1:Ib2,2)+miut(Ib1:Ib2,2))*Ux(Ib1:Ib2,2)
else if(scalar=='muy') then
F=(miu+miut)*Uy
Fwall=1.5*(miu(Ib1:Ib2,1)+miut(Ib1:Ib2,1))*Uy(Ib1:Ib2,1)-0.5*(miu(Ib1:Ib2,2)+miut(Ib1:Ib2,2))*Uy(Ib1:Ib2,2)
else if(scalar=='mvx') then
F=(miu+miut)*Vx
Fwall=1.5*(miu(Ib1:Ib2,1)+miut(Ib1:Ib2,1))*Vx(Ib1:Ib2,1)-0.5*(miu(Ib1:Ib2,2)+miut(Ib1:Ib2,2))*Vx(Ib1:Ib2,2)
else if(scalar=='mvy') then
F=(miu+miut)*Vy
Fwall=1.5*(miu(Ib1:Ib2,1)+miut(Ib1:Ib2,1))*Vy(Ib1:Ib2,1)-0.5*(miu(Ib1:Ib2,2)+miut(Ib1:Ib2,2))*Vy(Ib1:Ib2,2)
end if
Fx=0
Fy=0
DO j=1,Jc-1
 DO i=2,Ic-1
 Fw=0.5*(F(i,j)+F(i-1,j))
 Xgaw=0.5*(Xga(i,j)+Xga(i-1,j))
 Ygaw=0.5*(Yga(i,j)+Yga(i-1,j))
 Fe=0.5*(F(i,j)+F(i+1,j))
 Xgae=0.5*(Xga(i,j)+Xga(i+1,j))
 Ygae=0.5*(Yga(i,j)+Yga(i+1,j))
 Fn=0.5*(F(i,j)+F(i,j+1))
 Xgkn=0.5*(Xgk(i,j)+Xgk(i,j+1))
 Ygkn=0.5*(Ygk(i,j)+Ygk(i,j+1))
 if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
  Fs=Fwall(i)
  Xgks=Xgk(i,j)
  Ygks=Ygk(i,j)
 else if(j==1) then
  Fs=0.5*(F(i,j)+F(Ic+1-i,j))
  Xgks=0.5*(Xgk(i,j)-Xgk(Ic+1-i,j))
  Ygks=0.5*(Ygk(i,j)-Ygk(Ic+1-i,j))
 else
  Fs=0.5*(F(i,j)+F(i,j-1))
  Xgks=0.5*(Xgk(i,j)+Xgk(i,j-1))
  Ygks=0.5*(Ygk(i,j)+Ygk(i,j-1))
 end if
  Fx(i,j)=((Fe*Ygae-Fw*Ygaw)/dx-(Fn*Ygkn-Fs*Ygks)/dy)/Jg(i,j)
  Fy(i,j)=(-(Fe*Xgae-Fw*Xgaw)/dx+(Fn*Xgkn-Fs*Xgks)/dy)/Jg(i,j)
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
else if(scalar=='rou') then
roux=Fx
rouy=Fy
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
