Subroutine Initial
use Aero2DCOM
implicit none
integer i,j,ioerr
real(8) Tni0
ca=1006.43
ka=0.00008*(Ta-Tf)+0.0244
Ui=Vfar*cos(AoA*Pi/180)
Vi=Vfar*sin(AoA*Pi/180)
roui=Po*Ma/(R*Ta)
miui=miu0*(Ta/Ti)**1.5*(Ti+Si)/(Ta+Si)
Re=roui*Vfar*c/miui
Mach=Vfar/sqrt(gama*R*Ta/Ma)
if(Turmod=='sa') then
Tni=1d-3
DO i=1,100
 Tni0=Tni
 Tni=Tni-((roui/miui)**3*Tni**4-tvr*(roui/miui)**2*Tni**3-miui*tvr*Cniu1**3/roui)/(4*(roui/miui)**3*Tni**3-&
 3*tvr*(roui/miui)**2*Tni**2)
 if(abs(Tni-Tni0)<1d-8) exit
end DO
else if(Turmod=='ke') then
Tki=1.5*(Vfar*Itur)**2
Tei=roui*Cu*Tki**2/(miui*tvr)
else if(Turmod=='sst') then
Tki=1.5*(Vfar*Itur)**2
Twi=roui*Tki/(miui*tvr)
end if
if(Init=='Y') then
open(unit=1,file=filename(9),form='unformatted',status='old')
read(1,iostat=ioerr) U
read(1,iostat=ioerr) U
read(1,iostat=ioerr) U
read(1,iostat=ioerr) V
read(1,iostat=ioerr) P
read(1,iostat=ioerr) T
read(1,iostat=ioerr) rou
if(Turmod=='sa') then
read(1,iostat=ioerr) Tn
else if(Turmod=='ke') then
read(1,iostat=ioerr) Tk
read(1,iostat=ioerr) Te
else if(Turmod=='sst') then
read(1,iostat=ioerr) Tk
read(1,iostat=ioerr) Tw
end if
close(1)
print *,'Read initial data completed!'
else if(Init=='N') then
if(Stag=='N') then
U=Ui+1d-15
V=Vi+1d-15
P=0
T=Ta
else if(Stag=='Y') then
U=0.1
V=0.1
P=Po*(1+(gama-1)/2*Mach**2)**(gama/(gama-1))-Po
T=Ta*(1+(gama-1)/2*Mach**2)
end if
U(:,Jc)=Ui
U(1,:)=Ui
U(Ic,:)=Ui
V(:,Jc)=Vi
V(1,:)=Vi
V(Ic,:)=Vi
P(:,Jc)=0
P(1,:)=0
P(Ic,:)=0
T(:,Jc)=Ta
T(1,:)=Ta
T(Ic,:)=Ta
rou=(Po+P)*Ma/(R*T)
if(Turmod=='sa') then
Tn=Tni
Tn(:,Jc)=Tni
Tn(1,:)=Tni
Tn(Ic,:)=Tni
else if(Turmod=='ke') then
Tk=Tki
Tk(:,Jc)=Tki
Tk(1,:)=Tki
Tk(Ic,:)=Tki
Te=Tei
Te(:,Jc)=Tei
Te(1,:)=Tei
Te(Ic,:)=Tei
else if(Turmod=='sst') then
Tk=Tki
Tk(:,Jc)=Tki
Tk(1,:)=Tki
Tk(Ic,:)=Tki
Tw=Twi
Tw(:,Jc)=Twi
Tw(1,:)=Twi
Tw(Ic,:)=Twi
end if
else if(Init=='A') then
U(:,Jc)=Ui
U(1,:)=Ui
U(Ic,:)=Ui
V(:,Jc)=Vi
V(1,:)=Vi
V(Ic,:)=Vi
P(:,Jc)=0
P(1,:)=0
P(Ic,:)=0
T(:,Jc)=Ta
T(1,:)=Ta
T(Ic,:)=Ta
rou=(Po+P)*Ma/(R*T)
if(Turmod=='sa') then
Tn(:,Jc)=Tni
Tn(1,:)=Tni
Tn(Ic,:)=Tni
else if(Turmod=='ke') then
Tk(:,Jc)=Tki
Tk(1,:)=Tki
Tk(Ic,:)=Tki
Te(:,Jc)=Tei
Te(1,:)=Tei
Te(Ic,:)=Tei
else if(Turmod=='sst') then
Tk(:,Jc)=Tki
Tk(1,:)=Tki
Tk(Ic,:)=Tki
Tw(:,Jc)=Twi
Tw(1,:)=Twi
Tw(Ic,:)=Twi
end if
end if
ks=ksi
if(Turmod=='inv') then
miu=0
else
miu=miu0*(T/Ti)**1.5*(Ti+Si)/(T+Si)
end if
Pr=miu*ca/ka
Pc=9.24*((Pr/Prt)**0.75-1)*(1+0.28*exp(-0.007*Pr/Prt))
Call Derivatives('U')
Call Derivatives('V')
Call Derivatives('P')
if(Turmod=='sa') then
if(Walltreat=='lr') then
d=d+0.03*ksi
Xi=rou*Tn/miu+Cks*ksi/d
else if(Walltreat=='wf') then
Xi=rou*Tn/miu
end if
fniu1=Xi**3/(Xi**3+Cniu1**3)
miut=rou*Tn*fniu1
else if(Turmod=='ke') then
miut=(rou*Cu*Tk**2)/Te
else if(Turmod=='sst') then
St=sqrt(2*(Ux**2+Vy**2)+(Uy+Vx)**2)
fai2=max(2*sqrt(Tk)/(0.09*Tw*d),500*miu/(rou*d**2*Tw))
F2=tanh(fai2**2)
Ret=rou*Tk/(miu*Tw)
alphastar=alphastarf*(0.024+Ret/Rk)/(1+Ret/Rk)
miut=(rou*Tk)/(Tw*max(1./alphastar,St*F2/(alpha1*Tw)))
else if(Turmod=='lam'.or.Turmod=='inv') then
miut=0
end if
dP=0
bp=0
DO j=1,Jc
   DO i=1,Ic
   Un(i,j)=U(i,j)*Yga(i,j)-V(i,j)*Xga(i,j)
   end DO
end DO
DO j=1,Jc
   DO i=1,Ic
   Vn(i,j)=V(i,j)*Xgk(i,j)-U(i,j)*Ygk(i,j)
   end DO
end DO
DO j=1,Jc-1
   DO i=2,Ic-1
   Unw(i,j)=0.5*(Un(i,j)+Un(i-1,j))
   Une(i,j)=0.5*(Un(i,j)+Un(i+1,j))
   Vnn(i,j)=0.5*(Vn(i,j)+Vn(i,j+1))
   if(j==1.and.i>=Ib1.and.i<=Ib2) then
   Vns(i,j)=0
   else if(j==1) then
   Vns(i,j)=0.5*(Vn(i,j)-Vn(Ic+1-i,j))
   else
   Vns(i,j)=0.5*(Vn(i,j)+Vn(i,j-1))
   end if
   end DO
end DO
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
end Subroutine Initial
