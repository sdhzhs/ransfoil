Subroutine Initial
use Aero2DCOM
implicit none
integer i,j,ioerr
character(128) ioerrmsg
real(8),external:: interpl
real(8) Tni0
ca=1006.43
if(Tmptype=='fixed') then
  ka=0.00008*(Ta-Tf)+0.0244
else if(Tmptype=='flux') then
  ka=0.0244
end if
Ui=Vfar*cos(AoA*Pi/180)
Vi=Vfar*sin(AoA*Pi/180)
rhoi=Po*Ma/(R*Ta)
mui=mu0*(Ta/Ti)**1.5*(Ti+Si)/(Ta+Si)
Re=rhoi*Vfar*c/mui
Mach=Vfar/sqrt(gama*R*Ta/Ma)
if(Turmod=='sa') then
 Tni=1
 DO i=1,100
  Tni0=Tni
  Tni=Tni-((rhoi/mui)**3*Tni**4-tvr*(rhoi/mui)**2*Tni**3-mui*tvr*Cnu1**3/rhoi)/(4*(rhoi/mui)**3*Tni**3-&
  3*tvr*(rhoi/mui)**2*Tni**2)
  if(abs(Tni-Tni0)<1e-8) exit
 end DO
else if(Turmod=='ke') then
 Tki=1.5*(Vfar*Itur)**2
 Tei=rhoi*Cu*Tki**2/(mui*tvr)
else if(Turmod=='sst') then
 Tki=1.5*(Vfar*Itur)**2
 Twi=rhoi*Tki/(mui*tvr)
end if
if(Init=='Y') then
 open(unit=1,file=filename(9),form='unformatted',status='old',IOSTAT=ioerr,IOMSG=ioerrmsg)
 if(ioerr>0) stop ioerrmsg
  read(1,iostat=ioerr,IOMSG=ioerrmsg) U
  read(1,iostat=ioerr,IOMSG=ioerrmsg) U
  read(1,iostat=ioerr,IOMSG=ioerrmsg) U
  read(1,iostat=ioerr,IOMSG=ioerrmsg) V
  read(1,iostat=ioerr,IOMSG=ioerrmsg) P
  read(1,iostat=ioerr,IOMSG=ioerrmsg) T
  read(1,iostat=ioerr,IOMSG=ioerrmsg) rho
  if(Turmod=='sa') then
   read(1,iostat=ioerr,IOMSG=ioerrmsg) Tn
  else if(Turmod=='ke') then
   read(1,iostat=ioerr,IOMSG=ioerrmsg) Tk
   read(1,iostat=ioerr,IOMSG=ioerrmsg) Te
  else if(Turmod=='sst') then
   read(1,iostat=ioerr,IOMSG=ioerrmsg) Tk
   read(1,iostat=ioerr,IOMSG=ioerrmsg) Tw
  end if
  if(ioerr>0) STOP ioerrmsg
 close(1)
 print *,'Read initial data completed!'
else if(Init=='N') then
 if(Stag=='N') then
  U=Ui
  V=Vi
  P=0
  T=Ta
 else if(Stag=='Y') then
  U=0.1
  V=0.1
  P=Po*(1+(gama-1)/2*Mach**2)**(gama/(gama-1))-Po
  T=Ta*(1+(gama-1)/2*Mach**2)
 end if
 if(Turmod=='sa') then
  Tn=Tni
 else if(Turmod=='ke') then
  Tk=Tki
  Te=Tei
 else if(Turmod=='sst') then
  Tk=Tki
  Tw=Twi
 end if
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
if(Init=='N') then
 rho=(Po+P)*Ma/(R*T)
end if
ks=ksi
if(Turmod=='inv') then
 mu=0
else
 mu=mu0*(T/Ti)**1.5*(Ti+Si)/(T+Si)
end if
Pr=mu*ca/ka
Pc=9.24*((Pr/Prt)**0.75-1)*(1+0.28*exp(-0.007*Pr/Prt))
dP=0
Call Derivatives('U')
Call Derivatives('V')
Call Derivatives('P')
if(Turmod=='sa'.and.Walltreat=='lr') then
 d=d+0.03*ksi
end if
Call Turvis
Un=U*Yga-V*Xga
Vn=V*Xgk-U*Ygk
DO j=1,Jc-1
  DO i=2,Ic
  Unk(i,j)=interpl(Un(i,j),Un(i-1,j),dk(i,j),dk(i-1,j))
  end DO
end DO
DO j=1,Jc
  DO i=2,Ic-1
  if(j==1.and.i>=Ib1.and.i<=Ib2) then
   Vna(i,j)=0
  else if(j==1) then
   Vna(i,j)=interpl(Vn(i,j),-Vn(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
  else
   Vna(i,j)=interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1))
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
