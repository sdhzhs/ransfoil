Subroutine PCEcoe
use Aero2DCOM
implicit none
integer i,j
real(8) Up,Vp,Unpk,Vnpa,ww,we,ws,wn,dwk,dwa,Pak,Pka,cor,cc,noc,Xf,Yf
real(8) aP,aW,aE,aS,aN
real(8),external::interpl
real(8) du(Ic,Jc),dv(Ic,Jc),dw(Ic,Jc),Unp(Ic,Jc),Vnp(Ic,Jc)
logical(1) isSimp,isSimpC,isCom,isInOut
noc=1.d+0
cc=1.d+0
isSimp=solctrl=='SIMPLE'
isSimpC=solctrl=='SIMPLEC'
isCom=Proctrl=='com'
isInOut=Fstype=='vinpout'

!$OMP PARALLEL
if(isSimp) then
  !$OMP DO PRIVATE(i,Up,Vp)
  DO j=1,Jc-1
    DO i=Is,Ie
     du(i,j)=Rau*a1(i,j)*Jg(i,j)*dy/auP(i,j)
     dv(i,j)=Rau*y1(i,j)*Jg(i,j)*dx/auP(i,j)
     dw(i,j)=Rau*b1(i,j)*Jg(i,j)/auP(i,j)
     Up=U(i,j)+Rau*Px(i,j)*Jg(i,j)*dx*dy/auP(i,j)
     Vp=V(i,j)+Rau*Py(i,j)*Jg(i,j)*dx*dy/auP(i,j)
     Unp(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
     Vnp(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
    end DO
  end DO
  !$OMP END DO
else if(isSimpC) then
  !$OMP DO PRIVATE(i,Up,Vp)
  DO j=1,Jc-1
    DO i=Is,Ie
     du(i,j)=Rau*a1(i,j)*Jg(i,j)*dy/(auP(i,j)-Rau*auNB(i,j))
     dv(i,j)=Rau*y1(i,j)*Jg(i,j)*dx/(auP(i,j)-Rau*auNB(i,j))
     dw(i,j)=Rau*b1(i,j)*Jg(i,j)/(auP(i,j)-Rau*auNB(i,j))
     Up=U(i,j)+Rau*Px(i,j)*Jg(i,j)*dx*dy/(auP(i,j)-Rau*auNB(i,j))
     Vp=V(i,j)+Rau*Py(i,j)*Jg(i,j)*dx*dy/(auP(i,j)-Rau*auNB(i,j))
     Unp(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
     Vnp(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
    end DO
  end DO
  !$OMP END DO
end if
!$OMP DO PRIVATE(i,Unpk,dwk,Pak,cor)
DO j=1,Jc-1
 DO i=Is,Ie+1
  if(i==1) then
   duk(i,j)=interpl(du(i,j),du(Ic,j),dk(i,j),dk(Ic,j))
   dwk=interpl(dw(i,j)*dy,dw(Ic,j)*dy,dk(i,j),dk(Ic,j))
   Unpk=interpl(Unp(i,j),Unp(Ic,j),dk(i,j),dk(Ic,j))
  else if(Is>1.and.i==2) then
   if(isInOut) then
    duk(i,j)=interpl(du(i,j),du(i,j),dk(i,j),dk(i-1,j))
    dwk=interpl(dw(i,j)*dy,dw(i,j)*dy,dk(i,j),dk(i-1,j))
    Unpk=interpl(Unp(i,j),Unp(i,j),dk(i,j),dk(i-1,j))
   else
    duk(i,j)=interpl(du(i,j),0.0,dk(i,j),dk(i-1,j))
    dwk=interpl(dw(i,j)*dy,0.0,dk(i,j),dk(i-1,j))
    Unpk=interpl(Unp(i,j),Un(i-1,j),dk(i,j),dk(i-1,j))
   end if
  else if(Is>1.and.i==Ic) then
   if(isInOut) then
    duk(i,j)=interpl(du(i-1,j),du(i-1,j),dk(i,j),dk(i-1,j))
    dwk=interpl(dw(i-1,j)*dy,dw(i-1,j)*dy,dk(i,j),dk(i-1,j))
    Unpk=interpl(Unp(i-1,j),Unp(i-1,j),dk(i,j),dk(i-1,j))
   else
    duk(i,j)=interpl(0.0,du(i-1,j),dk(i,j),dk(i-1,j))
    dwk=interpl(0.0,dw(i-1,j)*dy,dk(i,j),dk(i-1,j))
    Unpk=interpl(Unp(i-1,j),Un(i,j),dk(i-1,j),dk(i,j))
   end if
  else if(i==Ip) then
   duk(i,j)=interpl(du(1,j),du(i-1,j),dk(1,j),dk(i-1,j))
   dwk=interpl(dw(1,j)*dy,dw(i-1,j)*dy,dk(1,j),dk(i-1,j))
   Unpk=interpl(Unp(1,j),Unp(i-1,j),dk(1,j),dk(i-1,j))
  else
   duk(i,j)=interpl(du(i,j),du(i-1,j),dk(i,j),dk(i-1,j))
   dwk=interpl(dw(i,j)*dy,dw(i-1,j)*dy,dk(i,j),dk(i-1,j))
   Unpk=interpl(Unp(i,j),Unp(i-1,j),dk(i,j),dk(i-1,j))
  end if
  if(j==1.and.(i<Ib1.or.i>Ib2+1)) then
   Pak=(P(i-1,j+1)+P(i,j+1)-P(Ic+2-i,j)-P(Ic+1-i,j))/(4*dy)
  else if(j==1) then
   if(i==1) then
    Pak=(P(Ic,j+1)+P(i,j+1)-P(Ic,j)-P(i,j))/(4*dy)
   else if(i==Ip) then
    Pak=(P(i-1,j+1)+P(1,j+1)-P(i-1,j)-P(1,j))/(4*dy)
   else
    Pak=(P(i-1,j+1)+P(i,j+1)-P(i-1,j)-P(i,j))/(4*dy)
   end if
  else
   if(i==1) then
    Pak=(P(Ic,j+1)+P(i,j+1)-P(Ic,j-1)-P(i,j-1))/(4*dy)
   else if(i==Ip) then
    Pak=(P(i-1,j+1)+P(1,j+1)-P(i-1,j-1)-P(1,j-1))/(4*dy)
   else
    Pak=(P(i-1,j+1)+P(i,j+1)-P(i-1,j-1)-P(i,j-1))/(4*dy)
   end if
  end if
  if(i==1) then
   cor=(1-Rau)*(Unk(i,j)-interpl(Un(i,j),Un(Ic,j),dk(i,j),dk(Ic,j)))
   Unk(i,j)=Unpk+duk(i,j)*(P(Ic,j)-P(i,j))+cc*cor+noc*dwk*Pak
  else if(i==Ip) then
   cor=(1-Rau)*(Unk(i,j)-interpl(Un(1,j),Un(i-1,j),dk(1,j),dk(i-1,j)))
   Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(1,j))+cc*cor+noc*dwk*Pak
  else
   cor=(1-Rau)*(Unk(i,j)-interpl(Un(i,j),Un(i-1,j),dk(i,j),dk(i-1,j)))
   Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(i,j))+cc*cor+noc*dwk*Pak
  end if
 end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i,Vnpa,dwa,Pka,cor)
DO j=1,Jc
 DO i=Is,Ie
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   dva(i,j)=interpl(dv(i,j),dv(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
   dwa=interpl(dw(i,j)*dx,dw(Ic+1-i,j)*dx,da(i,j),da(Ic+1-i,j))
   Pka=(P(Ic-i,j)+P(i+1,j)-P(Ic+2-i,j)-P(i-1,j))/(4*dx)
   Vnpa=interpl(Vnp(i,j),-Vnp(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
   cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),-Vn(Ic+1-i,j),da(i,j),da(Ic+1-i,j)))
   Vna(i,j)=Vnpa+dva(i,j)*(P(Ic+1-i,j)-P(i,j))+cc*cor+noc*dwa*Pka
  else if(j==1) then
   dva(i,j)=0
   Vna(i,j)=0
  else if(j==Jc) then
   if(i==1) then
    Pka=(P(i+1,j-1)+P(i+1,j)-P(Ic,j-1)-P(Ic,j))/(4*dx)
   else if(i==Ic) then
    Pka=(P(1,j-1)+P(1,j)-P(i-1,j-1)-P(i-1,j))/(4*dx)
   else
    Pka=(P(i+1,j-1)+P(i+1,j)-P(i-1,j-1)-P(i-1,j))/(4*dx)
   end if
   if(isInOut) then
    dva(i,j)=0.0
    dwa=0.0
    Xf=interpl(Xgk(i,j-1),Xgk(i,j),da(i,j-1),da(i,j))
    Yf=interpl(Ygk(i,j-1),Ygk(i,j),da(i,j-1),da(i,j))
    Vnpa=V(i,j)*Xf-U(i,j)*Yf
    cor=0.0
   else
    dva(i,j)=interpl(0.0,dv(i,j-1),da(i,j),da(i,j-1))
    dwa=interpl(0.0,dw(i,j-1)*dx,da(i,j),da(i,j-1))
    Vnpa=interpl(Vnp(i,j-1),Vn(i,j),da(i,j-1),da(i,j))
    cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1)))
   end if
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))+cc*cor+noc*dwa*Pka
  else
   dva(i,j)=interpl(dv(i,j),dv(i,j-1),da(i,j),da(i,j-1))
   dwa=interpl(dw(i,j)*dx,dw(i,j-1)*dx,da(i,j),da(i,j-1))
   if(i==1) then
    Pka=(P(i+1,j-1)+P(i+1,j)-P(Ic,j-1)-P(Ic,j))/(4*dx)
   else if(i==Ic) then
    Pka=(P(1,j-1)+P(1,j)-P(i-1,j-1)-P(i-1,j))/(4*dx)
   else
    Pka=(P(i+1,j-1)+P(i+1,j)-P(i-1,j-1)-P(i-1,j))/(4*dx)
   end if
   Vnpa=interpl(Vnp(i,j),Vnp(i,j-1),da(i,j),da(i,j-1))
   cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1)))
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))+cc*cor+noc*dwa*Pka
  end if
 end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i,ww,we,ws,wn,aP,aW,aE,aS,aN)
DO j=1,Jc-1
  DO i=Is,Ie
   aE=rhok(i+1,j)*duk(i+1,j)*dy
   aW=rhok(i,j)*duk(i,j)*dy
   aN=rhoa(i,j+1)*dva(i,j+1)*dx
   aS=rhoa(i,j)*dva(i,j)*dx
   aP=aE+aW+aN+aS
   if(isCom) then
    ww=sign(0.5,Unk(i,j))
    we=sign(0.5,Unk(i+1,j))
    ws=sign(0.5,Vna(i,j))
    wn=sign(0.5,Vna(i,j+1))
    if(i==Ic) then
     aE=aE-Rap*(0.5-we)*Unk(i+1,j)*dy/(R*T(1,j)/Ma)
    else
     aE=aE-Rap*(0.5-we)*Unk(i+1,j)*dy/(R*T(i+1,j)/Ma)
    end if
    if(i==1) then
     aW=aW+Rap*(0.5+ww)*Unk(i,j)*dy/(R*T(Ic,j)/Ma)
    else
     aW=aW+Rap*(0.5+ww)*Unk(i,j)*dy/(R*T(i-1,j)/Ma)
    end if
    aN=aN-Rap*(0.5-wn)*Vna(i,j+1)*dx/(R*T(i,j+1)/Ma)
    if(j==1.and.(i>Ib2.or.i<Ib1)) then
     aS=aS+Rap*(0.5+ws)*Vna(i,j)*dx/(R*T(Ic+1-i,j)/Ma)
    else if(j>1) then
     aS=aS+Rap*(0.5+ws)*Vna(i,j)*dx/(R*T(i,j-1)/Ma)
    end if
    aP=aP+Rap*((0.5+we)*Unk(i+1,j)*dy-(0.5-ww)*Unk(i,j)*dy+(0.5+wn)*Vna(i,j+1)*dx-(0.5-ws)*Vna(i,j)*dx)/(R*T(i,j)/Ma)
    if(isInOut.and.j==Jc-1) then
     aN=aN+Rap*(0.5-wn)*Vna(i,j+1)*dx/(R*T(i,j+1)/Ma)
     aP=aP-Rap*(0.5+wn)*Vna(i,j+1)*dx/(R*T(i,j)/Ma)
    end if
   end if
   aM(1,i,j)=aP
   aM(2,i,j)=aW
   aM(3,i,j)=aE
   aM(4,i,j)=aS
   aM(5,i,j)=aN
  end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
b=0
!$OMP END WORKSHARE
!$OMP DO PRIVATE(i)
DO j=1,Jc-1
  DO i=Is,Ie
   b(i,j)=rhok(i,j)*Unk(i,j)*dy-rhok(i+1,j)*Unk(i+1,j)*dy+rhoa(i,j)*Vna(i,j)*dx-rhoa(i,j+1)*Vna(i,j+1)*dx
  end DO
end DO
!$OMP END DO
!$OMP END PARALLEL
end Subroutine PCEcoe
