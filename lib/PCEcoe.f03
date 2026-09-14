Subroutine PCEcoe
use Aero2DCOM
implicit none
integer i,j
real(8) Up,Vp,Unpk,Vnpa,Unk0,Vna0,Unf0,Vnf0,ww,we,ws,wn,dwk,dwa,df0,Pak,Pka,cor,cc,noc,Xf,Yf
real(8) aP,aW,aE,aS,aN
real(8) du(Ic,Jc),dv(Ic,Jc),Unp(Ic,Jc),Vnp(Ic,Jc),d0(Ic,Jc),Un0(Ic,Jc),Vn0(Ic,Jc)
logical(1) isSimp,isSimpC,isCom,isPseudo,isInOut

noc=1.d+0
cc=1.d+0

isSimp=solctrlFlag==SIMPLE
isSimpC=solctrlFlag==SIMPLEC
isCom=ProctrlFlag==COM
isPseudo=EvolveFlag==PSEUDO
isInOut=FstypeFlag==VINPOUT

!$OMP PARALLEL
if(isSimp) then
  !$OMP DO PRIVATE(i,Up,Vp)
  DO j=1,Jc-1
    DO i=Is,Ie
     du(i,j)=Rau*a1(i,j)*Jg(i,j)*dy/auM(1,i,j)
     dv(i,j)=Rau*y1(i,j)*Jg(i,j)*dx/auM(1,i,j)
     dwno(i,j)=Rau*b1(i,j)*Jg(i,j)/auM(1,i,j)
     Up=U(i,j)+Rau*Px(i,j)*Jg(i,j)*dx*dy/auM(1,i,j)
     Vp=V(i,j)+Rau*Py(i,j)*Jg(i,j)*dx*dy/auM(1,i,j)
     Unp(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
     Vnp(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
     if(isPseudo) then
      Up=Rau*rho0(i,j)*U0(i,j)*Jg(i,j)*dx*dy/deltat/auM(1,i,j)
      Vp=Rau*rho0(i,j)*V0(i,j)*Jg(i,j)*dx*dy/deltat/auM(1,i,j)
      Un0(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
      Vn0(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
      d0(i,j)=Rau*rho0(i,j)*Jg(i,j)*dx*dy/deltat/auM(1,i,j)
     end if
    end DO
  end DO
  !$OMP END DO
else if(isSimpC) then
  !$OMP DO PRIVATE(i,Up,Vp)
  DO j=1,Jc-1
    DO i=Is,Ie
     du(i,j)=Rau*a1(i,j)*Jg(i,j)*dy/(auM(1,i,j)-Rau*auM(2,i,j))
     dv(i,j)=Rau*y1(i,j)*Jg(i,j)*dx/(auM(1,i,j)-Rau*auM(2,i,j))
     dwno(i,j)=Rau*b1(i,j)*Jg(i,j)/(auM(1,i,j)-Rau*auM(2,i,j))
     Up=U(i,j)+Rau*Px(i,j)*Jg(i,j)*dx*dy/(auM(1,i,j)-Rau*auM(2,i,j))
     Vp=V(i,j)+Rau*Py(i,j)*Jg(i,j)*dx*dy/(auM(1,i,j)-Rau*auM(2,i,j))
     Unp(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
     Vnp(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
     if(isPseudo) then
      Up=Rau*rho0(i,j)*U0(i,j)*Jg(i,j)*dx*dy/deltat/(auM(1,i,j)-Rau*auM(2,i,j))
      Vp=Rau*rho0(i,j)*V0(i,j)*Jg(i,j)*dx*dy/deltat/(auM(1,i,j)-Rau*auM(2,i,j))
      Un0(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
      Vn0(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
      d0(i,j)=Rau*rho0(i,j)*Jg(i,j)*dx*dy/deltat/(auM(1,i,j)-Rau*auM(2,i,j))
     end if
    end DO
  end DO
  !$OMP END DO
end if
!$OMP DO PRIVATE(i,Unpk,dwk,Pak,cor)
DO j=1,Jc-1
 DO i=Is,Ie+1
  if(i==1) then
   duk(i,j)=interpl(du(i,j),du(Ic,j),dk(i,j),dk(Ic,j))
   dwk=interpl(dwno(i,j)*dy,dwno(Ic,j)*dy,dk(i,j),dk(Ic,j))
   Unpk=interpl(Unp(i,j),Unp(Ic,j),dk(i,j),dk(Ic,j))
   if(isPseudo) then
    df0=interpl(d0(i,j),d0(Ic,j),dk(i,j),dk(Ic,j))
    Unf0=interpl(Un0(i,j),Un0(Ic,j),dk(i,j),dk(Ic,j))
   end if
  else if(Is>1.and.i==2) then
   if(isInOut) then
    duk(i,j)=interpl(du(i,j),du(i,j),dk(i,j),dk(i-1,j))
    dwk=interpl(dwno(i,j)*dy,dwno(i,j)*dy,dk(i,j),dk(i-1,j))
    Unpk=interpl(Unp(i,j),Unp(i,j),dk(i,j),dk(i-1,j))
    if(isPseudo) then
     df0=interpl(d0(i,j),d0(i,j),dk(i,j),dk(i-1,j))
     Unf0=interpl(Un0(i,j),Un0(i,j),dk(i,j),dk(i-1,j))
    end if
   else
    duk(i,j)=interpl(du(i,j),0.0,dk(i,j),dk(i-1,j))
    dwk=interpl(dwno(i,j)*dy,0.0,dk(i,j),dk(i-1,j))
    Unpk=interpl(Unp(i,j),Un(i-1,j),dk(i,j),dk(i-1,j))
    if(isPseudo) then
     df0=interpl(d0(i,j),0.0,dk(i,j),dk(i-1,j))
     Unf0=interpl(Un0(i,j),0.0,dk(i,j),dk(i-1,j))
    end if
   end if
  else if(Is>1.and.i==Ic) then
   if(isInOut) then
    duk(i,j)=interpl(du(i-1,j),du(i-1,j),dk(i,j),dk(i-1,j))
    dwk=interpl(dwno(i-1,j)*dy,dwno(i-1,j)*dy,dk(i,j),dk(i-1,j))
    Unpk=interpl(Unp(i-1,j),Unp(i-1,j),dk(i,j),dk(i-1,j))
    if(isPseudo) then
     df0=interpl(d0(i-1,j),d0(i-1,j),dk(i,j),dk(i-1,j))
     Unf0=interpl(Un0(i-1,j),Un0(i-1,j),dk(i,j),dk(i-1,j))
    end if
   else
    duk(i,j)=interpl(0.0,du(i-1,j),dk(i,j),dk(i-1,j))
    dwk=interpl(0.0,dwno(i-1,j)*dy,dk(i,j),dk(i-1,j))
    Unpk=interpl(Un(i,j),Unp(i-1,j),dk(i,j),dk(i-1,j))
    if(isPseudo) then
     df0=interpl(0.0,d0(i-1,j),dk(i,j),dk(i-1,j))
     Unf0=interpl(0.0,Un0(i-1,j),dk(i,j),dk(i-1,j))
    end if
   end if
  else if(i==Ip) then
   duk(i,j)=interpl(du(1,j),du(i-1,j),dk(1,j),dk(i-1,j))
   dwk=interpl(dwno(1,j)*dy,dwno(i-1,j)*dy,dk(1,j),dk(i-1,j))
   Unpk=interpl(Unp(1,j),Unp(i-1,j),dk(1,j),dk(i-1,j))
   if(isPseudo) then
    df0=interpl(d0(1,j),d0(i-1,j),dk(1,j),dk(i-1,j))
    Unf0=interpl(Un0(1,j),Un0(i-1,j),dk(1,j),dk(i-1,j))
   end if
  else
   duk(i,j)=interpl(du(i,j),du(i-1,j),dk(i,j),dk(i-1,j))
   dwk=interpl(dwno(i,j)*dy,dwno(i-1,j)*dy,dk(i,j),dk(i-1,j))
   Unpk=interpl(Unp(i,j),Unp(i-1,j),dk(i,j),dk(i-1,j))
   if(isPseudo) then
    df0=interpl(d0(i,j),d0(i-1,j),dk(i,j),dk(i-1,j))
    Unf0=interpl(Un0(i,j),Un0(i-1,j),dk(i,j),dk(i-1,j))
   end if
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
  Unk0=Unk(i,j)
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
  if(isPseudo) Unk(i,j)=Unk(i,j)-Unf0+df0*Unk0
 end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i,Vnpa,dwa,Pka,cor)
DO j=1,Jc
 DO i=Is,Ie
  Vna0=Vna(i,j)
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   dva(i,j)=interpl(dv(i,j),dv(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
   dwa=interpl(dwno(i,j)*dx,dwno(Ic+1-i,j)*dx,da(i,j),da(Ic+1-i,j))
   Pka=(P(Ic-i,j)+P(i+1,j)-P(Ic+2-i,j)-P(i-1,j))/(4*dx)
   Vnpa=interpl(Vnp(i,j),-Vnp(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
   cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),-Vn(Ic+1-i,j),da(i,j),da(Ic+1-i,j)))
   Vna(i,j)=Vnpa+dva(i,j)*(P(Ic+1-i,j)-P(i,j))+cc*cor+noc*dwa*Pka
   if(isPseudo) then
    df0=interpl(d0(i,j),d0(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
    Vnf0=interpl(Vn0(i,j),-Vn0(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
    Vna(i,j)=Vna(i,j)-Vnf0+df0*Vna0
   end if
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
    if(Vn(i,j)<0.0) then
     dva(i,j)=0.0
     dwa=0.0
     Xf=interpl(Xgk(i,j-1),Xgk(i,j),da(i,j-1),da(i,j))
     Yf=interpl(Ygk(i,j-1),Ygk(i,j),da(i,j-1),da(i,j))
     Vnpa=V(i,j)*Xf-U(i,j)*Yf
     cor=0.0
     if(isPseudo) then
      df0=0.0
      Vnf0=0.0
     end if
    else
     dva(i,j)=interpl(dv(i,j-1),dv(i,j-1),da(i,j),da(i,j-1))
     dwa=interpl(dwno(i,j-1)*dx,dwno(i,j-1)*dx,da(i,j),da(i,j-1))
     Vnpa=interpl(Vnp(i,j-1),Vnp(i,j-1),da(i,j),da(i,j-1))
     cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1)))
     if(isPseudo) then
      df0=interpl(d0(i,j-1),d0(i,j-1),da(i,j),da(i,j-1))
      Vnf0=interpl(Vn0(i,j-1),Vn0(i,j-1),da(i,j),da(i,j-1))
     end if
    end if
   else
    dva(i,j)=interpl(0.0,dv(i,j-1),da(i,j),da(i,j-1))
    dwa=interpl(0.0,dwno(i,j-1)*dx,da(i,j),da(i,j-1))
    Vnpa=interpl(Vn(i,j),Vnp(i,j-1),da(i,j),da(i,j-1))
    cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1)))
    if(isPseudo) then
     df0=interpl(0.0,d0(i,j-1),da(i,j),da(i,j-1))
     Vnf0=interpl(0.0,Vn0(i,j-1),da(i,j),da(i,j-1))
    end if
   end if
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))+cc*cor+noc*dwa*Pka
   if(isPseudo) Vna(i,j)=Vna(i,j)-Vnf0+df0*Vna0
  else
   dva(i,j)=interpl(dv(i,j),dv(i,j-1),da(i,j),da(i,j-1))
   dwa=interpl(dwno(i,j)*dx,dwno(i,j-1)*dx,da(i,j),da(i,j-1))
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
   if(isPseudo) then
    df0=interpl(d0(i,j),d0(i,j-1),da(i,j),da(i,j-1))
    Vnf0=interpl(Vn0(i,j),Vn0(i,j-1),da(i,j),da(i,j-1))
    Vna(i,j)=Vna(i,j)-Vnf0+df0*Vna0
   end if
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
    if(isInOut.and.j==Jc-1.and.Vn(i,j+1)<0.0) then
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
!$OMP WORKSHARE
rmsm=sum(abs(b))/(Ic*Jc)
!$OMP END WORKSHARE
if(isInOut) then
 !$OMP WORKSHARE
 b(:,Jc)=Vn(:,Jc)
 !$OMP END WORKSHARE
end if
!$OMP END PARALLEL
end Subroutine PCEcoe