Subroutine PCEcoe
use Aero2DCOM
implicit none
integer i,j
real(8) Sf,Uf,Vf,Unpk,Vnpa,ww,we,ws,wn,cor
real(8) aP,aW,aE,aS,aN
real(8),external::interpl
real(8) du(Ic,Jc),dv(Ic,Jc),Up(Ic,Jc),Vp(Ic,Jc)
logical(1) isSimp,isSimpC,isCom

isSimp=solctrl=='SIMPLE'
isSimpC=solctrl=='SIMPLEC'
isCom=Proctrl=='com'

!$OMP PARALLEL
if(isSimp) then
  !$OMP DO PRIVATE(i)
  DO j=1,Jc-1
    DO i=Is,Ie
     du(i,j)=Rau*Vol(i,j)/auP(i,j)
     dv(i,j)=Rau*Vol(i,j)/auP(i,j)
     Up(i,j)=U(i,j)+Rau*Px(i,j)*Vol(i,j)/auP(i,j)
     Vp(i,j)=V(i,j)+Rau*Py(i,j)*Vol(i,j)/auP(i,j)
    end DO
  end DO
  !$OMP END DO
else if(isSimpC) then
  !$OMP DO PRIVATE(i)
  DO j=1,Jc-1
    DO i=Is,Ie
     du(i,j)=Rau*Vol(i,j)/(auP(i,j)-Rau*auNB(i,j))
     dv(i,j)=Rau*Vol(i,j)/(auP(i,j)-Rau*auNB(i,j))
     Up(i,j)=U(i,j)+Rau*Px(i,j)*Vol(i,j)/(auP(i,j)-Rau*auNB(i,j))
     Vp(i,j)=V(i,j)+Rau*Py(i,j)*Vol(i,j)/(auP(i,j)-Rau*auNB(i,j))
    end DO
  end DO
  !$OMP END DO
end if
!$OMP DO PRIVATE(i,Sf,Uf,Vf,Unpk,cor)
DO j=1,Jc-1
 DO i=Is,Ie+1
  if(i==1) then
   duk(i,j)=interpl(du(Ic,j),du(i,j),dkw(i,j))
   Uf=interpl(Up(Ic,j),Up(i,j),dkw(i,j))
   Vf=interpl(Vp(Ic,j),Vp(i,j),dkw(i,j))
  else if(i==2) then
   if(Is>1) then
    duk(i,j)=interpl(0.0,du(i,j),dkw(i,j))
    Uf=interpl(U(i-1,j),Up(i,j),dkw(i,j))
    Vf=interpl(V(i-1,j),Vp(i,j),dkw(i,j))
   else
    duk(i,j)=interpl(du(i-1,j),du(i,j),dkw(i,j))
    Uf=interpl(Up(i-1,j),Up(i,j),dkw(i,j))
    Vf=interpl(Vp(i-1,j),Vp(i,j),dkw(i,j))
   end if
  else if(i==Ic) then
   if(Is>1) then
    duk(i,j)=interpl(du(i-1,j),0.0,dkw(i,j))
    Uf=interpl(Up(i-1,j),U(i,j),dkw(i,j))
    Vf=interpl(Vp(i-1,j),V(i,j),dkw(i,j))
   else
    duk(i,j)=interpl(du(i-1,j),du(i,j),dkw(i,j))
    Uf=interpl(Up(i-1,j),Up(i,j),dkw(i,j))
    Vf=interpl(Vp(i-1,j),Vp(i,j),dkw(i,j))
   end if
  else if(i==Ip) then
   duk(i,j)=interpl(du(i-1,j),du(1,j),dkw(i,j))
   Uf=interpl(Up(i-1,j),Up(1,j),dkw(i,j))
   Vf=interpl(Vp(i-1,j),Vp(1,j),dkw(i,j))
  else
   duk(i,j)=interpl(du(i-1,j),du(i,j),dkw(i,j))
   Uf=interpl(Up(i-1,j),Up(i,j),dkw(i,j))
   Vf=interpl(Vp(i-1,j),Vp(i,j),dkw(i,j))
  end if
  Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
  Unpk=Uf*Xfk(i,j)+Vf*Yfk(i,j)
  if(i==1) then
   Uf=interpl(U(Ic,j),U(i,j),dkw(i,j))
   Vf=interpl(V(Ic,j),V(i,j),dkw(i,j))
   cor=(1-Rau)*(Unk(i,j)-(Uf*Xfk(i,j)+Vf*Yfk(i,j)))
   Unk(i,j)=Unpk+duk(i,j)*(P(Ic,j)-P(i,j))*Sf/dkd(i,j)+cor
  else if(i==Ip) then
   Uf=interpl(U(i-1,j),U(1,j),dkw(i,j))
   Vf=interpl(V(i-1,j),V(1,j),dkw(i,j))
   cor=(1-Rau)*(Unk(i,j)-(Uf*Xfk(i,j)+Vf*Yfk(i,j)))
   Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(1,j))*Sf/dkd(i,j)+cor
  else
   Uf=interpl(U(i-1,j),U(i,j),dkw(i,j))
   Vf=interpl(V(i-1,j),V(i,j),dkw(i,j))
   cor=(1-Rau)*(Unk(i,j)-(Uf*Xfk(i,j)+Vf*Yfk(i,j)))
   Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(i,j))*Sf/dkd(i,j)+cor
  end if
 end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i,Sf,Uf,Vf,Vnpa,cor)
DO j=1,Jc
 DO i=Is,Ie
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   dva(i,j)=interpl(dv(Ic+1-i,j),dv(i,j),daw(i,j))
   Uf=interpl(Up(Ic+1-i,j),Up(i,j),daw(i,j))
   Vf=interpl(Vp(Ic+1-i,j),Vp(i,j),daw(i,j))
  else if(j==1) then
   dva(i,j)=0
   Uf=0
   Vf=0
  else if(j==Jc) then
   dva(i,j)=interpl(dv(i,j-1),0.0,daw(i,j))
   Uf=interpl(Up(i,j-1),U(i,j),daw(i,j))
   Vf=interpl(Vp(i,j-1),V(i,j),daw(i,j))
  else
   dva(i,j)=interpl(dv(i,j-1),dv(i,j),daw(i,j))
   Uf=interpl(Up(i,j-1),Up(i,j),daw(i,j))
   Vf=interpl(Vp(i,j-1),Vp(i,j),daw(i,j))
  end if
  Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
  Vnpa=Uf*Xfa(i,j)+Vf*Yfa(i,j)
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   Uf=interpl(U(Ic+1-i,j),U(i,j),daw(i,j))
   Vf=interpl(V(Ic+1-i,j),V(i,j),daw(i,j))
   cor=(1-Rau)*(Vna(i,j)-(Uf*Xfa(i,j)+Vf*Yfa(i,j)))
   Vna(i,j)=Vnpa+dva(i,j)*(P(Ic+1-i,j)-P(i,j))*Sf/dad(i,j)+cor
  else if(j==1) then
   Vna(i,j)=0
  else
   Uf=interpl(U(i,j-1),U(i,j),daw(i,j))
   Vf=interpl(V(i,j-1),V(i,j),daw(i,j))
   cor=(1-Rau)*(Vna(i,j)-(Uf*Xfa(i,j)+Vf*Yfa(i,j)))
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))*Sf/dad(i,j)+cor
  end if
 end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i,Sf,ww,we,ws,wn,aP,aW,aE,aS,aN)
DO j=1,Jc-1
  DO i=Is,Ie
   Sf=sqrt(Xfk(i+1,j)**2+Yfk(i+1,j)**2)
   aE=rhok(i+1,j)*duk(i+1,j)*Sf/dkd(i+1,j)
   Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
   aW=rhok(i,j)*duk(i,j)*Sf/dkd(i,j)
   Sf=sqrt(Xfa(i,j+1)**2+Yfa(i,j+1)**2)
   aN=rhoa(i,j+1)*dva(i,j+1)*Sf/dad(i,j+1)
   Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
   aS=rhoa(i,j)*dva(i,j)*Sf/dad(i,j)
   aP=aE+aW+aN+aS
   if(isCom) then
    ww=sign(0.5,Unk(i,j))
    we=sign(0.5,Unk(i+1,j))
    ws=sign(0.5,Vna(i,j))
    wn=sign(0.5,Vna(i,j+1))
    if(i==Ic) then
     aE=aE-Rap*(0.5-we)*Unk(i+1,j)/(R*T(1,j)/Ma)
    else
     aE=aE-Rap*(0.5-we)*Unk(i+1,j)/(R*T(i+1,j)/Ma)
    end if
    if(i==1) then
     aW=aW+Rap*(0.5+ww)*Unk(i,j)/(R*T(Ic,j)/Ma)
    else
     aW=aW+Rap*(0.5+ww)*Unk(i,j)/(R*T(i-1,j)/Ma)
    end if
    aN=aN-Rap*(0.5-wn)*Vna(i,j+1)/(R*T(i,j+1)/Ma)
    if(j==1.and.(i>Ib2.or.i<Ib1)) then
     aS=aS+Rap*(0.5+ws)*Vna(i,j)/(R*T(Ic+1-i,j)/Ma)
    else if(j>1) then
     aS=aS+Rap*(0.5+ws)*Vna(i,j)/(R*T(i,j-1)/Ma)
    end if
    aP=aP+Rap*((0.5+we)*Unk(i+1,j)-(0.5-ww)*Unk(i,j)+(0.5+wn)*Vna(i,j+1)-(0.5-ws)*Vna(i,j))/(R*T(i,j)/Ma)
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
   b(i,j)=rhok(i,j)*Unk(i,j)-rhok(i+1,j)*Unk(i+1,j)+rhoa(i,j)*Vna(i,j)-rhoa(i,j+1)*Vna(i,j+1)
  end DO
end DO
!$OMP END DO
!$OMP END PARALLEL
!print *,'Completed assembling the coefficient matrix:','P',minval(aM(1,Is:Ie,1:Jc-1)),maxval(aM(1,Is:Ie,1:Jc-1))
end Subroutine PCEcoe