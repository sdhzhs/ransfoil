Subroutine FluxI
use Aero2DCOM
implicit none
integer i,j
real(8) Uf,Vf

DO j=1,Jc-1
 DO i=Is,Ie+1
  if(i==1) then
   Uf=interpl(U(Ic,j),U(i,j),dkw(i,j))
   Vf=interpl(V(Ic,j),V(i,j),dkw(i,j))
  else if(i==Ip) then
   Uf=interpl(U(i-1,j),U(1,j),dkw(i,j))
   Vf=interpl(V(i-1,j),V(1,j),dkw(i,j))
  else
   Uf=interpl(U(i-1,j),U(i,j),dkw(i,j))
   Vf=interpl(V(i-1,j),V(i,j),dkw(i,j))
  end if
  Unk(i,j)=Uf*Xfk(i,j)+Vf*Yfk(i,j)
 end DO
end DO
DO j=1,Jc
 DO i=Is,Ie
  if(j==1.and.i>=Ib1.and.i<=Ib2) then
   Uf=0
   Vf=0
  else if(j==1) then
   Uf=interpl(U(Ic+1-i,j),U(i,j),daw(i,j))
   Vf=interpl(V(Ic+1-i,j),V(i,j),daw(i,j))
  else
   Uf=interpl(U(i,j-1),U(i,j),daw(i,j))
   Vf=interpl(V(i,j-1),V(i,j),daw(i,j))
  end if
  Vna(i,j)=Uf*Xfa(i,j)+Vf*Yfa(i,j)
 end DO
end DO
end Subroutine FluxI

Subroutine FluxRC
use Aero2DCOM
implicit none
integer i,j
real(8) Uf,Vf,Sf,Unpk,Vnpa,noc,dncx,dncy,Tfx,Tfy,Pno
real(8) Up(Ic,Jc),Vp(Ic,Jc)

noc=1.d+0

DO j=1,Jc-1
  DO i=Is,Ie
   Up(i,j)=U(i,j)+Rau*Px(i,j)*Vol(i,j)/auM(1,i,j)
   Vp(i,j)=V(i,j)+Rau*Py(i,j)*Vol(i,j)/auM(1,i,j)
  end DO
end DO

DO j=1,Jc-1
 DO i=Is,Ie+1
  if(i==1) then
   Uf=interpl(Up(Ic,j),Up(i,j),dkw(i,j))
   Vf=interpl(Vp(Ic,j),Vp(i,j),dkw(i,j))
  else if(Is>1.and.i==2) then
   Uf=interpl(U(i-1,j),Up(i,j),dkw(i,j))
   Vf=interpl(V(i-1,j),Vp(i,j),dkw(i,j))
  else if(Is>1.and.i==Ic) then
   Uf=interpl(Up(i-1,j),U(i,j),dkw(i,j))
   Vf=interpl(Vp(i-1,j),V(i,j),dkw(i,j))
  else if(i==Ip) then
   Uf=interpl(Up(i-1,j),Up(1,j),dkw(i,j))
   Vf=interpl(Vp(i-1,j),Vp(1,j),dkw(i,j))
  else
   Uf=interpl(Up(i-1,j),Up(i,j),dkw(i,j))
   Vf=interpl(Vp(i-1,j),Vp(i,j),dkw(i,j))
  end if
  Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
  Unpk=Uf*Xfk(i,j)+Vf*Yfk(i,j)
  if(i==1) then
   Uf=interpl(Px(Ic,j),Px(i,j),dkw(i,j))
   Vf=interpl(Py(Ic,j),Py(i,j),dkw(i,j))
   dncx=Xc(i,j)-Xc(Ic,j)
   dncy=Yc(i,j)-Yc(Ic,j)
   Tfx=-(Xfk(i,j)-Sf*dncx/dkd(i,j))
   Tfy=-(Yfk(i,j)-Sf*dncy/dkd(i,j))
   Pno=duk(i,j)*(Uf*Tfx+Vf*Tfy)
   Unk(i,j)=Unpk+duk(i,j)*(P(Ic,j)-P(i,j))*Sf/dkd(i,j)+noc*Pno
  else if(i==Ip) then
   Uf=interpl(Px(i-1,j),Px(1,j),dkw(i,j))
   Vf=interpl(Py(i-1,j),Py(1,j),dkw(i,j))
   dncx=Xc(1,j)-Xc(i-1,j)
   dncy=Yc(1,j)-Yc(i-1,j)
   Tfx=-(Xfk(i,j)-Sf*dncx/dkd(i,j))
   Tfy=-(Yfk(i,j)-Sf*dncy/dkd(i,j))
   Pno=duk(i,j)*(Uf*Tfx+Vf*Tfy)
   Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(1,j))*Sf/dkd(i,j)+noc*Pno
  else
   Uf=interpl(Px(i-1,j),Px(i,j),dkw(i,j))
   Vf=interpl(Py(i-1,j),Py(i,j),dkw(i,j))
   dncx=Xc(i,j)-Xc(i-1,j)
   dncy=Yc(i,j)-Yc(i-1,j)
   Tfx=-(Xfk(i,j)-Sf*dncx/dkd(i,j))
   Tfy=-(Yfk(i,j)-Sf*dncy/dkd(i,j))
   Pno=duk(i,j)*(Uf*Tfx+Vf*Tfy)
   Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(i,j))*Sf/dkd(i,j)+noc*Pno
  end if
 end DO
end DO

DO j=1,Jc
 DO i=Is,Ie
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   Uf=interpl(Up(Ic+1-i,j),Up(i,j),daw(i,j))
   Vf=interpl(Vp(Ic+1-i,j),Vp(i,j),daw(i,j))
  else if(j==1) then
   Uf=0
   Vf=0
  else if(j==Jc) then
   Uf=interpl(Up(i,j-1),U(i,j),daw(i,j))
   Vf=interpl(Vp(i,j-1),V(i,j),daw(i,j))
  else
   Uf=interpl(Up(i,j-1),Up(i,j),daw(i,j))
   Vf=interpl(Vp(i,j-1),Vp(i,j),daw(i,j))
  end if
  Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
  Vnpa=Uf*Xfa(i,j)+Vf*Yfa(i,j)
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   Uf=interpl(Px(Ic+1-i,j),Px(i,j),daw(i,j))
   Vf=interpl(Py(Ic+1-i,j),Py(i,j),daw(i,j))
   dncx=Xc(i,j)-Xc(Ic+1-i,j)
   dncy=Yc(i,j)-Yc(Ic+1-i,j)
   Tfx=-(Xfa(i,j)-Sf*dncx/dad(i,j))
   Tfy=-(Yfa(i,j)-Sf*dncy/dad(i,j))
   Pno=dva(i,j)*(Uf*Tfx+Vf*Tfy)
   Vna(i,j)=Vnpa+dva(i,j)*(P(Ic+1-i,j)-P(i,j))*Sf/dad(i,j)+noc*Pno
  else if(j==1) then
   Vna(i,j)=0
  else
   Uf=interpl(Px(i,j-1),Px(i,j),daw(i,j))
   Vf=interpl(Py(i,j-1),Py(i,j),daw(i,j))
   dncx=Xc(i,j)-Xc(i,j-1)
   dncy=Yc(i,j)-Yc(i,j-1)
   Tfx=-(Xfa(i,j)-Sf*dncx/dad(i,j))
   Tfy=-(Yfa(i,j)-Sf*dncy/dad(i,j))
   Pno=dva(i,j)*(Uf*Tfx+Vf*Tfy) 
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))*Sf/dad(i,j)+noc*Pno
  end if
 end DO
end DO
end Subroutine FluxRC