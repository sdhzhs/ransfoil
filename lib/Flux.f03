Subroutine FluxI
use Aero2DCOM
implicit none
integer i,j
real(8),external::interpl

Un=U*Yga-V*Xga
Vn=V*Xgk-U*Ygk
DO j=1,Jc-1
  DO i=Is,Ie+1
  if(i==1) then
   Unk(i,j)=interpl(Un(i,j),Un(Ic,j),dk(i,j),dk(Ic,j))
  else if(i==Ip) then
   Unk(i,j)=interpl(Un(1,j),Un(i-1,j),dk(1,j),dk(i-1,j))
  else
   Unk(i,j)=interpl(Un(i,j),Un(i-1,j),dk(i,j),dk(i-1,j))
  end if
  end DO
end DO
DO j=1,Jc
  DO i=Is,Ie
  if(j==1.and.i>=Ib1.and.i<=Ib2) then
   Vna(i,j)=0
  else if(j==1) then
   Vna(i,j)=interpl(Vn(i,j),-Vn(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
  else
   Vna(i,j)=interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1))
  end if
  end DO
end DO
end Subroutine FluxI

Subroutine FluxRC
use Aero2DCOM
implicit none
integer i,j
real(8) Up,Vp,Unpk,Vnpa
real(8),external::interpl
real(8) Unp(Ic,Jc),Vnp(Ic,Jc)
logical(1) isCoup

isCoup = solctrlFlag==COUP

Un=U*Yga-V*Xga
Vn=V*Xgk-U*Ygk

DO j=1,Jc-1
 DO i=Is,Ie
  Up=U(i,j)+Rau*Px(i,j)*Jg(i,j)*dx*dy/auM(1,i,j)
  Vp=V(i,j)+Rau*Py(i,j)*Jg(i,j)*dx*dy/auM(1,i,j)
  Unp(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
  Vnp(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
 end DO
end DO
DO j=1,Jc-1
 DO i=Is,Ie+1
  if(isCoup) then
   if(i==1) then
    Unpk=interpl(Unp(i,j),Unp(Ic,j),dk(i,j),dk(i,j))
   else if(Is>1.and.i==2) then
    Unpk=interpl(Unp(i,j),Un(i-1,j),dk(i,j),dk(i,j))
   else if(Is>1.and.i==Ic) then
    Unpk=interpl(Unp(i-1,j),Un(i,j),dk(i,j),dk(i,j))
   else if(i==Ip) then
    Unpk=interpl(Unp(1,j),Unp(i-1,j),dk(i-1,j),dk(i-1,j)) 
   else
    Unpk=interpl(Unp(i,j),Unp(i-1,j),dk(i,j),dk(i,j))
   end if
  else
   if(i==1) then
    Unpk=interpl(Unp(i,j),Unp(Ic,j),dk(i,j),dk(Ic,j))
   else if(Is>1.and.i==2) then
    Unpk=interpl(Unp(i,j),Un(i-1,j),dk(i,j),dk(i-1,j))
   else if(Is>1.and.i==Ic) then
    Unpk=interpl(Unp(i-1,j),Un(i,j),dk(i-1,j),dk(i,j))
   else if(i==Ip) then
    Unpk=interpl(Unp(1,j),Unp(i-1,j),dk(1,j),dk(i-1,j)) 
   else
    Unpk=interpl(Unp(i,j),Unp(i-1,j),dk(i,j),dk(i-1,j))
   end if
  end if
  if(i==1) then
   Unk(i,j)=Unpk+duk(i,j)*(P(Ic,j)-P(i,j))
  else if(i==Ip) then
   Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(1,j))
  else
   Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(i,j))
  end if
 end DO
end DO
DO j=1,Jc
 DO i=Is,Ie
  if(isCoup) then
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    Vnpa=interpl(Vnp(i,j),-Vnp(Ic+1-i,j),da(i,j),da(i,j))
    Vna(i,j)=Vnpa+dva(i,j)*(P(Ic+1-i,j)-P(i,j))
   else if(j==1) then
    Vna(i,j)=0
   else if(j==Jc) then
    Vnpa=interpl(Vnp(i,j-1),Vn(i,j),da(i,j),da(i,j))
    Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))
   else
    Vnpa=interpl(Vnp(i,j),Vnp(i,j-1),da(i,j),da(i,j))
    Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))
   end if
  else
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    Vnpa=interpl(Vnp(i,j),-Vnp(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
    Vna(i,j)=Vnpa+dva(i,j)*(P(Ic+1-i,j)-P(i,j))
   else if(j==1) then
    Vna(i,j)=0
   else if(j==Jc) then
    Vnpa=interpl(Vnp(i,j-1),Vn(i,j),da(i,j-1),da(i,j))
    Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))
   else
    Vnpa=interpl(Vnp(i,j),Vnp(i,j-1),da(i,j),da(i,j-1))
    Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))
   end if
  end if
 end DO
end DO
end Subroutine FluxRC