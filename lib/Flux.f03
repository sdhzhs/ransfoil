Subroutine FluxI
use Aero2DCOM
implicit none
integer i,j

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
real(8) Up,Vp,Unpk,Vnpa,cor,cc,dwk,dwa,noc,Pak,Pka,Psw,Psc,Pwn,Pnc,Pse,Pec,Pwc
real(8) Unp(Ic,Jc),Vnp(Ic,Jc)

noc=1.d+0
cc=0.d+0

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
  if(i==1) then
   Unpk=interpl(Unp(i,j),Unp(Ic,j),dk(i,j),dk(Ic,j))
   dwk=interpl(dwno(i,j)*dy,dwno(Ic,j)*dy,dk(i,j),dk(Ic,j))
  else if(Is>1.and.i==2) then
   Unpk=interpl(Unp(i,j),Un(i-1,j),dk(i,j),dk(i-1,j))
   dwk=interpl(dwno(i,j)*dy,0.0,dk(i,j),dk(i-1,j))
  else if(Is>1.and.i==Ic) then
   Unpk=interpl(Unp(i-1,j),Un(i,j),dk(i-1,j),dk(i,j))
   dwk=interpl(0.0,dwno(i-1,j)*dy,dk(i,j),dk(i-1,j))
  else if(i==Ip) then
   Unpk=interpl(Unp(1,j),Unp(i-1,j),dk(1,j),dk(i-1,j))
   dwk=interpl(dwno(1,j)*dy,dwno(i-1,j)*dy,dk(1,j),dk(i-1,j))
  else
   Unpk=interpl(Unp(i,j),Unp(i-1,j),dk(i,j),dk(i-1,j))
   dwk=interpl(dwno(i,j)*dy,dwno(i-1,j)*dy,dk(i,j),dk(i-1,j))
  end if
  if(i==1) then
   Pwn=P(Ic,j+1)
   if(j==1) then
    Psw=P(Ic,j)
   else
    Psw=P(Ic,j-1)
   end if
  else
   Pwn=P(i-1,j+1)
   if(j==1.and.(i<Ib1.or.i>Ib2+1)) then
    Psw=P(Ic+2-i,j)
   else if(j==1) then
    Psw=P(i-1,j)
   else
    Psw=P(i-1,j-1)
   end if
  end if
  if(i==Ip) then
   Pnc=P(1,j+1)
   if(j==1) then
    Psc=P(1,j)
   else
    Psc=P(1,j-1)
   end if
  else
   Pnc=P(i,j+1)
   if(j==1.and.(i<Ib1.or.i>Ib2+1)) then
    Psc=P(Ic+1-i,j)
   else if(j==1) then
    Psc=P(i,j)
   else
    Psc=P(i,j-1)
   end if
  end if
  Pak=(Pwn+Pnc-Psw-Psc)/(4*dy)
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
DO j=1,Jc
 DO i=Is,Ie
  if(i==1) then
   Pwc=P(Ic,j)
   if(j==1) then
    Psw=P(Ic,j)
   else
    Psw=P(Ic,j-1)
   end if
  else
   Pwc=P(i-1,j)
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
    Psw=P(Ic+2-i,j)
   else if(j==1) then
    Psw=P(i-1,j)
   else
    Psw=P(i-1,j-1)
   end if
  end if
  if(i==Ic) then
   Pec=P(1,j)
   if(j==1) then
    Pse=P(1,j)
   else
    Pse=P(1,j-1)
   end if
  else
   Pec=P(i+1,j)
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
    Pse=P(Ic-i,j)
   else if(j==1) then
    Pse=P(i+1,j)
   else
    Pse=P(i+1,j-1)
   end if
  end if
  Pka=(Pse+Pec-Psw-Pwc)/(4*dx)
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   dwa=interpl(dwno(i,j)*dx,dwno(Ic+1-i,j)*dx,da(i,j),da(Ic+1-i,j))
   Vnpa=interpl(Vnp(i,j),-Vnp(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
   cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),-Vn(Ic+1-i,j),da(i,j),da(Ic+1-i,j)))
   Vna(i,j)=Vnpa+dva(i,j)*(P(Ic+1-i,j)-P(i,j))+cc*cor+noc*dwa*Pka
  else if(j==1) then
   Vna(i,j)=0
  else if(j==Jc) then
   dwa=interpl(0.0,dwno(i,j-1)*dx,da(i,j),da(i,j-1))
   Vnpa=interpl(Vnp(i,j-1),Vn(i,j),da(i,j-1),da(i,j))
   cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1)))
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))+cc*cor+noc*dwa*Pka
  else
   dwa=interpl(dwno(i,j)*dx,dwno(i,j-1)*dx,da(i,j),da(i,j-1))
   Vnpa=interpl(Vnp(i,j),Vnp(i,j-1),da(i,j),da(i,j-1))
   cor=(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1)))
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))+cc*cor+noc*dwa*Pka
  end if
 end DO
end DO
end Subroutine FluxRC