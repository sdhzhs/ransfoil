Subroutine PCEcoe
use Aero2DCOM
implicit none
integer i,j
real(8) Up,Vp,Unpk,Vnpa,ww,we,ws,wn
real(8) aP,aW,aE,aS,aN
real(8),external:: interpl
real(8) du(Ic,Jc),dv(Ic,Jc),Unp(Ic,Jc),Vnp(Ic,Jc)
if(solctrl=='SIMPLE') then
  DO j=1,Jc-1
    DO i=2,Ic-1
     du(i,j)=Rau*(Yga(i,j)**2+Xga(i,j)**2)*dy/auP(i,j)
     dv(i,j)=Rau*(Ygk(i,j)**2+Xgk(i,j)**2)*dx/auP(i,j)
     Up=U(i,j)+Rau*Px(i,j)*Jg(i,j)*dx*dy/auP(i,j)
     Vp=V(i,j)+Rau*Py(i,j)*Jg(i,j)*dx*dy/auP(i,j)
     Unp(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
     Vnp(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
    end DO
  end DO
else if(solctrl=='SIMPLEC') then
  DO j=1,Jc-1
    DO i=2,Ic-1
     du(i,j)=Rau*(Yga(i,j)**2+Xga(i,j)**2)*dy/(auP(i,j)-Rau*auNB(i,j))
     dv(i,j)=Rau*(Ygk(i,j)**2+Xgk(i,j)**2)*dx/(auP(i,j)-Rau*auNB(i,j))
     Up=U(i,j)+Rau*Px(i,j)*Jg(i,j)*dx*dy/(auP(i,j)-Rau*auNB(i,j))
     Vp=V(i,j)+Rau*Py(i,j)*Jg(i,j)*dx*dy/(auP(i,j)-Rau*auNB(i,j))
     Unp(i,j)=Up*Yga(i,j)-Vp*Xga(i,j)
     Vnp(i,j)=Vp*Xgk(i,j)-Up*Ygk(i,j)
    end DO
  end DO
end if
DO j=1,Jc-1
 DO i=2,Ic
  if(i==2) then
   duk(i,j)=interpl(du(i,j),0.0,dk(i,j),dk(i-1,j))
   Unpk=interpl(Unp(i,j),Un(i-1,j),dk(i,j),dk(i-1,j))
  else if(i==Ic) then
   duk(i,j)=interpl(0.0,du(i-1,j),dk(i,j),dk(i-1,j))
   Unpk=interpl(Unp(i-1,j),Un(i,j),dk(i-1,j),dk(i,j))
  else
   duk(i,j)=interpl(du(i,j),du(i-1,j),dk(i,j),dk(i-1,j))
   Unpk=interpl(Unp(i,j),Unp(i-1,j),dk(i,j),dk(i-1,j))
  end if
  Unk(i,j)=Unpk+duk(i,j)*(P(i-1,j)-P(i,j))+(1-Rau)*(Unk(i,j)-interpl(Un(i,j),Un(i-1,j),dk(i,j),dk(i-1,j)))
 end DO
end DO
DO j=1,Jc
 DO i=2,Ic-1
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   dva(i,j)=interpl(dv(i,j),dv(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
   Vnpa=interpl(Vnp(i,j),-Vnp(Ic+1-i,j),da(i,j),da(Ic+1-i,j))
   Vna(i,j)=Vnpa+dva(i,j)*(P(Ic+1-i,j)-P(i,j))+(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),-Vn(Ic+1-i,j),da(i,j),da(Ic+1-i,j)))
  else if(j==1) then
   dva(i,j)=0
   Vnpa=0
   Vna(i,j)=0
  else if(j==Jc) then
   dva(i,j)=interpl(0.0,dv(i,j-1),da(i,j),da(i,j-1))
   Vnpa=interpl(Vnp(i,j-1),Vn(i,j),da(i,j-1),da(i,j))
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))+(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1)))
  else
   dva(i,j)=interpl(dv(i,j),dv(i,j-1),da(i,j),da(i,j-1))
   Vnpa=interpl(Vnp(i,j),Vnp(i,j-1),da(i,j),da(i,j-1))
   Vna(i,j)=Vnpa+dva(i,j)*(P(i,j-1)-P(i,j))+(1-Rau)*(Vna(i,j)-interpl(Vn(i,j),Vn(i,j-1),da(i,j),da(i,j-1)))
  end if
 end DO
end DO
DO j=1,Jc-1
  DO i=2,Ic-1
   aE=rhok(i+1,j)*duk(i+1,j)*dy
   aW=rhok(i,j)*duk(i,j)*dy
   aN=rhoa(i,j+1)*dva(i,j+1)*dx
   aS=rhoa(i,j)*dva(i,j)*dx
   aP=aE+aW+aN+aS
   if(Proctrl=='com') then
    ww=sign(0.5,Unk(i,j))
    we=sign(0.5,Unk(i+1,j))
    ws=sign(0.5,Vna(i,j))
    wn=sign(0.5,Vna(i,j+1))
    aE=aE-Rap*(0.5-we)*Unk(i+1,j)*dy/(R*T(i+1,j)/Ma)
    aW=aW+Rap*(0.5+ww)*Unk(i,j)*dy/(R*T(i-1,j)/Ma)
    aN=aN-Rap*(0.5-wn)*Vna(i,j+1)*dx/(R*T(i,j+1)/Ma)
    if(j==1.and.(i>Ib2.or.i<Ib1)) then
     aS=aS+Rap*(0.5+ws)*Vna(i,j)*dx/(R*T(Ic+1-i,j)/Ma)
    else
     aS=aS+Rap*(0.5+ws)*Vna(i,j)*dx/(R*T(i,j-1)/Ma)
    end if
    aP=aP+Rap*((0.5+we)*Unk(i+1,j)*dy-(0.5-ww)*Unk(i,j)*dy+(0.5+wn)*Vna(i,j+1)*dx-(0.5-ws)*Vna(i,j)*dx)/(R*T(i,j)/Ma)
   end if
   aM(1,i,j)=aP
   aM(2,i,j)=aW
   aM(3,i,j)=aE
   aM(4,i,j)=aS
   aM(5,i,j)=aN
  end DO
end DO
b=0
DO j=1,Jc-1
  DO i=2,Ic-1
   b(i,j)=rhok(i,j)*Unk(i,j)*dy-rhok(i+1,j)*Unk(i+1,j)*dy+rhoa(i,j)*Vna(i,j)*dx-rhoa(i,j+1)*Vna(i,j+1)*dx
  end DO
end DO
end Subroutine PCEcoe
