Subroutine PCEcoe
use Aero2DCOM
implicit none
integer i,j
real(8) Up,Vp,Unpw,Unpe,Vnps,Vnpn
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
     DO i=2,Ic-1
      if(i==2) then
          wdu(i,j)=0.5*du(i,j)
          Unpw=0.5*(Unp(i,j)+Un(i-1,j))
      else
          wdu(i,j)=0.5*(du(i,j)+du(i-1,j))
          Unpw=0.5*(Unp(i,j)+Unp(i-1,j))
      end if
      if(i==Ic-1) then
          edu(i,j)=0.5*du(i,j)
          Unpe=0.5*(Unp(i,j)+Un(i+1,j))
      else
          edu(i,j)=0.5*(du(i,j)+du(i+1,j))
          Unpe=0.5*(Unp(i,j)+Unp(i+1,j))
      end if
      if(j==1.and.(i>Ib2.or.i<Ib1)) then
          sdv(i,j)=0.5*(dv(i,j)+dv(Ic+1-i,j))
          Vnps=0.5*(Vnp(i,j)-Vnp(Ic+1-i,j))
      else if(j==1) then
          sdv(i,j)=dv(i,j)
          Vnps=Vnp(i,j)
      else
          sdv(i,j)=0.5*(dv(i,j)+dv(i,j-1))
          Vnps=0.5*(Vnp(i,j)+Vnp(i,j-1))
      end if
      if(j==Jc-1) then
          ndv(i,j)=0.5*dv(i,j)
          Vnpn=0.5*(Vnp(i,j)+Vn(i,j+1))
      else
          ndv(i,j)=0.5*(dv(i,j)+dv(i,j+1))
          Vnpn=0.5*(Vnp(i,j)+Vnp(i,j+1))
      end if
      Unw(i,j)=Unpw+wdu(i,j)*(P(i-1,j)-P(i,j))+(1-Rau)*(Unw(i,j)-0.5*(Un(i,j)+Un(i-1,j)))
      Une(i,j)=Unpe+edu(i,j)*(P(i,j)-P(i+1,j))+(1-Rau)*(Une(i,j)-0.5*(Un(i,j)+Un(i+1,j)))
      Vnn(i,j)=Vnpn+ndv(i,j)*(P(i,j)-P(i,j+1))+(1-Rau)*(Vnn(i,j)-0.5*(Vn(i,j)+Vn(i,j+1)))
      if(j==1.and.(i>Ib2.or.i<Ib1)) then
      Vns(i,j)=Vnps+sdv(i,j)*(P(Ic+1-i,j)-P(i,j))+(1-Rau)*(Vns(i,j)-0.5*(Vn(i,j)-Vn(Ic+1-i,j)))
      else if(j==1) then
      Vns(i,j)=0
      else
      Vns(i,j)=Vnps+sdv(i,j)*(P(i,j-1)-P(i,j))+(1-Rau)*(Vns(i,j)-0.5*(Vn(i,j)+Vn(i,j-1)))
      end if
     end DO
  end DO
  DO j=1,Jc-1
    DO i=2,Ic-1
      if(Proctrl=='incom') then
      aE(i,j)=roue(i,j)*edu(i,j)*dy
      aW(i,j)=rouw(i,j)*wdu(i,j)*dy
      aN(i,j)=roun(i,j)*ndv(i,j)*dx
      if(j==1.and.i>=Ib1.and.i<=Ib2) then
      aS(i,j)=0
      else
      aS(i,j)=rous(i,j)*sdv(i,j)*dx
      end if
      aP(i,j)=aE(i,j)+aW(i,j)+aN(i,j)+aS(i,j)
      else if(Proctrl=='com') then
      aE(i,j)=roue(i,j)*edu(i,j)*dy-Rap*(0.5-we(i,j))*Une(i,j)*dy/(R*T(i+1,j)/Ma)
      aW(i,j)=rouw(i,j)*wdu(i,j)*dy+Rap*(0.5+ww(i,j))*Unw(i,j)*dy/(R*T(i-1,j)/Ma)
      aN(i,j)=roun(i,j)*ndv(i,j)*dx-Rap*(0.5-wn(i,j))*Vnn(i,j)*dx/(R*T(i,j+1)/Ma)
      if(j==1.and.i>=Ib1.and.i<=Ib2) then
      aS(i,j)=0
      aP(i,j)=roue(i,j)*edu(i,j)*dy+rouw(i,j)*wdu(i,j)*dy+roun(i,j)*ndv(i,j)*dx+Rap*((0.5+we(i,j))*Une(i,j)*dy-&
      (0.5-ww(i,j))*Unw(i,j)*dy+(0.5+wn(i,j))*Vnn(i,j)*dx)/(R*T(i,j)/Ma)
      else if(j==1) then
      aS(i,j)=rous(i,j)*sdv(i,j)*dx+Rap*(0.5+ws(i,j))*Vns(i,j)*dx/(R*T(Ic+1-i,j)/Ma)
      aP(i,j)=roue(i,j)*edu(i,j)*dy+rouw(i,j)*wdu(i,j)*dy+roun(i,j)*ndv(i,j)*dx+rous(i,j)*sdv(i,j)*dx+&
      Rap*((0.5+we(i,j))*Une(i,j)*dy-(0.5-ww(i,j))*Unw(i,j)*dy+(0.5+wn(i,j))*Vnn(i,j)*dx-(0.5-ws(i,j))*Vns(i,j)*dx)/(R*T(i,j)/Ma)
      else
      aS(i,j)=rous(i,j)*sdv(i,j)*dx+Rap*(0.5+ws(i,j))*Vns(i,j)*dx/(R*T(i,j-1)/Ma)
      aP(i,j)=roue(i,j)*edu(i,j)*dy+rouw(i,j)*wdu(i,j)*dy+roun(i,j)*ndv(i,j)*dx+rous(i,j)*sdv(i,j)*dx+&
      Rap*((0.5+we(i,j))*Une(i,j)*dy-(0.5-ww(i,j))*Unw(i,j)*dy+(0.5+wn(i,j))*Vnn(i,j)*dx-(0.5-ws(i,j))*Vns(i,j)*dx)/(R*T(i,j)/Ma)
      end if
      end if
    end DO
  end DO
  b=0
  DO j=1,Jc-1
    DO i=2,Ic-1
      b(i,j)=rouw(i,j)*Unw(i,j)*dy-roue(i,j)*Une(i,j)*dy+rous(i,j)*Vns(i,j)*dx-roun(i,j)*Vnn(i,j)*dx
    end DO
  end DO
end Subroutine PCEcoe
