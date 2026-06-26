Subroutine pUpcoe
use Aero2DCOM
implicit none
integer i,j
real(8) Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,Sf,Pxw,Pxe,Pxs,Pxn,Pyw,Pye,Pys,Pyn,Pw,Pe,Ps,Pn
real(8) ww,we,ws,wn
real(8),external::interpl
real(8) du(Ic,Jc),dv(Ic,Jc),rms(Ic,Jc)
logical(1) isCom

isCom=ProctrlFlag==COM

!$OMP PARALLEL
!$OMP DO PRIVATE(i)
DO j=1,Jc-1
  DO i=Is,Ie
   du(i,j)=Rau*Vol(i,j)/auM(1,i,j)
   dv(i,j)=Rau*Vol(i,j)/auM(1,i,j)
  end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i)
DO j=1,Jc-1
  DO i=Is,Ie+1
   if(i==1) then
    duk(i,j)=interpl(du(Ic,j),du(i,j),5d-1)
   else if(Is>1.and.i==2) then
    duk(i,j)=interpl(0.0,du(i,j),5d-1)
   else if(Is>1.and.i==Ic) then
    duk(i,j)=interpl(du(i-1,j),0.0,5d-1)
   else if(i==Ip) then
    duk(i,j)=interpl(du(i-1,j),du(1,j),5d-1) 
   else
    duk(i,j)=interpl(du(i-1,j),du(i,j),5d-1)
   end if
  end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i)
DO j=1,Jc
  DO i=Is,Ie
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    dva(i,j)=interpl(dv(i,j),dv(Ic+1-i,j),5d-1)
   else if(j==1) then
    dva(i,j)=0
   else if(j==Jc) then
    dva(i,j)=interpl(dv(i,j-1),0.0,5d-1)
   else
    dva(i,j)=interpl(dv(i,j-1),dv(i,j),5d-1)
   end if
  end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
aM=0
aM(1,:,:)=1
b=P
apuM=0
apvM=0
!$OMP END WORKSHARE
!$OMP DO PRIVATE(i,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,we,ww,ws,wn,Pxw,Pxe,Pxs,Pxn,Pyw,Pye,Pys,Pyn,Pw,Pe,Ps,Pn)
DO j=1,Jc-1
  DO i=Is,Ie
    Xgaw=-Yfk(i,j)
    Ygaw=Xfk(i,j)
    Xgae=-Yfk(i+1,j)
    Ygae=Xfk(i+1,j)
    Xgks=Yfa(i,j)
    Ygks=-Xfa(i,j)
    Xgkn=Yfa(i,j+1)
    Ygkn=-Xfa(i,j+1)
    Sf=sqrt(Xfk(i+1,j)**2+Yfk(i+1,j)**2)
    aM(3,i,j)=rhok(i+1,j)*duk(i+1,j)*Sf/dkd(i+1,j)
    Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
    aM(2,i,j)=rhok(i,j)*duk(i,j)*Sf/dkd(i,j)
    Sf=sqrt(Xfa(i,j+1)**2+Yfa(i,j+1)**2)
    aM(5,i,j)=rhoa(i,j+1)*dva(i,j+1)*Sf/dad(i,j+1)
    Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
    aM(4,i,j)=rhoa(i,j)*dva(i,j)*Sf/dad(i,j)
    aM(1,i,j)=aM(3,i,j)+aM(2,i,j)+aM(5,i,j)+aM(4,i,j)
    if(isCom) then
      we=sign(0.5d0,Unk(i+1,j))
      if(i==Ic) then
        aM(3,i,j)=aM(3,i,j)-(0.5d0-we)*Unk(i+1,j)/(R*T(1,j)/Ma)    
      else
        aM(3,i,j)=aM(3,i,j)-(0.5d0-we)*Unk(i+1,j)/(R*T(i+1,j)/Ma)
      end if
      ww=sign(0.5d0,Unk(i,j))
      if(i==1) then
        aM(2,i,j)=aM(2,i,j)+(0.5d0+ww)*Unk(i,j)/(R*T(Ic,j)/Ma)
      else
        aM(2,i,j)=aM(2,i,j)+(0.5d0+ww)*Unk(i,j)/(R*T(i-1,j)/Ma)
      end if
      wn=sign(0.5d0,Vna(i,j+1))
      aM(5,i,j)=aM(5,i,j)-(0.5d0-wn)*Vna(i,j+1)/(R*T(i,j+1)/Ma)
      ws=sign(0.5d0,Vna(i,j))
      if(j==1.and.(i<Ib1.or.i>Ib2)) then
        aM(4,i,j)=aM(4,i,j)+(0.5d0+ws)*Vna(i,j)/(R*T(Ic+1-i,j)/Ma)
      else if(j>1) then
        aM(4,i,j)=aM(4,i,j)+(0.5d0+ws)*Vna(i,j)/(R*T(i,j-1)/Ma)
      end if
      aM(1,i,j)=aM(1,i,j)+((0.5d0+we)*Unk(i+1,j)-(0.5d0-ww)*Unk(i,j)+&
      (0.5d0+wn)*Vna(i,j+1)-(0.5d0-ws)*Vna(i,j))/(R*T(i,j)/Ma)
    end if
    apuM(3,i,j)=-0.5*rhok(i+1,j)*Ygae
    apuM(2,i,j)=0.5*rhok(i,j)*Ygaw
    apuM(5,i,j)=0.5*rhoa(i,j+1)*Ygkn
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     apuM(4,i,j)=0
    else
     apuM(4,i,j)=-0.5*rhoa(i,j)*Ygks
    end if
    apuM(1,i,j)=-(apuM(3,i,j)+apuM(2,i,j)+apuM(5,i,j)+apuM(4,i,j))
    apvM(3,i,j)=0.5*rhok(i+1,j)*Xgae
    apvM(2,i,j)=-0.5*rhok(i,j)*Xgaw
    apvM(5,i,j)=-0.5*rhoa(i,j+1)*Xgkn
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     apvM(4,i,j)=0
    else
     apvM(4,i,j)=0.5*rhoa(i,j)*Xgks
    end if
    apvM(1,i,j)=-(apvM(3,i,j)+apvM(2,i,j)+apvM(5,i,j)+apvM(4,i,j))
    if(i==1) then
     Pxw=0.5*(Vol(i,j)*Px(i,j)/auM(1,i,j)+Vol(Ic,j)*Px(Ic,j)/auM(1,Ic,j))
     Pyw=0.5*(Vol(i,j)*Py(i,j)/auM(1,i,j)+Vol(Ic,j)*Py(Ic,j)/auM(1,Ic,j))
    else
     Pxw=0.5*(Vol(i,j)*Px(i,j)/auM(1,i,j)+Vol(i-1,j)*Px(i-1,j)/auM(1,i-1,j))
     Pyw=0.5*(Vol(i,j)*Py(i,j)/auM(1,i,j)+Vol(i-1,j)*Py(i-1,j)/auM(1,i-1,j))
    end if
    if(i==Ic) then
     Pxe=0.5*(Vol(i,j)*Px(i,j)/auM(1,i,j)+Vol(1,j)*Px(1,j)/auM(1,1,j))
     Pye=0.5*(Vol(i,j)*Py(i,j)/auM(1,i,j)+Vol(1,j)*Py(1,j)/auM(1,1,j))
    else
     Pxe=0.5*(Vol(i,j)*Px(i,j)/auM(1,i,j)+Vol(i+1,j)*Px(i+1,j)/auM(1,i+1,j))
     Pye=0.5*(Vol(i,j)*Py(i,j)/auM(1,i,j)+Vol(i+1,j)*Py(i+1,j)/auM(1,i+1,j))
    end if
    if(j==1.and.(i>Ib2.or.i<Ib1)) then
     Pxs=0.5*(Vol(i,j)*Px(i,j)/auM(1,i,j)+Vol(Ic+1-i,j)*Px(Ic+1-i,j)/auM(1,Ic+1-i,j))
     Pys=0.5*(Vol(i,j)*Py(i,j)/auM(1,i,j)+Vol(Ic+1-i,j)*Py(Ic+1-i,j)/auM(1,Ic+1-i,j))
    else if(j==1) then
     Pxs=0
     Pys=0
    else
     Pxs=0.5*(Vol(i,j)*Px(i,j)/auM(1,i,j)+Vol(i,j-1)*Px(i,j-1)/auM(1,i,j-1))
     Pys=0.5*(Vol(i,j)*Py(i,j)/auM(1,i,j)+Vol(i,j-1)*Py(i,j-1)/auM(1,i,j-1))
    end if
    Pxn=0.5*(Vol(i,j)*Px(i,j)/auM(1,i,j)+Vol(i,j+1)*Px(i,j+1)/auM(1,i,j+1))
    Pyn=0.5*(Vol(i,j)*Py(i,j)/auM(1,i,j)+Vol(i,j+1)*Py(i,j+1)/auM(1,i,j+1))
    b(i,j)=Rau*(rhok(i,j)*Ygaw*Pxw-rhok(i,j)*Xgaw*Pyw-rhok(i+1,j)*Ygae*Pxe+rhok(i+1,j)*Xgae*Pye+rhoa(i,j)*Xgks*Pys-&
    rhoa(i,j)*Ygks*Pxs-rhoa(i,j+1)*Xgkn*Pyn+rhoa(i,j+1)*Ygkn*Pxn)
    if(isCom) then
      if(i==Ic) then
        Pe=(0.5d0+we)*P(i,j)/(R*T(i,j)/Ma)+(0.5d0-we)*P(1,j)/(R*T(1,j)/Ma)
      else
        Pe=(0.5d0+we)*P(i,j)/(R*T(i,j)/Ma)+(0.5d0-we)*P(i+1,j)/(R*T(i+1,j)/Ma)
      end if
      if(i==1) then
        Pw=(0.5d0-ww)*P(i,j)/(R*T(i,j)/Ma)+(0.5d0+ww)*P(Ic,j)/(R*T(Ic,j)/Ma)
      else
        Pw=(0.5d0-ww)*P(i,j)/(R*T(i,j)/Ma)+(0.5d0+ww)*P(i-1,j)/(R*T(i-1,j)/Ma)
      end if
      Pn=(0.5d0+wn)*P(i,j)/(R*T(i,j)/Ma)+(0.5d0-wn)*P(i,j+1)/(R*T(i,j+1)/Ma)
      if(j==1.and.(i<Ib1.or.i>Ib2)) then
        Ps=(0.5d0-ws)*P(i,j)/(R*T(i,j)/Ma)+(0.5d0+ws)*P(Ic+1-i,j)/(R*T(Ic+1-i,j)/Ma)
      else if(j==1) then
        Ps=P(i,j)/(R*T(i,j)/Ma)
      else
        Ps=(0.5d0-ws)*P(i,j)/(R*T(i,j)/Ma)+(0.5d0+ws)*P(i,j-1)/(R*T(i,j-1)/Ma)
      end if
      b(i,j)=b(i,j)+Pe*Unk(i+1,j)-Pw*Unk(i,j)+Pn*Vna(i,j+1)-Ps*Vna(i,j)
    end if
  end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
rms=0
!$OMP END WORKSHARE
!$OMP DO PRIVATE(i)
DO j=1,Jc-1
  DO i=Is,Ie
    rms(i,j)=aM(5,i,j)*P(i,j+1)+b(i,j)-aM(1,i,j)*P(i,j)+apuM(5,i,j)*U(i,j+1)-apuM(1,i,j)*U(i,j)+apvM(5,i,j)*V(i,j+1)-&
    apvM(1,i,j)*V(i,j)
    if(j==1.and.(i>Ib2.or.i<Ib1)) then
      rms(i,j)=rms(i,j)+aM(4,i,j)*P(Ic+1-i,j)+apuM(4,i,j)*U(Ic+1-i,j)+apvM(4,i,j)*V(Ic+1-i,j)
    else if(j==1) then
      rms(i,j)=rms(i,j)+aM(4,i,j)*P(i,j)+apuM(4,i,j)*U(i,j)+apvM(4,i,j)*V(i,j)
    else
      rms(i,j)=rms(i,j)+aM(4,i,j)*P(i,j-1)+apuM(4,i,j)*U(i,j-1)+apvM(4,i,j)*V(i,j-1)
    end if
    if(i==1) then
      rms(i,j)=rms(i,j)+aM(2,i,j)*P(Ic,j)+apuM(2,i,j)*U(Ic,j)+apvM(2,i,j)*V(Ic,j)
    else
      rms(i,j)=rms(i,j)+aM(2,i,j)*P(i-1,j)+apuM(2,i,j)*U(i-1,j)+apvM(2,i,j)*V(i-1,j)
    end if
    if(i==Ic) then
      rms(i,j)=rms(i,j)+aM(3,i,j)*P(1,j)+apuM(3,i,j)*U(1,j)+apvM(3,i,j)*V(1,j)
    else
      rms(i,j)=rms(i,j)+aM(3,i,j)*P(i+1,j)+apuM(3,i,j)*U(i+1,j)+apvM(3,i,j)*V(i+1,j)
    end if
  end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
rmsm=sum(abs(rms))/(Ic*Jc)
!$OMP END WORKSHARE
!$OMP END PARALLEL
end Subroutine pUpcoe
