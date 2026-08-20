Subroutine pUpcoe
use Aero2DCOM
implicit none
integer i,j
real(8) Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,Sf,Pxw,Pxe,Pxs,Pxn,Pyw,Pye,Pys,Pyn,Pw,Pe,Ps,Pn,dkcw,dkce,dacs,dacn
real(8) ww,we,ws,wn,noc,dncx,dncy,Tfx,Tfy,Pno
real(8) du(Ic,Jc),dv(Ic,Jc),rms(Ic,Jc)
logical(1) isCom

isCom=ProctrlFlag==COM
noc=1.d+0

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
    duk(i,j)=interpl(du(Ic,j),du(i,j),dkw(i,j))
   else if(Is>1.and.i==2) then
    duk(i,j)=interpl(0.0,du(i,j),dkw(i,j))
   else if(Is>1.and.i==Ic) then
    duk(i,j)=interpl(du(i-1,j),0.0,dkw(i,j))
   else if(i==Ip) then
    duk(i,j)=interpl(du(i-1,j),du(1,j),dkw(i,j))
   else
    duk(i,j)=interpl(du(i-1,j),du(i,j),dkw(i,j))
   end if
  end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i)
DO j=1,Jc
  DO i=Is,Ie
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    dva(i,j)=interpl(dv(Ic+1-i,j),dv(i,j),daw(i,j))
   else if(j==1) then
    dva(i,j)=0
   else if(j==Jc) then
    dva(i,j)=interpl(dv(i,j-1),0.0,daw(i,j))
   else
    dva(i,j)=interpl(dv(i,j-1),dv(i,j),daw(i,j))
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
!$OMP DO PRIVATE(i,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,we,ww,ws,wn,dkcw,dkce,dacs,dacn,&
!$OMP Pxw,Pxe,Pxs,Pxn,Pyw,Pye,Pys,Pyn,Pw,Pe,Ps,Pn,dncx,dncy,Tfx,Tfy,Pno)
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
    dkcw=dkw(i,j)
    dkce=dkw(i+1,j)
    dacs=daw(i,j)
    dacn=daw(i,j+1)
    apuM(3,i,j)=-(1-dkce)*rhok(i+1,j)*Ygae
    apuM(2,i,j)=dkcw*rhok(i,j)*Ygaw
    apuM(5,i,j)=(1-dacn)*rhoa(i,j+1)*Ygkn
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     apuM(4,i,j)=0
     apuM(1,i,j)=-(-dkce*rhok(i+1,j)*Ygae+(1-dkcw)*rhok(i,j)*Ygaw+dacn*rhoa(i,j+1)*Ygkn)
    else
     apuM(4,i,j)=-dacs*rhoa(i,j)*Ygks
     apuM(1,i,j)=-(-dkce*rhok(i+1,j)*Ygae+(1-dkcw)*rhok(i,j)*Ygaw+dacn*rhoa(i,j+1)*Ygkn-(1-dacs)*rhoa(i,j)*Ygks)
    end if
    apvM(3,i,j)=(1-dkce)*rhok(i+1,j)*Xgae
    apvM(2,i,j)=-dkcw*rhok(i,j)*Xgaw
    apvM(5,i,j)=-(1-dacn)*rhoa(i,j+1)*Xgkn
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     apvM(4,i,j)=0
     apvM(1,i,j)=-(dkce*rhok(i+1,j)*Xgae-(1-dkcw)*rhok(i,j)*Xgaw-dacn*rhoa(i,j+1)*Xgkn)
    else
     apvM(4,i,j)=dacs*rhoa(i,j)*Xgks
     apvM(1,i,j)=-(dkce*rhok(i+1,j)*Xgae-(1-dkcw)*rhok(i,j)*Xgaw-dacn*rhoa(i,j+1)*Xgkn+(1-dacs)*rhoa(i,j)*Xgks)
    end if
    if(i==1) then
     Pxw=interpl(Vol(Ic,j)*Px(Ic,j)/auM(1,Ic,j),Vol(i,j)*Px(i,j)/auM(1,i,j),dkcw)
     Pyw=interpl(Vol(Ic,j)*Py(Ic,j)/auM(1,Ic,j),Vol(i,j)*Py(i,j)/auM(1,i,j),dkcw)
    else
     Pxw=interpl(Vol(i-1,j)*Px(i-1,j)/auM(1,i-1,j),Vol(i,j)*Px(i,j)/auM(1,i,j),dkcw)
     Pyw=interpl(Vol(i-1,j)*Py(i-1,j)/auM(1,i-1,j),Vol(i,j)*Py(i,j)/auM(1,i,j),dkcw)
    end if
    if(i==Ic) then
     Pxe=interpl(Vol(i,j)*Px(i,j)/auM(1,i,j),Vol(1,j)*Px(1,j)/auM(1,1,j),dkce)
     Pye=interpl(Vol(i,j)*Py(i,j)/auM(1,i,j),Vol(1,j)*Py(1,j)/auM(1,1,j),dkce)
    else
     Pxe=interpl(Vol(i,j)*Px(i,j)/auM(1,i,j),Vol(i+1,j)*Px(i+1,j)/auM(1,i+1,j),dkce)
     Pye=interpl(Vol(i,j)*Py(i,j)/auM(1,i,j),Vol(i+1,j)*Py(i+1,j)/auM(1,i+1,j),dkce)
    end if
    if(j==1.and.(i>Ib2.or.i<Ib1)) then
     Pxs=interpl(Vol(Ic+1-i,j)*Px(Ic+1-i,j)/auM(1,Ic+1-i,j),Vol(i,j)*Px(i,j)/auM(1,i,j),dacs)
     Pys=interpl(Vol(Ic+1-i,j)*Py(Ic+1-i,j)/auM(1,Ic+1-i,j),Vol(i,j)*Py(i,j)/auM(1,i,j),dacs)
    else if(j==1) then
     Pxs=0
     Pys=0
    else
     Pxs=interpl(Vol(i,j-1)*Px(i,j-1)/auM(1,i,j-1),Vol(i,j)*Px(i,j)/auM(1,i,j),dacs)
     Pys=interpl(Vol(i,j-1)*Py(i,j-1)/auM(1,i,j-1),Vol(i,j)*Py(i,j)/auM(1,i,j),dacs)
    end if
    Pxn=interpl(Vol(i,j)*Px(i,j)/auM(1,i,j),Vol(i,j+1)*Px(i,j+1)/auM(1,i,j+1),dacn)
    Pyn=interpl(Vol(i,j)*Py(i,j)/auM(1,i,j),Vol(i,j+1)*Py(i,j+1)/auM(1,i,j+1),dacn)
    b(i,j)=Rau*(rhok(i,j)*Ygaw*Pxw-rhok(i,j)*Xgaw*Pyw-rhok(i+1,j)*Ygae*Pxe+rhok(i+1,j)*Xgae*Pye+rhoa(i,j)*Xgks*Pys-&
    rhoa(i,j)*Ygks*Pxs-rhoa(i,j+1)*Xgkn*Pyn+rhoa(i,j+1)*Ygkn*Pxn)
    Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
    if(i==1) then
     Pxw=interpl(Px(Ic,j),Px(i,j),dkcw)
     Pyw=interpl(Py(Ic,j),Py(i,j),dkcw)
     dncx=Xc(i,j)-Xc(Ic,j)
     dncy=Yc(i,j)-Yc(Ic,j)
     Tfx=-(Xfk(i,j)-Sf*dncx/dkd(i,j))
     Tfy=-(Yfk(i,j)-Sf*dncy/dkd(i,j))
     Pno=duk(i,j)*(Pxw*Tfx+Pyw*Tfy)
    else
     Pxw=interpl(Px(i-1,j),Px(i,j),dkcw)
     Pyw=interpl(Py(i-1,j),Py(i,j),dkcw)
     dncx=Xc(i,j)-Xc(i-1,j)
     dncy=Yc(i,j)-Yc(i-1,j)
     Tfx=-(Xfk(i,j)-Sf*dncx/dkd(i,j))
     Tfy=-(Yfk(i,j)-Sf*dncy/dkd(i,j))
     Pno=duk(i,j)*(Pxw*Tfx+Pyw*Tfy)
    end if
    b(i,j)=b(i,j)+rhok(i,j)*noc*Pno
    Sf=sqrt(Xfk(i+1,j)**2+Yfk(i+1,j)**2)
    if(i==Ic) then
     Pxe=interpl(Px(i,j),Px(1,j),dkcw)
     Pye=interpl(Py(i,j),Py(1,j),dkcw)
     dncx=Xc(1,j)-Xc(i,j)
     dncy=Yc(1,j)-Yc(i,j)
     Tfx=-(Xfk(i+1,j)-Sf*dncx/dkd(i+1,j))
     Tfy=-(Yfk(i+1,j)-Sf*dncy/dkd(i+1,j))
     Pno=duk(i+1,j)*(Pxe*Tfx+Pye*Tfy)
    else
     Pxe=interpl(Px(i,j),Px(i+1,j),dkcw)
     Pye=interpl(Py(i,j),Py(i+1,j),dkcw)
     dncx=Xc(i+1,j)-Xc(i,j)
     dncy=Yc(i+1,j)-Yc(i,j)
     Tfx=-(Xfk(i+1,j)-Sf*dncx/dkd(i+1,j))
     Tfy=-(Yfk(i+1,j)-Sf*dncy/dkd(i+1,j))
     Pno=duk(i+1,j)*(Pxe*Tfx+Pye*Tfy)
    end if
    b(i,j)=b(i,j)-rhok(i+1,j)*noc*Pno
    Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
    if(j==1.and.(i>Ib2.or.i<Ib1)) then
     Pxs=interpl(Px(Ic+1-i,j),Px(i,j),dacs)
     Pys=interpl(Py(Ic+1-i,j),Py(i,j),dacs)
     dncx=Xc(i,j)-Xc(Ic+1-i,j)
     dncy=Yc(i,j)-Yc(Ic+1-i,j)
     Tfx=-(Xfa(i,j)-Sf*dncx/dad(i,j))
     Tfy=-(Yfa(i,j)-Sf*dncy/dad(i,j))
     Pno=dva(i,j)*(Pxs*Tfx+Pys*Tfy)
    else if(j==1) then
     Pno=0
    else
     Pxs=interpl(Px(i,j-1),Px(i,j),dacs)
     Pys=interpl(Py(i,j-1),Py(i,j),dacs)
     dncx=Xc(i,j)-Xc(i,j-1)
     dncy=Yc(i,j)-Yc(i,j-1)
     Tfx=-(Xfa(i,j)-Sf*dncx/dad(i,j))
     Tfy=-(Yfa(i,j)-Sf*dncy/dad(i,j))
     Pno=dva(i,j)*(Pxs*Tfx+Pys*Tfy)
    end if
    b(i,j)=b(i,j)+rhoa(i,j)*noc*Pno
    Sf=sqrt(Xfa(i,j+1)**2+Yfa(i,j+1)**2)
    Pxn=interpl(Px(i,j),Px(i,j+1),dacn)
    Pyn=interpl(Py(i,j),Py(i,j+1),dacn)
    dncx=Xc(i,j+1)-Xc(i,j)
    dncy=Yc(i,j+1)-Yc(i,j)
    Tfx=-(Xfa(i,j+1)-Sf*dncx/dad(i,j+1))
    Tfy=-(Yfa(i,j+1)-Sf*dncy/dad(i,j+1))
    Pno=dva(i,j+1)*(Pxn*Tfx+Pyn*Tfy)
    b(i,j)=b(i,j)-rhoa(i,j+1)*noc*Pno
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
