Subroutine pUpcoe
use Aero2DCOM
implicit none
integer i,j
real(8) Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,Pxw,Pxe,Pxs,Pxn,Pyw,Pye,Pys,Pyn,Pw,Pe,Ps,Pn,dkc,dkw,dke,dac,das,dan
real(8) ww,we,ws,wn,dwk,dwa,noc,Pak,Pka,Psw,Psc,Pwn,Pnc,Pse,Pec,Pwc
real(8) du(Ic,Jc),dv(Ic,Jc),rms(Ic,Jc),Pkno(Ip,Jc),Pano(Ic,Jp)
logical(1) isCom

isCom=ProctrlFlag==COM
noc=1.d+0

!$OMP PARALLEL
!$OMP DO PRIVATE(i)
DO j=1,Jc-1
  DO i=Is,Ie
   du(i,j)=Rau*(Yga(i,j)**2+Xga(i,j)**2)*dy/auM(1,i,j)
   dv(i,j)=Rau*(Ygk(i,j)**2+Xgk(i,j)**2)*dx/auM(1,i,j)
   dwno(i,j)=Rau*b1(i,j)*Jg(i,j)/auM(1,i,j)
  end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i,dwk,Pak,Psw,Psc,Pwn,Pnc)
DO j=1,Jc-1
  DO i=Is,Ie+1
   if(i==1) then
    duk(i,j)=interpl(du(i,j),du(Ic,j),dk(i,j),dk(i,j))
    dwk=interpl(dwno(i,j)*dy,dwno(Ic,j)*dy,dk(i,j),dk(Ic,j))
   else if(Is>1.and.i==2) then
    duk(i,j)=interpl(du(i,j),0.0,dk(i,j),dk(i,j))
    dwk=interpl(dwno(i,j)*dy,0.0,dk(i,j),dk(i-1,j))
   else if(Is>1.and.i==Ic) then
    duk(i,j)=interpl(0.0,du(i-1,j),dk(i,j),dk(i,j))
    dwk=interpl(0.0,dwno(i-1,j)*dy,dk(i,j),dk(i-1,j))
   else if(i==Ip) then
    duk(i,j)=interpl(du(1,j),du(i-1,j),dk(i-1,j),dk(i-1,j))
    dwk=interpl(dwno(1,j)*dy,dwno(i-1,j)*dy,dk(1,j),dk(i-1,j))
   else
    duk(i,j)=interpl(du(i,j),du(i-1,j),dk(i,j),dk(i,j))
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
   Pkno(i,j)=noc*dwk*Pak
  end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i,dwa,Pka,Psw,Pse,Pec,Pwc)
DO j=1,Jc
  DO i=Is,Ie
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    dva(i,j)=interpl(dv(i,j),dv(Ic+1-i,j),da(i,j),da(i,j))
    dwa=interpl(dwno(i,j)*dx,dwno(Ic+1-i,j)*dx,da(i,j),da(Ic+1-i,j))
   else if(j==1) then
    dva(i,j)=0
    dwa=0
   else if(j==Jc) then
    dva(i,j)=interpl(0.0,dv(i,j-1),da(i,j),da(i,j))
    dwa=interpl(0.0,dwno(i,j-1)*dx,da(i,j),da(i,j-1))
   else
    dva(i,j)=interpl(dv(i,j),dv(i,j-1),da(i,j),da(i,j))
    dwa=interpl(dwno(i,j)*dx,dwno(i,j-1)*dx,da(i,j),da(i,j-1))
   end if
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
   Pano(i,j)=noc*dwa*Pka
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
!$OMP DO PRIVATE(i,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,dkc,dkw,dke,dac,das,dan,we,ww,ws,wn,&
!$OMP Pxw,Pxe,Pxs,Pxn,Pyw,Pye,Pys,Pyn,Pw,Pe,Ps,Pn)
DO j=1,Jc-1
  DO i=Is,Ie
    Xgaw=(Xg(i,j+1)-Xg(i,j))/dy
    Ygaw=(Yg(i,j+1)-Yg(i,j))/dy
    Xgae=(Xg(i+1,j+1)-Xg(i+1,j))/dy
    Ygae=(Yg(i+1,j+1)-Yg(i+1,j))/dy
    Xgks=(Xg(i+1,j)-Xg(i,j))/dx
    Ygks=(Yg(i+1,j)-Yg(i,j))/dx
    Xgkn=(Xg(i+1,j+1)-Xg(i,j+1))/dx
    Ygkn=(Yg(i+1,j+1)-Yg(i,j+1))/dx
    aM(3,i,j)=rhok(i+1,j)*duk(i+1,j)*dy
    aM(2,i,j)=rhok(i,j)*duk(i,j)*dy
    aM(5,i,j)=rhoa(i,j+1)*dva(i,j+1)*dx
    aM(4,i,j)=rhoa(i,j)*dva(i,j)*dx
    aM(1,i,j)=aM(3,i,j)+aM(2,i,j)+aM(5,i,j)+aM(4,i,j)
    if(isCom) then
      we=sign(0.5d0,Unk(i+1,j))
      if(i==Ic) then
        aM(3,i,j)=aM(3,i,j)-(0.5d0-we)*Unk(i+1,j)*dy/(R*T(1,j)/Ma)    
      else
        aM(3,i,j)=aM(3,i,j)-(0.5d0-we)*Unk(i+1,j)*dy/(R*T(i+1,j)/Ma)
      end if
      ww=sign(0.5d0,Unk(i,j))
      if(i==1) then
        aM(2,i,j)=aM(2,i,j)+(0.5d0+ww)*Unk(i,j)*dy/(R*T(Ic,j)/Ma)
      else
        aM(2,i,j)=aM(2,i,j)+(0.5d0+ww)*Unk(i,j)*dy/(R*T(i-1,j)/Ma)
      end if
      wn=sign(0.5d0,Vna(i,j+1))
      aM(5,i,j)=aM(5,i,j)-(0.5d0-wn)*Vna(i,j+1)*dx/(R*T(i,j+1)/Ma)
      ws=sign(0.5d0,Vna(i,j))
      if(j==1.and.(i<Ib1.or.i>Ib2)) then
        aM(4,i,j)=aM(4,i,j)+(0.5d0+ws)*Vna(i,j)*dx/(R*T(Ic+1-i,j)/Ma)
      else if(j>1) then
        aM(4,i,j)=aM(4,i,j)+(0.5d0+ws)*Vna(i,j)*dx/(R*T(i,j-1)/Ma)
      end if
      aM(1,i,j)=aM(1,i,j)+((0.5d0+we)*Unk(i+1,j)*dy-(0.5d0-ww)*Unk(i,j)*dy+&
      (0.5d0+wn)*Vna(i,j+1)*dx-(0.5d0-ws)*Vna(i,j)*dx)/(R*T(i,j)/Ma)
    end if
    dkc=dk(i,j)
    if(i==1) then
     dkw=dk(Ic,j)
    else
     dkw=dk(i-1,j)
    end if
    if(i==Ic) then
     dke=dk(1,j)
    else
     dke=dk(i+1,j)
    end if
    dac=da(i,j)
    dan=da(i,j+1)
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     das=da(i,j)
    else if(j==1) then
     das=da(Ic+1-i,j)
    else
     das=da(i,j-1)
    end if
    apuM(3,i,j)=-dkc*rhok(i+1,j)*Ygae*dy/(dkc+dke)
    apuM(2,i,j)=dkc*rhok(i,j)*Ygaw*dy/(dkc+dkw)
    apuM(5,i,j)=dac*rhoa(i,j+1)*Ygkn*dx/(dac+dan)
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     apuM(4,i,j)=0
     apuM(1,i,j)=-(-dke*rhok(i+1,j)*Ygae*dy/(dkc+dke)+dkw*rhok(i,j)*Ygaw*dy/(dkc+dkw)+dan*rhoa(i,j+1)*Ygkn*dx/(dac+dan))
    else
     apuM(4,i,j)=-dac*rhoa(i,j)*Ygks*dx/(dac+das)
     apuM(1,i,j)=-(-dke*rhok(i+1,j)*Ygae*dy/(dkc+dke)+dkw*rhok(i,j)*Ygaw*dy/(dkc+dkw)+dan*rhoa(i,j+1)*Ygkn*dx/(dac+dan)-das*rhoa(i,j)*Ygks*dx/(dac+das))
    end if
    apvM(3,i,j)=dkc*rhok(i+1,j)*Xgae*dy/(dkc+dke)
    apvM(2,i,j)=-dkc*rhok(i,j)*Xgaw*dy/(dkc+dkw)
    apvM(5,i,j)=-dac*rhoa(i,j+1)*Xgkn*dx/(dac+dan)
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     apvM(4,i,j)=0
     apvM(1,i,j)=-(dke*rhok(i+1,j)*Xgae*dy/(dkc+dke)-dkw*rhok(i,j)*Xgaw*dy/(dkc+dkw)-dan*rhoa(i,j+1)*Xgkn*dx/(dac+dan))
    else
     apvM(4,i,j)=dac*rhoa(i,j)*Xgks*dx/(dac+das)
     apvM(1,i,j)=-(dke*rhok(i+1,j)*Xgae*dy/(dkc+dke)-dkw*rhok(i,j)*Xgaw*dy/(dkc+dkw)-dan*rhoa(i,j+1)*Xgkn*dx/(dac+dan)+das*rhoa(i,j)*Xgks*dx/(dac+das))
    end if
    if(i==1) then
     Pxw=interpl(Jg(i,j)*Px(i,j)/auM(1,i,j),Jg(Ic,j)*Px(Ic,j)/auM(1,Ic,j),dkc,dkw)
     Pyw=interpl(Jg(i,j)*Py(i,j)/auM(1,i,j),Jg(Ic,j)*Py(Ic,j)/auM(1,Ic,j),dkc,dkw)
    else
     Pxw=interpl(Jg(i,j)*Px(i,j)/auM(1,i,j),Jg(i-1,j)*Px(i-1,j)/auM(1,i-1,j),dkc,dkw)
     Pyw=interpl(Jg(i,j)*Py(i,j)/auM(1,i,j),Jg(i-1,j)*Py(i-1,j)/auM(1,i-1,j),dkc,dkw)
    end if
    if(i==Ic) then
     Pxe=interpl(Jg(i,j)*Px(i,j)/auM(1,i,j),Jg(1,j)*Px(1,j)/auM(1,1,j),dkc,dke)
     Pye=interpl(Jg(i,j)*Py(i,j)/auM(1,i,j),Jg(1,j)*Py(1,j)/auM(1,1,j),dkc,dke)
    else
     Pxe=interpl(Jg(i,j)*Px(i,j)/auM(1,i,j),Jg(i+1,j)*Px(i+1,j)/auM(1,i+1,j),dkc,dke)
     Pye=interpl(Jg(i,j)*Py(i,j)/auM(1,i,j),Jg(i+1,j)*Py(i+1,j)/auM(1,i+1,j),dkc,dke)
    end if
    if(j==1.and.(i>Ib2.or.i<Ib1)) then
     Pxs=interpl(Jg(i,j)*Px(i,j)/auM(1,i,j),Jg(Ic+1-i,j)*Px(Ic+1-i,j)/auM(1,Ic+1-i,j),dac,das)
     Pys=interpl(Jg(i,j)*Py(i,j)/auM(1,i,j),Jg(Ic+1-i,j)*Py(Ic+1-i,j)/auM(1,Ic+1-i,j),dac,das)
    else if(j==1) then
     Pxs=0
     Pys=0
    else
     Pxs=interpl(Jg(i,j)*Px(i,j)/auM(1,i,j),Jg(i,j-1)*Px(i,j-1)/auM(1,i,j-1),dac,das)
     Pys=interpl(Jg(i,j)*Py(i,j)/auM(1,i,j),Jg(i,j-1)*Py(i,j-1)/auM(1,i,j-1),dac,das)
    end if
    Pxn=interpl(Jg(i,j)*Px(i,j)/auM(1,i,j),Jg(i,j+1)*Px(i,j+1)/auM(1,i,j+1),dac,dan)
    Pyn=interpl(Jg(i,j)*Py(i,j)/auM(1,i,j),Jg(i,j+1)*Py(i,j+1)/auM(1,i,j+1),dac,dan)
    b(i,j)=dx*dy*Rau*(rhok(i,j)*Ygaw*Pxw-rhok(i,j)*Xgaw*Pyw-rhok(i+1,j)*Ygae*Pxe+rhok(i+1,j)*Xgae*Pye+rhoa(i,j)*Xgks*Pys-&
    rhoa(i,j)*Ygks*Pxs-rhoa(i,j+1)*Xgkn*Pyn+rhoa(i,j+1)*Ygkn*Pxn)
    b(i,j)=b(i,j)+rhok(i,j)*Pkno(i,j)*dy-rhok(i+1,j)*Pkno(i+1,j)*dy+rhoa(i,j)*Pano(i,j)*dx-rhoa(i,j+1)*Pano(i,j+1)*dx
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
      b(i,j)=b(i,j)+Pe*Unk(i+1,j)*dy-Pw*Unk(i,j)*dy+Pn*Vna(i,j+1)*dx-Ps*Vna(i,j)*dx
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
