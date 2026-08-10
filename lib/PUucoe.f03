Subroutine pUucoe(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,DF,Sf,dkcw,dkce,dacs,dacn,Gxf,Gyf,dncx,dncy,Tfx,Tfy
real(8) F(Ic,Jc),Ga(Ic,Jc),Fx(Ic,Jc),Fy(Ic,Jc),Fwall(Ib1:Ib2)
real(8) Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc),Dw(Ic,Jc),De(Ic,Jc),Ds(Ic,Jc),Dn(Ic,Jc),bno(Ic,Jc),cor(Ic,Jc)
integer scalar
logical(1) isU,isV,isKe,isSst,isSa,isLam,isInv,isCom,isIncom,isWf,isLr,isParvel

isU = scalar==VELX
isV = scalar==VELY
isKe = TurmodFlag==KE
isSst = TurmodFlag==SST
isSa = TurmodFlag==SA
isLam = TurmodFlag==LAM
isInv = TurmodFlag==INV
isCom = ProctrlFlag==COM
isIncom = ProctrlFlag==INCOM
isWf = WalltreatFlag==WF
isLr = WalltreatFlag==LR
isParvel = wallfunutype==PARVEL

!$OMP PARALLEL
!$OMP WORKSHARE
Fwall=0
Ga=mu+mut
!$OMP END WORKSHARE
if(isU) then
 !$OMP WORKSHARE
 F=U
 Fx=Ux
 Fy=Uy
 !$OMP END WORKSHARE
else if(isV) then
 !$OMP WORKSHARE
 F=V
 Fx=Vx
 Fy=Vy
 !$OMP END WORKSHARE
end if

!$OMP DO PRIVATE(i,Sf,Gxf,Gyf,dncx,dncy,Tfx,Tfy)
DO j=1,Jc-1
 DO i=Is,Ie
  Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
  if(i==1) then
   Dw(i,j)=Sf*interpl(Ga(Ic,j),Ga(i,j),dkw(i,j))/dkd(i,j)
   Gxf=interpl(Ga(Ic,j)*Fx(Ic,j),Ga(i,j)*Fx(i,j),dkw(i,j))
   Gyf=interpl(Ga(Ic,j)*Fy(Ic,j),Ga(i,j)*Fy(i,j),dkw(i,j))
   dncx=Xc(Ic,j)-Xc(i,j)
   dncy=Yc(Ic,j)-Yc(i,j)
  else
   Dw(i,j)=Sf*interpl(Ga(i-1,j),Ga(i,j),dkw(i,j))/dkd(i,j)
   Gxf=interpl(Ga(i-1,j)*Fx(i-1,j),Ga(i,j)*Fx(i,j),dkw(i,j))
   Gyf=interpl(Ga(i-1,j)*Fy(i-1,j),Ga(i,j)*Fy(i,j),dkw(i,j))
   dncx=Xc(i-1,j)-Xc(i,j)
   dncy=Yc(i-1,j)-Yc(i,j)
  end if
  Tfx=-Xfk(i,j)-Sf**2*dncx/(-dncx*Xfk(i,j)-dncy*Yfk(i,j))
  Tfy=-Yfk(i,j)-Sf**2*dncy/(-dncx*Xfk(i,j)-dncy*Yfk(i,j))
  bno(i,j)=0.0
  bno(i,j)=bno(i,j)+Gxf*Tfx+Gyf*Tfy
  Sf=sqrt(Xfk(i+1,j)**2+Yfk(i+1,j)**2)
  if(i==Ic) then
   De(i,j)=Sf*interpl(Ga(i,j),Ga(1,j),dkw(i+1,j))/dkd(i+1,j)
   Gxf=interpl(Ga(i,j)*Fx(i,j),Ga(1,j)*Fx(1,j),dkw(i+1,j))
   Gyf=interpl(Ga(i,j)*Fy(i,j),Ga(1,j)*Fy(1,j),dkw(i+1,j))
   dncx=Xc(1,j)-Xc(i,j)
   dncy=Yc(1,j)-Yc(i,j)
  else
   De(i,j)=Sf*interpl(Ga(i,j),Ga(i+1,j),dkw(i+1,j))/dkd(i+1,j)
   Gxf=interpl(Ga(i,j)*Fx(i,j),Ga(i+1,j)*Fx(i+1,j),dkw(i+1,j))
   Gyf=interpl(Ga(i,j)*Fy(i,j),Ga(i+1,j)*Fy(i+1,j),dkw(i+1,j))
   dncx=Xc(i+1,j)-Xc(i,j)
   dncy=Yc(i+1,j)-Yc(i,j)
  end if
  Tfx=Xfk(i+1,j)-Sf**2*dncx/(dncx*Xfk(i+1,j)+dncy*Yfk(i+1,j))
  Tfy=Yfk(i+1,j)-Sf**2*dncy/(dncx*Xfk(i+1,j)+dncy*Yfk(i+1,j))
  bno(i,j)=bno(i,j)+Gxf*Tfx+Gyf*Tfy
  Sf=sqrt(Xfa(i,j+1)**2+Yfa(i,j+1)**2)
  Dn(i,j)=Sf*interpl(Ga(i,j),Ga(i,j+1),daw(i,j+1))/dad(i,j+1)
  Gxf=interpl(Ga(i,j)*Fx(i,j),Ga(i,j+1)*Fx(i,j+1),daw(i,j+1))
  Gyf=interpl(Ga(i,j)*Fy(i,j),Ga(i,j+1)*Fy(i,j+1),daw(i,j+1))
  dncx=Xc(i,j+1)-Xc(i,j)
  dncy=Yc(i,j+1)-Yc(i,j)
  Tfx=Xfa(i,j+1)-Sf**2*dncx/(dncx*Xfa(i,j+1)+dncy*Yfa(i,j+1))
  Tfy=Yfa(i,j+1)-Sf**2*dncy/(dncx*Xfa(i,j+1)+dncy*Yfa(i,j+1))
  bno(i,j)=bno(i,j)+Gxf*Tfx+Gyf*Tfy
  Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   Ds(i,j)=Sf*interpl(Ga(Ic+1-i,j),Ga(i,j),daw(i,j))/dad(i,j)
   Gxf=interpl(Ga(Ic+1-i,j)*Fx(Ic+1-i,j),Ga(i,j)*Fx(i,j),daw(i,j))
   Gyf=interpl(Ga(Ic+1-i,j)*Fy(Ic+1-i,j),Ga(i,j)*Fy(i,j),daw(i,j))
   dncx=Xc(Ic+1-i,j)-Xc(i,j)
   dncy=Yc(Ic+1-i,j)-Yc(i,j)
  else if(j==1) then
   Ds(i,j)=Sf*Ga(i,j)/dad(i,j)
   Gxf=Ga(i,j)*Fx(i,j)
   Gyf=Ga(i,j)*Fy(i,j)
   dncx=Xw(i)-Xc(i,j)
   dncy=Yw(i)-Yc(i,j)
  else
   Ds(i,j)=Sf*interpl(Ga(i,j-1),Ga(i,j),daw(i,j))/dad(i,j)
   Gxf=interpl(Ga(i,j-1)*Fx(i,j-1),Ga(i,j)*Fx(i,j),daw(i,j))
   Gyf=interpl(Ga(i,j-1)*Fy(i,j-1),Ga(i,j)*Fy(i,j),daw(i,j))
   dncx=Xc(i,j-1)-Xc(i,j)
   dncy=Yc(i,j-1)-Yc(i,j)
  end if
  Tfx=-Xfa(i,j)-Sf**2*dncx/(-dncx*Xfa(i,j)-dncy*Yfa(i,j))
  Tfy=-Yfa(i,j)-Sf**2*dncy/(-dncx*Xfa(i,j)-dncy*Yfa(i,j))
  bno(i,j)=bno(i,j)+Gxf*Tfx+Gyf*Tfy
 end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
aM(1,:,:)=1
aM(2,:,:)=0
aM(3,:,:)=0
aM(4,:,:)=0
aM(5,:,:)=0
!$OMP END WORKSHARE
if(isU) then
 !$OMP WORKSHARE
 aupM=0
 !$OMP END WORKSHARE
else if(isV) then
 !$OMP WORKSHARE
 avpM=0
 !$OMP END WORKSHARE
end if
!$OMP DO PRIVATE(i,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,dkcw,dkce,dacs,dacn,DF)
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
    Fw(i,j)=rhok(i,j)*Unk(i,j)
    Fe(i,j)=rhok(i+1,j)*Unk(i+1,j)
    Fs(i,j)=rhoa(i,j)*Vna(i,j)
    Fn(i,j)=rhoa(i,j+1)*Vna(i,j+1)
    DF=Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)
    dkcw=dkw(i,j)
    dkce=dkw(i+1,j)
    dacs=daw(i,j)
    dacn=daw(i,j+1)
    if(isU) then
     aupM(2,i,j)=dkcw*Ygaw
     aupM(3,i,j)=-(1-dkce)*Ygae
     aupM(5,i,j)=(1-dacn)*Ygkn
    else if(isV) then
     avpM(2,i,j)=-dkcw*Xgaw
     avpM(3,i,j)=(1-dkce)*Xgae
     avpM(5,i,j)=-(1-dacn)*Xgkn
    end if
    aM(2,i,j)=Dw(i,j)+max(Fw(i,j),0.0)
    aM(3,i,j)=De(i,j)+max(-Fe(i,j),0.0)
    aM(5,i,j)=Dn(i,j)+max(-Fn(i,j),0.0)
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     if(isU) then
      aupM(4,i,j)=0
      aupM(1,i,j)=-((1-dkcw)*Ygaw-dkce*Ygae+dacn*Ygkn-Ygks)
     else if(isV) then
      avpM(4,i,j)=0
      avpM(1,i,j)=-(-(1-dkcw)*Xgaw+dkce*Xgae-dacn*Xgkn+Xgks)
     end if
     aM(4,i,j)=0
     aM(1,i,j)=aM(2,i,j)+aM(3,i,j)+aM(4,i,j)+aM(5,i,j)
     if((isSa.and.isLr).or.(isSst.and.isLr).or.isLam.or.isInv) then
      aM(1,i,j)=aM(1,i,j)+2*Ds(i,j)
     else if(isSa.and.isWf.or.(isSst.and.isWf).or.isKe) then
      if(isParvel) then
       if(isU) aM(1,i,j)=aM(1,i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*(Yfa(i,j)/DR(i))**2
       if(isV) aM(1,i,j)=aM(1,i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*(Xfa(i,j)/DR(i))**2
      else
       aM(1,i,j)=aM(1,i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)
      end if
     end if
    else
     if(isU) then
      aupM(4,i,j)=-dacs*Ygks
      aupM(1,i,j)=-((1-dkcw)*Ygaw-dkce*Ygae+dacn*Ygkn-(1-dacs)*Ygks)
     else if(isV) then
      avpM(4,i,j)=dacs*Xgks
      avpM(1,i,j)=-(-(1-dkcw)*Xgaw+dkce*Xgae-dacn*Xgkn+(1-dacs)*Xgks)
     end if
     aM(4,i,j)=Ds(i,j)+max(Fs(i,j),0.0)
     aM(1,i,j)=aM(2,i,j)+aM(3,i,j)+aM(4,i,j)+aM(5,i,j)
    end if
    !if(DF>0) aM(1,i,j)=aM(1,i,j)+DF
  end DO
end DO
!$OMP END DO
Call Defercorrect(F,Fx,Fy,Fwall,cor,Fw,Fe,Fs,Fn)
!$OMP WORKSHARE
b=F
!$OMP END WORKSHARE
!$OMP DO PRIVATE(i,DF)
DO j=1,Jc-1
  DO i=Is,Ie
  if(isU) then
   if(isCom) then
    b(i,j)=Vol(i,j)*(-2*(muxx(i,j)+mvyx(i,j))/3+muxx(i,j)+mvxy(i,j))
   else if(isIncom) then
    b(i,j)=Vol(i,j)*(muxx(i,j)+mvxy(i,j))
   end if
   if(isParvel.and.j==1.and.(i>=Ib1.and.i<=Ib2).and.(isSa.and.isWf.or.(isSst.and.isWf).or.isKe)) then
    b(i,j)=b(i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*Xfa(i,j)*Yfa(i,j)*V(i,j)/DR(i)**2
   end if 
  else if(isV) then
   if(isCom) then
    b(i,j)=Vol(i,j)*(-2*(muxy(i,j)+mvyy(i,j))/3+muyx(i,j)+mvyy(i,j))
   else if(isIncom) then
    b(i,j)=Vol(i,j)*(muyx(i,j)+mvyy(i,j))
   end if
   if(isParvel.and.j==1.and.(i>=Ib1.and.i<=Ib2).and.(isSa.and.isWf.or.(isSst.and.isWf).or.isKe)) then
    b(i,j)=b(i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*Xfa(i,j)*Yfa(i,j)*U(i,j)/DR(i)**2
   end if
  end if
  b(i,j)=b(i,j)+bno(i,j)+cor(i,j)
  !DF=Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)
  !if(DF<0) b(i,j)=b(i,j)-DF*F(i,j)
  end DO
end DO
!$OMP END DO
if(isU) then
 !$OMP WORKSHARE
 auM(1:5,:,:)=aM
 bu=b
 !$OMP END WORKSHARE
else if(isV) then
 !$OMP WORKSHARE
 auM(6,:,:)=aM(1,:,:)
 bv=b
 !$OMP END WORKSHARE
end if
!$OMP END PARALLEL
end Subroutine pUucoe
