Subroutine pUucoe(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,DF,Sf
real(8),external:: interpl
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

!$OMP DO PRIVATE(i,Sf)
DO j=1,Jc-1
 DO i=Is,Ie
  Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
  if(i==1) then
   Dw(i,j)=Sf*interpl(Ga(Ic,j),Ga(i,j),dkw(i,j))/dkd(i,j)
  else
   Dw(i,j)=Sf*interpl(Ga(i-1,j),Ga(i,j),dkw(i,j))/dkd(i,j)
  end if
  Sf=sqrt(Xfk(i+1,j)**2+Yfk(i+1,j)**2)
  if(i==Ic) then
   De(i,j)=Sf*interpl(Ga(i,j),Ga(1,j),dkw(i+1,j))/dkd(i+1,j)
  else
   De(i,j)=Sf*interpl(Ga(i,j),Ga(i+1,j),dkw(i+1,j))/dkd(i+1,j)
  end if
  Sf=sqrt(Xfa(i,j+1)**2+Yfa(i,j+1)**2)
  Dn(i,j)=Sf*interpl(Ga(i,j),Ga(i,j+1),daw(i,j+1))/dad(i,j+1)
  Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   Ds(i,j)=Sf*interpl(Ga(Ic+1-i,j),Ga(i,j),daw(i,j))/dad(i,j)
  else if(j==1) then
   Ds(i,j)=Sf*Ga(i,j)/dad(i,j)
  else
   Ds(i,j)=Sf*interpl(Ga(i,j-1),Ga(i,j),daw(i,j))/dad(i,j)
  end if
  bno(i,j)=0.0
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
!$OMP DO PRIVATE(i,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,DF)
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
    if(isU) then
     aupM(2,i,j)=0.5*Ygaw
     aupM(3,i,j)=-0.5*Ygae
     aupM(5,i,j)=0.5*Ygkn
    else if(isV) then
     avpM(2,i,j)=-0.5*Xgaw
     avpM(3,i,j)=0.5*Xgae
     avpM(5,i,j)=-0.5*Xgkn
    end if
    aM(2,i,j)=Dw(i,j)+max(Fw(i,j),0.0)
    aM(3,i,j)=De(i,j)+max(-Fe(i,j),0.0)
    aM(5,i,j)=Dn(i,j)+max(-Fn(i,j),0.0)
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     if(isU) then
      aupM(4,i,j)=0
      aupM(1,i,j)=-(aupM(2,i,j)+aupM(3,i,j)+aupM(4,i,j)+aupM(5,i,j)-Ygks)
     else if(isV) then
      avpM(4,i,j)=0
      avpM(1,i,j)=-(avpM(2,i,j)+avpM(3,i,j)+avpM(4,i,j)+avpM(5,i,j)+Xgks)
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
      aupM(4,i,j)=-0.5*Ygks
      aupM(1,i,j)=-(aupM(2,i,j)+aupM(3,i,j)+aupM(4,i,j)+aupM(5,i,j))
     else if(isV) then
      avpM(4,i,j)=0.5*Xgks
      avpM(1,i,j)=-(avpM(2,i,j)+avpM(3,i,j)+avpM(4,i,j)+avpM(5,i,j))
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
