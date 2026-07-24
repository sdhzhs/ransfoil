Subroutine pUucoe(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Dnow,Dnoe,Dnos,Dnon,Faw,Fae,Fks,Fkn,Fwallw,Fwalle,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,DF,dkc,dkw,dke,dac,das,dan
real(8) F(Ic,Jc),Ga(Ic,Jc),Fwall(Ib1:Ib2)
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
 !$OMP END WORKSHARE
else if(isV) then
 !$OMP WORKSHARE
 F=V
 !$OMP END WORKSHARE
end if

!$OMP DO PRIVATE(i,Dnow,Dnoe,Dnos,Dnon,Faw,Fae,Fks,Fkn,Fwallw,Fwalle)
DO j=1,Jc-1
  DO i=Is,Ie
   if(i==1) then
    Dw(i,j)=interpl(a1(i,j)*Ga(i,j),a1(Ic,j)*Ga(Ic,j),dk(i,j),dk(Ic,j))*dy/dx
    Dnow=interpl(b1(i,j)*Ga(i,j),b1(Ic,j)*Ga(Ic,j),dk(i,j),dk(Ic,j))*dy
   else
    Dw(i,j)=interpl(a1(i,j)*Ga(i,j),a1(i-1,j)*Ga(i-1,j),dk(i,j),dk(i-1,j))*dy/dx
    Dnow=interpl(b1(i,j)*Ga(i,j),b1(i-1,j)*Ga(i-1,j),dk(i,j),dk(i-1,j))*dy
   end if
   if(i==Ic) then
    De(i,j)=interpl(a1(i,j)*Ga(i,j),a1(1,j)*Ga(1,j),dk(i,j),dk(1,j))*dy/dx
    Dnoe=interpl(b1(i,j)*Ga(i,j),b1(1,j)*Ga(1,j),dk(i,j),dk(1,j))*dy
   else
    De(i,j)=interpl(a1(i,j)*Ga(i,j),a1(i+1,j)*Ga(i+1,j),dk(i,j),dk(i+1,j))*dy/dx
    Dnoe=interpl(b1(i,j)*Ga(i,j),b1(i+1,j)*Ga(i+1,j),dk(i,j),dk(i+1,j))*dy
   end if
   Dn(i,j)=interpl(y1(i,j)*Ga(i,j),y1(i,j+1)*Ga(i,j+1),da(i,j),da(i,j+1))*dx/dy
   Dnon=interpl(b1(i,j)*Ga(i,j),b1(i,j+1)*Ga(i,j+1),da(i,j),da(i,j+1))*dx
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
    Ds(i,j)=interpl(y1(i,j)*Ga(i,j),y1(Ic+1-i,j)*Ga(Ic+1-i,j),da(i,j),da(Ic+1-i,j))*dx/dy
    Dnos=interpl(b1(i,j)*Ga(i,j),b1(Ic+1-i,j)*Ga(Ic+1-i,j),da(i,j),da(Ic+1-i,j))*dx
    Faw=0.5*(F(i-1,j+1)+F(i,j+1)-F(Ic+2-i,j)-F(Ic+1-i,j))/(2*dy)
    Fae=0.5*(F(i+1,j+1)+F(i,j+1)-F(Ic-i,j)-F(Ic+1-i,j))/(2*dy)
    Fks=0.5*(F(Ic-i,j)+F(i+1,j)-F(Ic+2-i,j)-F(i-1,j))/(2*dx)
    Fkn=0.5*(F(i+1,j+1)+F(i+1,j)-F(i-1,j+1)-F(i-1,j))/(2*dx)
   else if(j==1) then
    Ds(i,j)=y1(i,j)*Ga(i,j)*dx/dy
    Dnos=b1(i,j)*Ga(i,j)*dx
    if(i==Ib1) then
     Fwallw=Fwall(i)
    else
     Fwallw=0.5*(Fwall(i)+Fwall(i-1))
    end if
    if(i==Ib2) then
     Fwalle=Fwall(i)
    else
     Fwalle=0.5*(Fwall(i)+Fwall(i+1))
    end if
    if(i==1) then
     Faw=(0.25*(F(Ic,j+1)+F(i,j+1)+F(Ic,j)+F(i,j))-Fwallw)/dy
    else
     Faw=(0.25*(F(i-1,j+1)+F(i,j+1)+F(i-1,j)+F(i,j))-Fwallw)/dy
    end if
    if(i==Ic) then
     Fae=(0.25*(F(1,j+1)+F(i,j+1)+F(1,j)+F(i,j))-Fwalle)/dy
    else
     Fae=(0.25*(F(i+1,j+1)+F(i,j+1)+F(i+1,j)+F(i,j))-Fwalle)/dy
    end if
    if(i==1) then
     Fkn=(F(i+1,j+1)+F(i+1,j)-F(Ic,j+1)-F(Ic,j))/(4*dx)
    else if(i==Ic) then
     Fkn=(F(1,j+1)+F(1,j)-F(i-1,j+1)-F(i-1,j))/(4*dx)
    else
     Fkn=(F(i+1,j+1)+F(i+1,j)-F(i-1,j+1)-F(i-1,j))/(4*dx)
    end if
    Fks=(Fwalle-Fwallw)/dx
   else if(i==1) then
    Ds(i,j)=interpl(y1(i,j)*Ga(i,j),y1(i,j-1)*Ga(i,j-1),da(i,j),da(i,j-1))*dx/dy
    Dnos=interpl(b1(i,j)*Ga(i,j),b1(i,j-1)*Ga(i,j-1),da(i,j),da(i,j-1))*dx
    Faw=(F(Ic,j+1)+F(i,j+1)-F(Ic,j-1)-F(i,j-1))/(4*dy)
    Fae=(F(i+1,j+1)+F(i,j+1)-F(i+1,j-1)-F(i,j-1))/(4*dy)
    Fks=(F(i+1,j-1)+F(i+1,j)-F(Ic,j-1)-F(Ic,j))/(4*dx)
    Fkn=(F(i+1,j+1)+F(i+1,j)-F(Ic,j+1)-F(Ic,j))/(4*dx)
   else if(i==Ic) then
    Ds(i,j)=interpl(y1(i,j)*Ga(i,j),y1(i,j-1)*Ga(i,j-1),da(i,j),da(i,j-1))*dx/dy
    Dnos=interpl(b1(i,j)*Ga(i,j),b1(i,j-1)*Ga(i,j-1),da(i,j),da(i,j-1))*dx
    Faw=(F(i-1,j+1)+F(i,j+1)-F(i-1,j-1)-F(i,j-1))/(4*dy)
    Fae=(F(1,j+1)+F(i,j+1)-F(1,j-1)-F(i,j-1))/(4*dy)
    Fks=(F(1,j-1)+F(1,j)-F(i-1,j-1)-F(i-1,j))/(4*dx)
    Fkn=(F(1,j+1)+F(1,j)-F(i-1,j+1)-F(i-1,j))/(4*dx)
   else
    Ds(i,j)=interpl(y1(i,j)*Ga(i,j),y1(i,j-1)*Ga(i,j-1),da(i,j),da(i,j-1))*dx/dy
    Dnos=interpl(b1(i,j)*Ga(i,j),b1(i,j-1)*Ga(i,j-1),da(i,j),da(i,j-1))*dx
    Faw=(F(i-1,j+1)+F(i,j+1)-F(i-1,j-1)-F(i,j-1))/(4*dy)
    Fae=(F(i+1,j+1)+F(i,j+1)-F(i+1,j-1)-F(i,j-1))/(4*dy)
    Fks=(F(i+1,j-1)+F(i+1,j)-F(i-1,j-1)-F(i-1,j))/(4*dx)
    Fkn=(F(i+1,j+1)+F(i+1,j)-F(i-1,j+1)-F(i-1,j))/(4*dx)
   end if
   bno(i,j)=Dnow*Faw-Dnoe*Fae+Dnos*Fks-Dnon*Fkn
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
!$OMP DO PRIVATE(i,Xgaw,Xgae,Ygaw,Ygae,Xgks,Xgkn,Ygks,Ygkn,dkc,dkw,dke,dac,das,dan,DF)
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
    Fw(i,j)=dy*rhok(i,j)*Unk(i,j)
    Fe(i,j)=dy*rhok(i+1,j)*Unk(i+1,j)
    Fs(i,j)=dx*rhoa(i,j)*Vna(i,j)
    Fn(i,j)=dx*rhoa(i,j+1)*Vna(i,j+1)
    DF=Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)
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
    if(isU) then
     aupM(2,i,j)=dkc*Ygaw*dy/(dkc+dkw)
     aupM(3,i,j)=-dkc*Ygae*dy/(dkc+dke)
     aupM(5,i,j)=dac*Ygkn*dx/(dac+dan)
    else if(isV) then
     avpM(2,i,j)=-dkc*Xgaw*dy/(dkc+dkw)
     avpM(3,i,j)=dkc*Xgae*dy/(dkc+dke)
     avpM(5,i,j)=-dac*Xgkn*dx/(dac+dan)
    end if
    aM(2,i,j)=Dw(i,j)+max(Fw(i,j),0.0)
    aM(3,i,j)=De(i,j)+max(-Fe(i,j),0.0)
    aM(5,i,j)=Dn(i,j)+max(-Fn(i,j),0.0)
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     if(isU) then
      aupM(4,i,j)=0
      aupM(1,i,j)=-(dkw*Ygaw*dy/(dkc+dkw)-dke*Ygae*dy/(dkc+dke)+dan*Ygkn*dx/(dac+dan)-Ygks*dx)
     else if(isV) then
      avpM(4,i,j)=0
      avpM(1,i,j)=-(-dkw*Xgaw*dy/(dkc+dkw)+dke*Xgae*dy/(dkc+dke)-dan*Xgkn*dx/(dac+dan)+Xgks*dx)
     end if
     aM(4,i,j)=0
     aM(1,i,j)=aM(2,i,j)+aM(3,i,j)+aM(4,i,j)+aM(5,i,j)
     if((isSa.and.isLr).or.(isSst.and.isLr).or.isLam.or.isInv) then
      aM(1,i,j)=aM(1,i,j)+2*Ds(i,j)
     else if(isSa.and.isWf.or.(isSst.and.isWf).or.isKe) then
      if(isParvel) then
       if(isU) aM(1,i,j)=aM(1,i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*(Yga(i,j)/da(i,j))**2
       if(isV) aM(1,i,j)=aM(1,i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*(Xga(i,j)/da(i,j))**2
      else
       aM(1,i,j)=aM(1,i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)
      end if
     end if
    else
     if(isU) then
      aupM(4,i,j)=-dac*Ygks*dx/(dac+das)
      aupM(1,i,j)=-(dkw*Ygaw*dy/(dkc+dkw)-dke*Ygae*dy/(dkc+dke)-das*Ygks*dx/(dac+das)+dan*Ygkn*dx/(dac+dan))
     else if(isV) then
      avpM(4,i,j)=dac*Xgks*dx/(dac+das)
      avpM(1,i,j)=-(-dkw*Xgaw*dy/(dkc+dkw)+dke*Xgae*dy/(dkc+dke)+das*Xgks*dx/(dac+das)-dan*Xgkn*dx/(dac+dan))
     end if
     aM(4,i,j)=Ds(i,j)+max(Fs(i,j),0.0)
     aM(1,i,j)=aM(2,i,j)+aM(3,i,j)+aM(4,i,j)+aM(5,i,j)
    end if
    !if(DF>0) aM(1,i,j)=aM(1,i,j)+DF
  end DO
end DO
!$OMP END DO
Call Defercorrect(F,Fwall,cor,Fw,Fe,Fs,Fn)
!$OMP WORKSHARE
b=F
!$OMP END WORKSHARE
!$OMP DO PRIVATE(i,DF)
DO j=1,Jc-1
  DO i=Is,Ie
  if(isU) then
   if(isCom) then
    b(i,j)=Jg(i,j)*(-2*(muxx(i,j)+mvyx(i,j))/3+muxx(i,j)+mvxy(i,j))*dx*dy
   else if(isIncom) then
    b(i,j)=Jg(i,j)*(muxx(i,j)+mvxy(i,j))*dx*dy
   end if
   if(isParvel.and.j==1.and.(i>=Ib1.and.i<=Ib2).and.(isSa.and.isWf.or.(isSst.and.isWf).or.isKe)) then
    b(i,j)=b(i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*Xga(i,j)*Yga(i,j)*V(i,j)/da(i,j)**2
   end if 
  else if(isV) then
   if(isCom) then
    b(i,j)=Jg(i,j)*(-2*(muxy(i,j)+mvyy(i,j))/3+muyx(i,j)+mvyy(i,j))*dx*dy
   else if(isIncom) then
    b(i,j)=Jg(i,j)*(muyx(i,j)+mvyy(i,j))*dx*dy
   end if
   if(isParvel.and.j==1.and.(i>=Ib1.and.i<=Ib2).and.(isSa.and.isWf.or.(isSst.and.isWf).or.isKe)) then
    b(i,j)=b(i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*Xga(i,j)*Yga(i,j)*U(i,j)/da(i,j)**2
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
