Subroutine Condiff(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Dnow,Dnoe,Dnos,Dnon,Faw,Fae,Fks,Fkn,Fwallw,Fwalle
real(8) Xi,fnu1,fnu2,Sv,rm,gm,Ret,Dwplus,phi1,F1,betaistar,Fmt,alphaf,betai,Tmin,Dampk,Ymax,Ym
real(8) aP,aW,aE,aS,aN,DF
real(8),external:: interpl
real(8) Fwall(Ib1:Ib2)
real(8) F(Ic,Jc),Ga(Ic,Jc)
real(8) Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc),Dw(Ic,Jc),De(Ic,Jc),Ds(Ic,Jc),Dn(Ic,Jc),bno(Ic,Jc),cor(Ic,Jc)
real(8) St(Ic,Jc),Sm(Ic,Jc),fw1(Ic,Jc),alphastar(Ic,Jc),betastar(Ic,Jc),alpha(Ic,Jc),beta(Ic,Jc),Dwt(Ic,Jc),C3e(Ic,Jc)
character(*) scalar
character(6) wallfunutype,wallfunktype
logical(1) productlimit,sstlowre,sstcom,saprodlimit

wallfunutype='parvel'
wallfunktype='loglaw'
productlimit=.false.
sstlowre=.false.
sstcom=.false.
saprodlimit=.false.
Ym=11.225

!$OMP PARALLEL
if(scalar=='U'.or.scalar=='V') then
 !$OMP WORKSHARE
 Fwall=0
 Ga=mu+mut
 !$OMP END WORKSHARE
 if(scalar=='U') then
  !$OMP WORKSHARE
  F=U
  !$OMP END WORKSHARE
 else if(scalar=='V') then
  !$OMP WORKSHARE
  F=V
  !$OMP END WORKSHARE
 end if
else if(scalar=='T') then
 !$OMP WORKSHARE
 F=T
 Ga=ka/ca+mut/Prt
 !$OMP END WORKSHARE
 if(Tmptype=='fixed') then
  !$OMP WORKSHARE
  Fwall=Tf
  !$OMP END WORKSHARE
 else if(Tmptype=='flux') then
  if(Walltreat=='wf') then
   !$OMP WORKSHARE
   Fwall=T(Ib1:Ib2,1)+Qf*Tplus/(ca*rho(Ib1:Ib2,1)*ustar)
   !$OMP END WORKSHARE
  else
   !$OMP WORKSHARE
   Fwall=T(Ib1:Ib2,1)+Qf*Yp/(ca*Ga(Ib1:Ib2,1))
   !$OMP END WORKSHARE
  end if
  !$OMP SINGLE
  Tf=Fwall((Ib1+Ib2)/2)
  !$OMP END SINGLE
 end if
else if(scalar=='Tn') then
 !$OMP WORKSHARE
 F=Tn
 Fwall=1.5*Tn(Ib1:Ib2,1)-0.5*Tn(Ib1:Ib2,2)
 !Fwall=0
 Ga=(mu+rho*Tn)/sigman
 !$OMP END WORKSHARE
else if(scalar=='Tk'.and.Turmod=='ke') then
 !$OMP WORKSHARE
 F=Tk
 Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
 Ga=mu+mut/sigmak
 !$OMP END WORKSHARE
else if(scalar=='Tk'.and.Turmod=='sst') then
 !$OMP WORKSHARE
 F=Tk
 Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
 !Fwall=Tk(Ib1:Ib2,1)
 !Fwall=0
 Ga=mu+mut/sigmatk
 !$OMP END WORKSHARE
else if(scalar=='Te') then
 !$OMP WORKSHARE
 F=Te
 Fwall=1.5*Te(Ib1:Ib2,1)-0.5*Te(Ib1:Ib2,2)
 Ga=mu+mut/sigmae
 !$OMP END WORKSHARE
else if(scalar=='Tw') then
 !$OMP WORKSHARE
 F=Tw
 Fwall=1.5*Tw(Ib1:Ib2,1)-0.5*Tw(Ib1:Ib2,2)
 !Fwall=Tw(Ib1:Ib2,1)
 Ga=mu+mut/sigmatw
 !$OMP END WORKSHARE
end if
if(Turmod=='ke') then
!$OMP WORKSHARE
  Ymax=maxval(Ystar)
!$OMP END WORKSHARE
else
!$OMP WORKSHARE
  Ymax=maxval(Yplus)
!$OMP END WORKSHARE
end if
if(scalar=='T'.or.scalar=='Tk'.or.scalar=='Tw'.or.scalar=='Te') then
 !$OMP WORKSHARE
 St=sqrt(2*(Ux**2+Vy**2)+(Uy+Vx)**2)
 !$OMP END WORKSHARE
end if
if(Turmod=='sst'.and.(scalar=='Tk'.or.scalar=='Tw')) then
 !$OMP DO PRIVATE(Ret,Dwplus,phi1,F1,betaistar,Fmt,alphaf,betai,i)
 DO j=1,Jc-1
   DO i=Is,Ie
    Dwplus=max(2*rho(i,j)*(Tkx(i,j)*Twx(i,j)+Tky(i,j)*Twy(i,j))/(sigmaw2*Tw(i,j)),1e-10)
    phi1=min(max(sqrt(Tk(i,j))/(betastarf*Tw(i,j)*d(i,j)),500*mu(i,j)/(rho(i,j)*d(i,j)**2*Tw(i,j))),&
    4*rho(i,j)*Tk(i,j)/(sigmaw2*Dwplus*d(i,j)**2))
    F1=tanh(phi1**4)
    betai=F1*betai1+(1-F1)*betai2
    if(sstlowre) then
      Ret=rho(i,j)*Tk(i,j)/(mu(i,j)*Tw(i,j))
      alphastar(i,j)=alphastarf*(0.024+Ret/Rk)/(1+Ret/Rk)
      betaistar=betastarf*(4./15+(Ret/Rbeta)**4)/(1+(Ret/Rbeta)**4)
    else
      alphastar(i,j)=alphastarf
      betaistar=betastarf
    end if
    if(sstcom) then
      if(sqrt(2*Tk(i,j)/(gama*R*T(i,j)/Ma))<=Mt0) then
       Fmt=0
      else
       Fmt=2*Tk(i,j)/(gama*R*T(i,j)/Ma)-Mt0**2
      end if
    else
      Fmt=0
    end if
    betastar(i,j)=betaistar*(1+Zetastar*Fmt)
    alphaf=F1*(betai1/betastarf-kapa**2/(sigmaw1*sqrt(betastarf)))+&
    (1-F1)*(betai2/betastarf-kapa**2/(sigmaw2*sqrt(betastarf)))
    if(sstlowre) then
      alpha(i,j)=alphaf/alphastar(i,j)*(alpha0+Ret/Rw)/(1+Ret/Rw)
    else
      alpha(i,j)=alphaf/alphastar(i,j)
    end if
    beta(i,j)=betai*(1-betaistar/betai*Zetastar*Fmt)
    Dwt(i,j)=2*(1-F1)*rho(i,j)*(Tkx(i,j)*Twx(i,j)+Tky(i,j)*Twy(i,j))/(Tw(i,j)*sigmaw2)
   end DO
 end DO
 !$OMP END DO
else if(scalar=='Tn') then
 !$OMP DO PRIVATE(Xi,fnu1,fnu2,Sv,rm,gm,i)
 DO j=1,Jc-1
   DO i=Is,Ie
    if(Walltreat=='lr') then
     Xi=rho(i,j)*Tn(i,j)/mu(i,j)+Cks*ksi/d(i,j)
    else
     Xi=rho(i,j)*Tn(i,j)/mu(i,j)
    end if
    fnu1=Xi**3/(Xi**3+Cnu1**3)
    !if(Ymax>Ym) fnu1=1.0
    fnu2=1-Tn(i,j)/(mu(i,j)/rho(i,j)+Tn(i,j)*fnu1)
    !if(Ymax>Ym) fnu2=0.0
    Sv=abs(Uy(i,j)-Vx(i,j))
    if(saprodlimit) Sv=Sv+2.0*min(0.0,St(i,j)-Sv)
    Sm(i,j)=Sv+Tn(i,j)*fnu2/(kapa*d(i,j))**2
    rm=Tn(i,j)/(Sm(i,j)*(kapa*d(i,j))**2)
    gm=rm+Cw2*(rm**6-rm)
    fw1(i,j)=gm*((1+Cw3**6)/(gm**6+Cw3**6))**(1./6)
   end DO
 end DO
 !$OMP END DO
else if(scalar=='Te') then
 !$OMP WORKSHARE
 C3e=tanh(abs(V/U))
 !$OMP END WORKSHARE
end if
!$OMP DO PRIVATE(Dnow,Dnoe,Dnos,Dnon,Faw,Fae,Fks,Fkn,Fwallw,Fwalle,i)
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
   Faw=(F(i-1,j+1)+F(i,j+1)-F(Ic+2-i,j)-F(Ic+1-i,j))/(4*dy)
   Fae=(F(i+1,j+1)+F(i,j+1)-F(Ic-i,j)-F(Ic+1-i,j))/(4*dy)
   Fks=(F(Ic-i,j)+F(i+1,j)-F(Ic+2-i,j)-F(i-1,j))/(4*dx)
   Fkn=(F(i+1,j+1)+F(i+1,j)-F(i-1,j+1)-F(i-1,j))/(4*dx)
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
Tmin=minval(Tplus)
!$OMP END WORKSHARE
!$OMP WORKSHARE
aM(1,:,:)=1
aM(2,:,:)=0
aM(3,:,:)=0
aM(4,:,:)=0
aM(5,:,:)=0
!$OMP END WORKSHARE
!$OMP DO PRIVATE(aP,aW,aE,aS,aN,DF,i)
DO j=1,Jc-1
  DO i=Is,Ie
   Fw(i,j)=dy*rhok(i,j)*Unk(i,j)
   Fe(i,j)=dy*rhok(i+1,j)*Unk(i+1,j)
   Fs(i,j)=dx*rhoa(i,j)*Vna(i,j)
   Fn(i,j)=dx*rhoa(i,j+1)*Vna(i,j+1)
   DF=Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)
   aW=Dw(i,j)+max(Fw(i,j),0.0)
   aE=De(i,j)+max(-Fe(i,j),0.0)
   aN=Dn(i,j)+max(-Fn(i,j),0.0)
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    aS=0
    aP=aW+aE+aS+aN
    if(DF>0) aP=aP+DF
    !aP=aP+DF
    if((Turmod=='sa'.and.Walltreat=='lr').or.(Turmod=='sst'.and.Walltreat=='lr').or.Turmod=='lam'.or.Turmod=='inv') then
     if(scalar=='Tn'.and.ks(i)>0) then
      aP=aP+2*Ds(i,j)*Yp(i)/d(i,j)
     else if(scalar=='Tw') then
      aW=0
      aE=0
      aN=0
      aP=1
     else
      aP=aP+2*Ds(i,j)
     end if
     if(scalar=='Tn'.and.fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Jg(i,j)*dx*dy
     if(scalar=='Tk') aP=aP+rho(i,j)*betastar(i,j)*Tw(i,j)*Jg(i,j)*dx*dy
     if(scalar=='T'.and.Tmptype=='flux') aP=aP-2*Ds(i,j)
    else if(Turmod=='sa'.and.Walltreat=='wf'.or.(Turmod=='sst'.and.Walltreat=='wf').or.Turmod=='ke') then
     if(scalar=='U') then
      if(wallfunutype=='parvel') then
       aP=aP+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*(Yga(i,j)/da(i,j))**2
      else
       aP=aP+rho(i,j)*ustar(i)*DR(i)/Uplus(i)
      end if
     else if(scalar=='V') then
      if(wallfunutype=='parvel') then
       aP=aP+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*(Xga(i,j)/da(i,j))**2
      else
       aP=aP+rho(i,j)*ustar(i)*DR(i)/Uplus(i)
      end if
     else if(scalar=='T') then
      if(Tmptype=='fixed'.and.Tmin>=0) aP=aP+rho(i,j)*ustar(i)*DR(i)/Tplus(i)
     else if(scalar=='Tn') then
      aP=aP+2*Ds(i,j)
      if(Ymax>Ym.and.fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/(kapa*d(i,j))**2*Jg(i,j)*dx*dy
      if(Ymax<=Ym.and.fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Jg(i,j)*dx*dy
      !if(fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Jg(i,j)*dx*dy
     else if(Turmod=='sst'.and.scalar=='Tk') then
      if(wallfunktype=='genlaw') then
       aP=aP+rho(i,j)*ustar(i)**3*Uplus(i)*Jg(i,j)*dx*dy/(Tk(i,j)*Yp(i))
      else
       aP=aP+rho(i,j)*ustar(i)**3*Jg(i,j)*dx*dy/(Tk(i,j)*kapa*Yp(i))
      end if
     else if(Turmod=='ke'.and.scalar=='Tk') then
      if(wallfunktype=='genlaw') then
       aP=aP+rho(i,j)*ustar(i)**3*Uplus(i)*Jg(i,j)*dx*dy/(Tk(i,j)*Yp(i))
      else
       aP=aP+rho(i,j)*ustar(i)**3*Jg(i,j)*dx*dy/(Tk(i,j)*kapa*Yp(i))
      end if
     else if(scalar=='Te'.or.scalar=='Tw') then
      aW=0
      aE=0
      aN=0
      aP=1
     end if
    end if
   else
    aS=Ds(i,j)+max(Fs(i,j),0.0)
    aP=aW+aE+aS+aN
    if(DF>0) aP=aP+DF
    !aP=aP+DF
    if(scalar=='Tk'.and.Turmod=='ke') aP=aP+rho(i,j)*Te(i,j)*Jg(i,j)*dx*dy/Tk(i,j)
    if(scalar=='Te'.and.Turmod=='ke') aP=aP+C2e*rho(i,j)*Te(i,j)*Jg(i,j)*dx*dy/Tk(i,j)
    if(scalar=='Tn'.and.Walltreat=='lr'.and.fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Jg(i,j)*dx*dy
    if(scalar=='Tn'.and.Walltreat=='wf'.and.fw1(i,j)>=0) then
     if(Ymax>Ym) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/(kapa*d(i,j))**2*Jg(i,j)*dx*dy
     if(Ymax<=Ym) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Jg(i,j)*dx*dy
     !aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Jg(i,j)*dx*dy
    end if
    if(scalar=='Tk'.and.Turmod=='sst') aP=aP+rho(i,j)*betastar(i,j)*Tw(i,j)*Jg(i,j)*dx*dy
    if(scalar=='Tw'.and.Turmod=='sst') aP=aP+rho(i,j)*beta(i,j)*Tw(i,j)*Jg(i,j)*dx*dy
   end if
   aM(1,i,j)=aP
   aM(2,i,j)=aW
   aM(3,i,j)=aE
   aM(4,i,j)=aS
   aM(5,i,j)=aN
  end DO
end DO
!$OMP END DO
Call Defercorrect(F,Fwall,cor,Fw,Fe,Fs,Fn)
!$OMP WORKSHARE
b=F
!$OMP END WORKSHARE
!$OMP DO PRIVATE(i)
DO j=1,Jc-1
  DO i=Is,Ie
   if(scalar=='U') then
    if(Proctrl=='com') then
     b(i,j)=Jg(i,j)*(-Px(i,j)-2*(muxx(i,j)+mvyx(i,j))/3+muxx(i,j)+mvxy(i,j))*dx*dy
    else if(Proctrl=='incom') then
     b(i,j)=Jg(i,j)*(-Px(i,j)+muxx(i,j)+mvxy(i,j))*dx*dy
     !if(Turmod=='ke'.or.Turmod=='sst') b(i,j)=b(i,j)-Jg(i,j)*2*rho(i,j)*Tkx(i,j)/3*dx*dy
    end if
    if(wallfunutype=='parvel'.and.j==1.and.(i>=Ib1.and.i<=Ib2).and.(Turmod=='sa'.and.Walltreat=='wf'.or.Turmod=='sst'.and.Walltreat=='wf'.or.Turmod=='ke')) then
     b(i,j)=b(i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*Xga(i,j)*Yga(i,j)*V(i,j)/da(i,j)**2
    end if
   else if(scalar=='V') then
    if(Proctrl=='com') then
     b(i,j)=Jg(i,j)*(-Py(i,j)-2*(muxy(i,j)+mvyy(i,j))/3+muyx(i,j)+mvyy(i,j))*dx*dy
    else if(Proctrl=='incom') then
     b(i,j)=Jg(i,j)*(-Py(i,j)+muyx(i,j)+mvyy(i,j))*dx*dy
     !if(Turmod=='ke'.or.Turmod=='sst') b(i,j)=b(i,j)-Jg(i,j)*2*rho(i,j)*Tky(i,j)/3*dx*dy
    end if
    if(wallfunutype=='parvel'.and.j==1.and.(i>=Ib1.and.i<=Ib2).and.(Turmod=='sa'.and.Walltreat=='wf'.or.Turmod=='sst'.and.Walltreat=='wf'.or.Turmod=='ke')) then
     b(i,j)=b(i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*Xga(i,j)*Yga(i,j)*U(i,j)/da(i,j)**2
    end if
   else if(scalar=='T') then
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     if(Tmptype=='fixed') then
       if((Turmod=='sa'.and.Walltreat=='lr').or.(Turmod=='sst'.and.Walltreat=='lr').or.Turmod=='lam'.or.Turmod=='inv') then
        if(Proctrl=='com') then
         if(visheat=='Y') then
          b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j)-2*(mu(i,j)+mut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+&
          (mu(i,j)+mut(i,j))*St(i,j)**2)*dx*dy/ca+2*Ds(i,j)*Tf
         else if(visheat=='N') then
          b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))*dx*dy/ca+2*Ds(i,j)*Tf
         end if
        else if(Proctrl=='incom') then
         if(visheat=='Y') then
          b(i,j)=Jg(i,j)*(mu(i,j)+mut(i,j))*St(i,j)**2*dx*dy/ca+2*Ds(i,j)*Tf
         else if(visheat=='N') then
          b(i,j)=2*Ds(i,j)*Tf
         end if
        end if
       else if(Turmod=='sst'.and.Walltreat=='wf'.or.(Turmod=='sa'.and.Walltreat=='wf').or.Turmod=='ke') then
        if(Proctrl=='com') then
         b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))*dx*dy/ca+rho(i,j)*ustar(i)*Tf*DR(i)/Tplus(i)
         !b(i,j)=rho(i,j)*ustar(i)*Tf*DR(i)/Tplus(i)
         if(visheat=='Y'.and.(Turmod=='sst'.or.Ymax<=Ym)) b(i,j)=b(i,j)+Jg(i,j)*(-2*(mu(i,j)+mut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+(mu(i,j)+mut(i,j))*St(i,j)**2)*dx*dy/ca
        else if(Proctrl=='incom') then
         b(i,j)=rho(i,j)*ustar(i)*Tf*DR(i)/Tplus(i)
         if(visheat=='Y'.and.(Turmod=='sst'.or.Ymax<=Ym)) b(i,j)=b(i,j)+Jg(i,j)*(mu(i,j)+mut(i,j))*St(i,j)**2*dx*dy/ca
        end if
        if(Tmin<0) b(i,j)=b(i,j)-rho(i,j)*ustar(i)*T(i,j)*DR(i)/Tplus(i)
       end if
     else if(Tmptype=='flux') then
       if((Turmod=='sa'.and.Walltreat=='lr').or.(Turmod=='sst'.and.Walltreat=='lr').or.Turmod=='lam'.or.Turmod=='inv') then
        if(Proctrl=='com') then
         if(visheat=='Y') then
          b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j)-2*(mu(i,j)+mut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+&
          (mu(i,j)+mut(i,j))*St(i,j)**2)*dx*dy/ca+Qf*DR(i)/ca
         else if(visheat=='N') then
          b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))*dx*dy/ca+Qf*DR(i)/ca
         end if
        else if(Proctrl=='incom') then
         if(visheat=='Y') then
          b(i,j)=Jg(i,j)*(mu(i,j)+mut(i,j))*St(i,j)**2*dx*dy/ca+Qf*DR(i)/ca
         else if(visheat=='N') then
          b(i,j)=Qf*DR(i)/ca
         end if
        end if
       else if(Turmod=='sst'.and.Walltreat=='wf'.or.(Turmod=='sa'.and.Walltreat=='wf').or.Turmod=='ke') then
         if(Proctrl=='com') then
          b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))*dx*dy/ca+Qf*DR(i)/ca
         else if(Proctrl=='incom') then
          b(i,j)=Qf*DR(i)/ca
         end if
       end if
     end if
    else
     if(Proctrl=='com') then
      if(visheat=='Y') then
       b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j)-2*(mu(i,j)+mut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+&
       (mu(i,j)+mut(i,j))*St(i,j)**2)*dx*dy/ca
      else if(visheat=='N') then
       b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))*dx*dy/ca
      end if
     else if(Proctrl=='incom') then
      if(visheat=='Y') then
       b(i,j)=Jg(i,j)*(mu(i,j)+mut(i,j))*St(i,j)**2*dx*dy/ca
      else if(visheat=='N') then
       b(i,j)=0
      end if
     end if
    end if
   else if(scalar=='Tn') then
    b(i,j)=Cb2*rho(i,j)*(Tnx(i,j)**2+Tny(i,j)**2)*Jg(i,j)*dx*dy/sigman+rho(i,j)*Cb1*Sm(i,j)*Tn(i,j)*Jg(i,j)*dx*dy
    if(fw1(i,j)<0) then
     if(Walltreat=='wf'.and.Ymax>Ym) then
      b(i,j)=b(i,j)-rho(i,j)*Cw1*fw1(i,j)*(Tn(i,j)/(kapa*d(i,j)))**2*Jg(i,j)*dx*dy
     else
      b(i,j)=b(i,j)-rho(i,j)*Cw1*fw1(i,j)*(Tn(i,j)/d(i,j))**2*Jg(i,j)*dx*dy
     end if
    end if
   else if(scalar=='Tk'.and.Turmod=='ke') then
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     if(wallfunktype=='genlaw') then
      if(wallfunutype=='parvel') then
       b(i,j)=rho(i,j)*ustar(i)*(Un(i,j)/da(i,j))**2*Jg(i,j)*dx*dy/(Uplus(i)*Yp(i))
      else
       b(i,j)=rho(i,j)*ustar(i)*(U(i,j)**2+V(i,j)**2)*Jg(i,j)*dx*dy/(Uplus(i)*Yp(i))
      end if       
     else
      if(wallfunutype=='parvel') then
       b(i,j)=rho(i,j)*ustar(i)*(Un(i,j)/da(i,j))**2*Jg(i,j)*dx*dy/(Uplus(i)**2*kapa*Yp(i))
      else
       b(i,j)=rho(i,j)*ustar(i)*(U(i,j)**2+V(i,j)**2)*Jg(i,j)*dx*dy/(Uplus(i)**2*kapa*Yp(i))
      end if
     end if
    else
     b(i,j)=(mu(i,j)+mut(i,j))*St(i,j)**2*Jg(i,j)*dx*dy
     !b(i,j)=b(i,j)-rho(i,j)*Te(i,j)*Jg(i,j)*dx*dy
    end if
    if(Proctrl=='com') b(i,j)=b(i,j)-2*rho(i,j)*Te(i,j)*Tk(i,j)*Jg(i,j)*dx*dy/(gama*R*T(i,j)/Ma)
   else if(scalar=='Tk'.and.Turmod=='sst') then
    Dampk=rho(i,j)*betastar(i,j)*Tk(i,j)*Tw(i,j)*Jg(i,j)*dx*dy
    if(j==1.and.(i>=Ib1.and.i<=Ib2).and.Walltreat=='wf') then
     if(wallfunktype=='genlaw') then
      if(wallfunutype=='parvel') then
       b(i,j)=rho(i,j)*ustar(i)*(Un(i,j)/da(i,j))**2*Jg(i,j)*dx*dy/(Uplus(i)*Yp(i))
      else
       b(i,j)=rho(i,j)*ustar(i)*(U(i,j)**2+V(i,j)**2)*Jg(i,j)*dx*dy/(Uplus(i)*Yp(i))
      end if
     else
      if(wallfunutype=='parvel') then
       b(i,j)=rho(i,j)*ustar(i)*(Un(i,j)/da(i,j))**2*Jg(i,j)*dx*dy/(Uplus(i)**2*kapa*Yp(i))
      else
       b(i,j)=rho(i,j)*ustar(i)*(U(i,j)**2+V(i,j)**2)*Jg(i,j)*dx*dy/(Uplus(i)**2*kapa*Yp(i))
      end if
     end if
    else
     if(productlimit) then
      b(i,j)=min(mut(i,j)*St(i,j)**2*Jg(i,j)*dx*dy,10*Dampk)
     else
      b(i,j)=mut(i,j)*St(i,j)**2*Jg(i,j)*dx*dy
     end if
     !b(i,j)=b(i,j)-Dampk
    end if
   else if(scalar=='Te') then
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     cycle
    else
     b(i,j)=Te(i,j)*C1e*(mu(i,j)+mut(i,j))*St(i,j)**2*Jg(i,j)*dx*dy/Tk(i,j)
     !b(i,j)=b(i,j)-C2e*rho(i,j)*Te(i,j)*Jg(i,j)*dx*dy/Tk(i,j)
    end if
   else if(scalar=='Tw') then
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     cycle
    else
     b(i,j)=rho(i,j)*alpha(i,j)*alphastar(i,j)*St(i,j)**2*Jg(i,j)*dx*dy+Dwt(i,j)*Jg(i,j)*dx*dy
     !b(i,j)=b(i,j)-rho(i,j)*beta(i,j)*Tw(i,j)**2*Jg(i,j)*dx*dy
    end if
   end if
   if(scalar=='Tk'.or.scalar=='Te'.or.scalar=='Tw') then
    b(i,j)=b(i,j)+bno(i,j)
   else
    b(i,j)=b(i,j)+bno(i,j)+cor(i,j)
   end if
   !DF=Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)
   !if(DF<0) b(i,j)=b(i,j)-DF*F(i,j)
  end DO
end DO
!$OMP END DO
if(scalar=='U') then
 !$OMP WORKSHARE
 auP=aM(1,:,:)
 auNB=aM(2,:,:)+aM(3,:,:)+aM(4,:,:)+aM(5,:,:)
 !$OMP END WORKSHARE
end if
!$OMP END PARALLEL
end Subroutine Condiff
