Subroutine Condiff(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Sf,Xi,fnu1,fnu2,Sv,rm,gm,Ret,Dwplus,phi1,F1,betaistar,Fmt,alphaf,betai,Tmin,Dampk,Ymax,Ym
real(8) aP,aW,aE,aS,aN,DF
real(8),external:: interpl
real(8) Fwall(Ib1:Ib2)
real(8) F(Ic,Jc),Ga(Ic,Jc),Fx(Ic,Jc),Fy(Ic,Jc)
real(8) Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc),Dw(Ic,Jc),De(Ic,Jc),Ds(Ic,Jc),Dn(Ic,Jc),bno(Ic,Jc),cor(Ic,Jc)
real(8) St(Ic,Jc),Sm(Ic,Jc),fw1(Ic,Jc),alphastar(Ic,Jc),betastar(Ic,Jc),alpha(Ic,Jc),beta(Ic,Jc),Dwt(Ic,Jc),C3e(Ic,Jc)
character(*) scalar
character(6) wallfunutype,wallfunktype
logical(1) productlimit,sstlowre,sstcom,saprodlimit
logical(1) isU,isV,isT,isTn,isTk,isTe,isTw,isKe,isSst,isSa,isLam,isInv,isFixed,isFlux,isWf,isLr,isParvel,isGenlaw,&
isCom,isIncom,isVisheatY

wallfunutype='parvel'
wallfunktype='loglaw'
productlimit=.false.
sstlowre=.false.
sstcom=.false.
saprodlimit=.false.
Ym=11.225

isU = scalar=='U'
isV = scalar=='V'
isT = scalar=='T'
isTn = scalar=='Tn'
isTk = scalar=='Tk'
isTe = scalar=='Te'
isTw = scalar=='Tw'
isKe = Turmod=='ke'
isSst = Turmod=='sst'
isSa = Turmod=='sa'
isLam = Turmod=='lam'
isInv = Turmod=='inv'
isFixed = Tmptype=='fixed'
isFlux = Tmptype=='flux'
isWf = Walltreat=='wf'
isLr = Walltreat=='lr'
isParvel = wallfunutype=='parvel'
isGenlaw = wallfunktype=='genlaw'
isCom = Proctrl=='com'
isIncom = Proctrl=='incom'
isVisheatY = visheat=='Y'

!$OMP PARALLEL
if(isU.or.isV) then
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
else if(isT) then
 !$OMP WORKSHARE
 F=T
 Fx=Tx
 Fy=Ty
 Ga=ka/ca+mut/Prt
 !$OMP END WORKSHARE
 if(isFixed) then
  !$OMP WORKSHARE
  Fwall=Tf
  !$OMP END WORKSHARE
 else if(isFlux) then
  if(isWf) then
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
else if(isTn) then
 !$OMP WORKSHARE
 F=Tn
 Fx=Tnx
 Fy=Tny
 Fwall=1.5*Tn(Ib1:Ib2,1)-0.5*Tn(Ib1:Ib2,2)
 !Fwall=0
 Ga=(mu+rho*Tn)/sigman
 !$OMP END WORKSHARE
else if(isTk.and.isKe) then
 !$OMP WORKSHARE
 F=Tk
 Fx=Tkx
 Fy=Tky
 Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
 Ga=mu+mut/sigmak
 !$OMP END WORKSHARE
else if(isTk.and.isSst) then
 !$OMP WORKSHARE
 F=Tk
 Fx=Tkx
 Fy=Tky
 Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
 !Fwall=Tk(Ib1:Ib2,1)
 !Fwall=0
 Ga=mu+mut/sigmatk
 !$OMP END WORKSHARE
else if(isTe) then
 !$OMP WORKSHARE
 F=Te
 Fx=Tex
 Fy=Tey
 Fwall=1.5*Te(Ib1:Ib2,1)-0.5*Te(Ib1:Ib2,2)
 Ga=mu+mut/sigmae
 !$OMP END WORKSHARE
else if(isTw) then
 !$OMP WORKSHARE
 F=Tw
 Fx=Twx 
 Fy=Twy
 Fwall=1.5*Tw(Ib1:Ib2,1)-0.5*Tw(Ib1:Ib2,2)
 !Fwall=Tw(Ib1:Ib2,1)
 Ga=mu+mut/sigmatw
 !$OMP END WORKSHARE
end if
if(isKe) then
!$OMP WORKSHARE
  Ymax=maxval(Ystar)
!$OMP END WORKSHARE
else
!$OMP WORKSHARE
  Ymax=maxval(Yplus)
!$OMP END WORKSHARE
end if
if(isT.or.isTk.or.isTw.or.isTe) then
 !$OMP WORKSHARE
 St=sqrt(2*(Ux**2+Vy**2)+(Uy+Vx)**2)
 !$OMP END WORKSHARE
end if
if(isSst.and.(isTk.or.isTw)) then
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
else if(isTn) then
 !$OMP DO PRIVATE(Xi,fnu1,fnu2,Sv,rm,gm,i)
 DO j=1,Jc-1
   DO i=Is,Ie
    if(isLr) then
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
else if(isTe) then
 !$OMP WORKSHARE
 C3e=tanh(abs(V/U))
 !$OMP END WORKSHARE
end if
!$OMP DO PRIVATE(Sf,i)
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
   Fw(i,j)=rhok(i,j)*Unk(i,j)
   Fe(i,j)=rhok(i+1,j)*Unk(i+1,j)
   Fs(i,j)=rhoa(i,j)*Vna(i,j)
   Fn(i,j)=rhoa(i,j+1)*Vna(i,j+1)
   DF=Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)
   aW=Dw(i,j)+max(Fw(i,j),0.0)
   aE=De(i,j)+max(-Fe(i,j),0.0)
   aN=Dn(i,j)+max(-Fn(i,j),0.0)
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    aS=0
    aP=aW+aE+aS+aN
    if(DF>0) aP=aP+DF
    !aP=aP+DF
    if((isSa.and.isLr).or.(isSst.and.isLr).or.isLam.or.isInv) then
     if(isTn.and.ks(i)>0) then
      aP=aP+Ds(i,j)*Yp(i)/d(i,j)
     else if(isTw) then
      aW=0
      aE=0
      aN=0
      aP=1
     else
      aP=aP+Ds(i,j)
     end if
     if(isTn.and.fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Vol(i,j)
     if(isTk) aP=aP+rho(i,j)*betastar(i,j)*Tw(i,j)*Vol(i,j)
     if(isT.and.isFlux) aP=aP-Ds(i,j)
    else if(isSa.and.isWf.or.(isSst.and.isWf).or.isKe) then
     if(isU) then
      if(isParvel) then
       aP=aP+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*(Yfa(i,j)/DR(i))**2
      else
       aP=aP+rho(i,j)*ustar(i)*DR(i)/Uplus(i)
      end if
     else if(isV) then
      if(isParvel) then
       aP=aP+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*(Xfa(i,j)/DR(i))**2
      else
       aP=aP+rho(i,j)*ustar(i)*DR(i)/Uplus(i)
      end if
     else if(isT) then
      if(isFixed.and.Tmin>=0) aP=aP+rho(i,j)*ustar(i)*DR(i)/Tplus(i)
     else if(isTn) then
      aP=aP+Ds(i,j)
      if(Ymax>Ym.and.fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/(kapa*d(i,j))**2*Vol(i,j)
      if(Ymax<=Ym.and.fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Vol(i,j)
     else if(isSst.and.isTk) then
      if(isGenlaw) then
       aP=aP+rho(i,j)*ustar(i)**3*Uplus(i)*Vol(i,j)/(Tk(i,j)*Yp(i))
      else
       aP=aP+rho(i,j)*ustar(i)**3*Vol(i,j)/(Tk(i,j)*kapa*Yp(i))
      end if
     else if(isKe.and.isTk) then
      if(isGenlaw) then
       aP=aP+rho(i,j)*ustar(i)**3*Uplus(i)*Vol(i,j)/(Tk(i,j)*Yp(i))
      else
       aP=aP+rho(i,j)*ustar(i)**3*Vol(i,j)/(Tk(i,j)*kapa*Yp(i))
      end if
     else if(isTe.or.isTw) then
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
    if(isTk.and.isKe) aP=aP+rho(i,j)*Te(i,j)*Vol(i,j)/Tk(i,j)
    if(isTe.and.isKe) aP=aP+C2e*rho(i,j)*Te(i,j)*Vol(i,j)/Tk(i,j)
    if(isTn.and.isLr.and.fw1(i,j)>=0) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Vol(i,j)
    if(isTn.and.isWf.and.fw1(i,j)>=0) then
     if(Ymax>Ym) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/(kapa*d(i,j))**2*Vol(i,j)
     if(Ymax<=Ym) aP=aP+rho(i,j)*Cw1*fw1(i,j)*Tn(i,j)/d(i,j)**2*Vol(i,j)
    end if
    if(isTk.and.isSst) aP=aP+rho(i,j)*betastar(i,j)*Tw(i,j)*Vol(i,j)
    if(isTw.and.isSst) aP=aP+rho(i,j)*beta(i,j)*Tw(i,j)*Vol(i,j)
   end if
   aM(1,i,j)=aP
   aM(2,i,j)=aW
   aM(3,i,j)=aE
   aM(4,i,j)=aS
   aM(5,i,j)=aN
  end DO
end DO
!$OMP END DO
Call Defercorrect(F,Fx,Fy,Fwall,cor,Fw,Fe,Fs,Fn)
!$OMP WORKSHARE
b=F
!$OMP END WORKSHARE
!$OMP DO PRIVATE(Dampk,i)
DO j=1,Jc-1
  DO i=Is,Ie
   if(isU) then
    if(isCom) then
     b(i,j)=Vol(i,j)*(-Px(i,j)-2*(muxx(i,j)+mvyx(i,j))/3+muxx(i,j)+mvxy(i,j))
    else if(isIncom) then
     b(i,j)=Vol(i,j)*(-Px(i,j)+muxx(i,j)+mvxy(i,j))
     !if(isKe.or.isSst) b(i,j)=b(i,j)-Vol(i,j)*2*rho(i,j)*Tkx(i,j)/3
    end if
    if(isParvel.and.j==1.and.(i>=Ib1.and.i<=Ib2).and.(isSa.and.isWf.or.(isSst.and.isWf).or.isKe)) then
     b(i,j)=b(i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*Xfa(i,j)*Yfa(i,j)*V(i,j)/DR(i)**2
    end if
   else if(isV) then
    if(isCom) then
     b(i,j)=Vol(i,j)*(-Py(i,j)-2*(muxy(i,j)+mvyy(i,j))/3+muyx(i,j)+mvyy(i,j))
    else if(isIncom) then
     b(i,j)=Vol(i,j)*(-Py(i,j)+muyx(i,j)+mvyy(i,j))
     !if(isKe.or.isSst) b(i,j)=b(i,j)-Vol(i,j)*2*rho(i,j)*Tky(i,j)/3
    end if
    if(isParvel.and.j==1.and.(i>=Ib1.and.i<=Ib2).and.(isSa.and.isWf.or.(isSst.and.isWf).or.isKe)) then
     b(i,j)=b(i,j)+rho(i,j)*ustar(i)*DR(i)/Uplus(i)*Xfa(i,j)*Yfa(i,j)*U(i,j)/DR(i)**2
    end if
   else if(isT) then
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     if(isFixed) then
       if((isSa.and.isLr).or.(isSst.and.isLr).or.isLam.or.isInv) then
        if(isCom) then
         if(isVisheatY) then
          b(i,j)=Vol(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j)-2*(mu(i,j)+mut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+&
          (mu(i,j)+mut(i,j))*St(i,j)**2)/ca+Ds(i,j)*Tf
         else
          b(i,j)=Vol(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))/ca+Ds(i,j)*Tf
         end if
        else if(isIncom) then
         if(isVisheatY) then
          b(i,j)=Vol(i,j)*(mu(i,j)+mut(i,j))*St(i,j)**2/ca+Ds(i,j)*Tf
         else
          b(i,j)=Ds(i,j)*Tf
         end if
        end if
       else if(isSst.and.isWf.or.(isSa.and.isWf).or.isKe) then
        if(isCom) then
         b(i,j)=Vol(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))/ca+rho(i,j)*ustar(i)*Tf*DR(i)/Tplus(i)
         !b(i,j)=rho(i,j)*ustar(i)*Tf*DR(i)/Tplus(i)
         if(isVisheatY.and.(isSst.or.Ymax<=Ym)) b(i,j)=b(i,j)+Vol(i,j)*(-2*(mu(i,j)+mut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+(mu(i,j)+mut(i,j))*St(i,j)**2)/ca
        else if(isIncom) then
         b(i,j)=rho(i,j)*ustar(i)*Tf*DR(i)/Tplus(i)
         if(isVisheatY.and.(isSst.or.Ymax<=Ym)) b(i,j)=b(i,j)+Vol(i,j)*(mu(i,j)+mut(i,j))*St(i,j)**2/ca
        end if
        if(Tmin<0) b(i,j)=b(i,j)-rho(i,j)*ustar(i)*T(i,j)*DR(i)/Tplus(i)
       end if
     else if(isFlux) then
       if((isSa.and.isLr).or.(isSst.and.isLr).or.isLam.or.isInv) then
        if(isCom) then
         if(isVisheatY) then
          b(i,j)=Vol(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j)-2*(mu(i,j)+mut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+&
          (mu(i,j)+mut(i,j))*St(i,j)**2)/ca+Qf*DR(i)/ca
         else
          b(i,j)=Vol(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))/ca+Qf*DR(i)/ca
         end if
        else if(isIncom) then
         if(isVisheatY) then
          b(i,j)=Vol(i,j)*(mu(i,j)+mut(i,j))*St(i,j)**2/ca+Qf*DR(i)/ca
         else
          b(i,j)=Qf*DR(i)/ca
         end if
        end if
       else if(isSst.and.isWf.or.(isSa.and.isWf).or.isKe) then
         if(isCom) then
          b(i,j)=Vol(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))/ca+Qf*DR(i)/ca
         else if(isIncom) then
          b(i,j)=Qf*DR(i)/ca
         end if
       end if
     end if
    else
     if(isCom) then
      if(isVisheatY) then
       b(i,j)=Vol(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j)-2*(mu(i,j)+mut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+&
       (mu(i,j)+mut(i,j))*St(i,j)**2)/ca
      else
       b(i,j)=Vol(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))/ca
      end if
     else if(isIncom) then
      if(isVisheatY) then
       b(i,j)=Vol(i,j)*(mu(i,j)+mut(i,j))*St(i,j)**2/ca
      else
       b(i,j)=0
      end if
     end if
    end if
   else if(isTn) then
    b(i,j)=Cb2*rho(i,j)*(Tnx(i,j)**2+Tny(i,j)**2)*Vol(i,j)/sigman+rho(i,j)*Cb1*Sm(i,j)*Tn(i,j)*Vol(i,j)
    if(fw1(i,j)<0) then
     if(isWf.and.Ymax>Ym) then
      b(i,j)=b(i,j)-rho(i,j)*Cw1*fw1(i,j)*(Tn(i,j)/(kapa*d(i,j)))**2*Vol(i,j)
     else
      b(i,j)=b(i,j)-rho(i,j)*Cw1*fw1(i,j)*(Tn(i,j)/d(i,j))**2*Vol(i,j)
     end if
    end if
   else if(isTk.and.isKe) then
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     if(isGenlaw) then
      if(isParvel) then
       b(i,j)=rho(i,j)*ustar(i)*((U(i,j)*Yfa(i,j)-V(i,j)*Xfa(i,j))/DR(i))**2*Vol(i,j)/(Uplus(i)*Yp(i))
      else
       b(i,j)=rho(i,j)*ustar(i)*(U(i,j)**2+V(i,j)**2)*Vol(i,j)/(Uplus(i)*Yp(i))
      end if       
     else
      if(isParvel) then
       b(i,j)=rho(i,j)*ustar(i)*((U(i,j)*Yfa(i,j)-V(i,j)*Xfa(i,j))/DR(i))**2*Vol(i,j)/(Uplus(i)**2*kapa*Yp(i))
      else
       b(i,j)=rho(i,j)*ustar(i)*(U(i,j)**2+V(i,j)**2)*Vol(i,j)/(Uplus(i)**2*kapa*Yp(i))
      end if
     end if
    else
     b(i,j)=(mu(i,j)+mut(i,j))*St(i,j)**2*Vol(i,j)
     !b(i,j)=b(i,j)-rho(i,j)*Te(i,j)*Vol(i,j)
    end if
    if(isCom) b(i,j)=b(i,j)-2*rho(i,j)*Te(i,j)*Tk(i,j)*Vol(i,j)/(gama*R*T(i,j)/Ma)
   else if(isTk.and.isSst) then
    Dampk=rho(i,j)*betastar(i,j)*Tk(i,j)*Tw(i,j)*Vol(i,j)
    if(j==1.and.(i>=Ib1.and.i<=Ib2).and.isWf) then
     if(isGenlaw) then
      if(isParvel) then
       b(i,j)=rho(i,j)*ustar(i)*((U(i,j)*Yfa(i,j)-V(i,j)*Xfa(i,j))/DR(i))**2*Vol(i,j)/(Uplus(i)*Yp(i))
      else
       b(i,j)=rho(i,j)*ustar(i)*(U(i,j)**2+V(i,j)**2)*Vol(i,j)/(Uplus(i)*Yp(i))
      end if
     else
      if(isParvel) then
       b(i,j)=rho(i,j)*ustar(i)*((U(i,j)*Yfa(i,j)-V(i,j)*Xfa(i,j))/DR(i))**2*Vol(i,j)/(Uplus(i)**2*kapa*Yp(i))
      else
       b(i,j)=rho(i,j)*ustar(i)*(U(i,j)**2+V(i,j)**2)*Vol(i,j)/(Uplus(i)**2*kapa*Yp(i))
      end if
     end if
    else
     if(productlimit) then
      b(i,j)=min(mut(i,j)*St(i,j)**2*Vol(i,j),10*Dampk)
     else
      b(i,j)=mut(i,j)*St(i,j)**2*Vol(i,j)
     end if
     !b(i,j)=b(i,j)-Dampk
    end if
   else if(isTe) then
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     cycle
    else
     b(i,j)=Te(i,j)*C1e*(mu(i,j)+mut(i,j))*St(i,j)**2*Vol(i,j)/Tk(i,j)
     !b(i,j)=b(i,j)-C2e*rho(i,j)*Te(i,j)*Vol(i,j)/Tk(i,j)
    end if
   else if(isTw) then
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     cycle
    else
     b(i,j)=rho(i,j)*alpha(i,j)*alphastar(i,j)*St(i,j)**2*Vol(i,j)+Dwt(i,j)*Vol(i,j)
     !b(i,j)=b(i,j)-rho(i,j)*beta(i,j)*Tw(i,j)**2*Vol(i,j)
    end if
   end if
   if(isTk.or.isTe.or.isTw) then
    b(i,j)=b(i,j)+bno(i,j)
   else
    b(i,j)=b(i,j)+bno(i,j)+cor(i,j)
   end if
   !DF=Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)
   !if(DF<0) b(i,j)=b(i,j)-DF*F(i,j)
  end DO
end DO
!$OMP END DO
if(isU) then
 !$OMP WORKSHARE
 auP=aM(1,:,:)
 auNB=aM(2,:,:)+aM(3,:,:)+aM(4,:,:)+aM(5,:,:)
 !$OMP END WORKSHARE
end if
!$OMP END PARALLEL
!print *,'Completed assembling the coefficient matrix:',scalar,minval(aM(1,Is:Ie,1:Jc-1)),maxval(aM(1,Is:Ie,1:Jc-1))
end Subroutine Condiff
