Subroutine Wallfunc
use Aero2DCOM
implicit none
integer i,j
real(8) Ym,Yt,Ymax
real(8) Dwplus,phi1,F1,Tplusl,Tplust,Tplusc,Dl,Dt
real(8) ksplus(Ib1:Ib2),deltaB(Ib1:Ib2),Prough(Ib1:Ib2),lamda(Ib1:Ib2),Uplusl(Ib1:Ib2),Uplust(Ib1:Ib2),&
Twplusl(Ib1:Ib2),Twplust(Ib1:Ib2),Twplus(Ib1:Ib2),betai(Ib1:Ib2)
character(6) wallfunutype,wallfunktype
logical(1) isKe,isSst,isSa,isFixed,isFlux,isWf,isLr,isParvel,isGenlaw,isLoglaw,isVisheatY

wallfunutype='parvel'
wallfunktype='loglaw'
Ym=11.225
Yt=11.8

isKe = Turmod=='ke'
isSst = Turmod=='sst'
isSa = Turmod=='sa'
isFixed = Tmptype=='fixed'
isFlux = Tmptype=='flux'
isWf = Walltreat=='wf'
isLr = Walltreat=='lr'
isParvel = wallfunutype=='parvel'
isGenlaw = wallfunktype=='genlaw'
isLoglaw = wallfunktype=='loglaw'
isVisheatY = visheat=='Y'

!$OMP PARALLEL
!$OMP DO
DO i=Ib1,Ib2
  if(isSst) then
   if(isParvel) then
    ustar(i)=((mu(i,1)*abs(U(i,1)*Yfa(i,1)-V(i,1)*Xfa(i,1))/DR(i)/(rho(i,1)*Yp(i)))**2+(alpha1*Tk(i,1))**2)**0.25
   else
    ustar(i)=((mu(i,1)*sqrt(U(i,1)**2+V(i,1)**2)/(rho(i,1)*Yp(i)))**2+(alpha1*Tk(i,1))**2)**0.25
   end if
  else if(isSa) then
   if(isParvel) then
    ustar(i)=((mu(i,1)*abs(U(i,1)*Yfa(i,1)-V(i,1)*Xfa(i,1))/DR(i)/(rho(i,1)*Yp(i)))**2+(mut(i,1)/(rho(i,1)*kapa*Yp(i)))**4)**0.25
   else
    ustar(i)=((mu(i,1)*sqrt(U(i,1)**2+V(i,1)**2)/(rho(i,1)*Yp(i)))**2+(mut(i,1)/(rho(i,1)*kapa*Yp(i)))**4)**0.25
   end if
  else if(isKe) then
   ustar(i)=Cu**0.25*Tk(i,1)**0.5
  end if
  ksplus(i)=rho(i,1)*ks(i)*ustar(i)/mu(i,1)
  if(isSa.or.isSst) then
   Yplus(i)=rho(i,1)*Yp(i)*ustar(i)/mu(i,1)
   if(isWf.and.ks(i)>0) then
    Yplus(i)=max(Yplus(i),ksplus(i)/2,2.5)
   end if
  else if(isKe) then
   Ystar(i)=rho(i,1)*Yp(i)*ustar(i)/mu(i,1)
   Ystar(i)=max(Ystar(i),ksplus(i)/2,Ym)
  end if
  if(ksplus(i)<=2.25) then
   deltaB(i)=0
  else if(ksplus(i)>2.25.and.ksplus(i)<=90) then
   deltaB(i)=log((ksplus(i)-2.25)/87.25+Cks*ksplus(i))*sin(0.4258*(log(ksplus(i))-0.811))/kapa
  else
   deltaB(i)=log(1+Cks*ksplus(i))/kapa
  end if
  if(ks(i)>0) then
   Prough(i)=3.15*Pr(i,1)**0.695*(exp(kapa*deltaB(i))/Ep-1/Ep)**0.359+(1/exp(kapa*deltaB(i)))**0.6*Pc(i,1)
  else
   Prough(i)=Pc(i,1)
  end if
end DO
!$OMP END DO
if((isSa.or.isSst).and.isWf) then
 !$OMP WORKSHARE
 lamda=-0.01*Yplus**4/(1+5.*Yplus)
 Uplusl=Yplus
 Uplust=log(Ep*Yplus)/kapa-deltaB
 Uplus=exp(lamda)*Uplusl+exp(1./lamda)*Uplust
 lamda=-0.01*(Pr(Ib1:Ib2,1)*Yplus)**4/(1+5.*Pr(Ib1:Ib2,1)**3*Yplus)
 Ymax=maxval(Yplus)
 !$OMP END WORKSHARE
 !$OMP DO PRIVATE(Tplusl,Tplust,Tplusc,Dl,Dt)
 DO i=Ib1,Ib2
  if(isVisheatY.and.Ymax>Ym) then
   Tplusl=Pr(i,1)*Uplusl(i)
   Tplust=Prt*(Uplust(i)+Prough(i))
   !Tplust=Prt*Uplust(i)
   Tplusc=exp(lamda(i))*Tplusl+exp(1./lamda(i))*Tplust
   if(isFixed) then
    Dl=0.5*rho(i,1)*ustar(i)*Pr(i,1)*(U(i,1)**2+V(i,1)**2)
    Dt=0.5*rho(i,1)*ustar(i)*(Prt*(U(i,1)**2+V(i,1)**2)+(Pr(i,1)-Prt)*(Yt*ustar(i))**2)
    !Dt=0.5*rho(i,1)*ustar(i)*Prt*(U(i,1)**2+V(i,1)**2)
    Q(i)=(ca*rho(i,1)*ustar(i)*(Tf-T(i,1))-exp(lamda(i))*Dl-exp(1./lamda(i))*Dt)/Tplusc
    Tplusl=Tplusl+Dl/Q(i)
    Tplust=Tplust+Dt/Q(i)
    Tplus(i)=exp(lamda(i))*Tplusl+exp(1./lamda(i))*Tplust
   else if(isFlux) then
    Q(i)=Qf
    Tplus(i)=Tplusc
   end if
  else
   Tplusl=Pr(i,1)*Uplusl(i)
   Tplust=Prt*(Uplust(i)+Prough(i))
   Tplus(i)=exp(lamda(i))*Tplusl+exp(1./lamda(i))*Tplust
   if(isFixed) then
    Q(i)=ca*rho(i,1)*ustar(i)*(Tf-T(i,1))/Tplus(i)
   else if(isFlux) then
    Q(i)=Qf
   end if
  end if
 end DO
 !$OMP END DO
else if(isKe) then
 !$OMP DO PRIVATE(Tplusc,Dt)
 DO i=Ib1,Ib2
  Uplus(i)=log(Ep*Ystar(i))/kapa-deltaB(i)
  !if(Ystar(i)<Ym) Uplus(i)=Ystar(i)
  if(isVisheatY) then
   Tplusc=Prt*(Uplus(i)+Prough(i))
   !Tplusc=Prt*Uplus(i)
   if(isFixed) then
    Dt=0.5*rho(i,1)*ustar(i)*(Prt*(U(i,1)**2+V(i,1)**2)+(Pr(i,1)-Prt)*(Yt*ustar(i))**2)
    !Dt=0.5*rho(i,1)*ustar(i)*Prt*(U(i,1)**2+V(i,1)**2)
    Q(i)=(ca*rho(i,1)*ustar(i)*(Tf-T(i,1))-Dt)/Tplusc
    Tplus(i)=Tplusc+Dt/Q(i)
   else if(isFlux) then
    Q(i)=Qf
    Tplus(i)=Tplusc
   end if
  else
   Tplus(i)=Prt*(Uplus(i)+Prough(i))
   !if(Ystar(i)<Yt) Tplus(i)=Pr(i,1)*Uplus(i)
   if(isFixed) then
    Q(i)=ca*rho(i,1)*ustar(i)*(Tf-T(i,1))/Tplus(i)
   else if(isFlux) then
    Q(i)=Qf
   end if
  end if
 end DO
 !$OMP END DO
end if
if(isSst) then
 !$OMP DO PRIVATE(Dwplus,phi1,F1)
 DO j=1,Jc
  DO i=1,Ic
   Dwplus=max(2*rho(i,j)*(Tkx(i,j)*Twx(i,j)+Tky(i,j)*Twy(i,j))/(sigmaw2*Tw(i,j)),1e-10)
   phi1=min(max(sqrt(Tk(i,j))/(betastarf*Tw(i,j)*d(i,j)),500*mu(i,j)/(rho(i,j)*d(i,j)**2*Tw(i,j))),&
   4*rho(i,j)*Tk(i,j)/(sigmaw2*Dwplus*d(i,j)**2))
   F1=tanh(phi1**4)
   sigmatk(i,j)=1/(F1/sigmak1+(1-F1)/sigmak2)
   sigmatw(i,j)=1/(F1/sigmaw1+(1-F1)/sigmaw2)
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    betai(i)=F1*betai1+(1-F1)*betai2
   end if
  end DO
 end DO
 !$OMP END DO
 if(isWf) then
  !$OMP WORKSHARE
  lamda=-0.01*Yplus**4/(1+5.*Yplus)
  Twplusl=6./(betai*Yplus**2)
  Twplust=1/(sqrt(betastarf)*kapa*Yplus)
  Twplus=exp(lamda)*Twplusl+exp(1./lamda)*Twplust
  !$OMP END WORKSHARE
 else if(isLr) then
  !$OMP DO
  DO i=Ib1,Ib2
   if(ksplus(i)>0) then
    ksplus(i)=max(1.0,ksplus(i))
    if(ksplus(i)<25) then
     Twplusl(i)=(50/ksplus(i))**2
    else if(ksplus(i)>=25) then
     Twplusl(i)=100/ksplus(i)
    end if
    Twplus(i)=min(Twplusl(i),60./(betai(i)*Yplus(i)**2))
   else
    Twplus(i)=60./(betai(i)*Yplus(i)**2)
   end if
  end DO
  !$OMP END DO
 end if
 !$OMP WORKSHARE
 Tw(Ib1:Ib2,1)=rho(Ib1:Ib2,1)*ustar**2*Twplus/mu(Ib1:Ib2,1)
 !$OMP END WORKSHARE
else if(isKe) then
 if(isGenlaw) then
  !$OMP WORKSHARE
  Te(Ib1:Ib2,1)=ustar**3*Uplus/Yp
  !$OMP END WORKSHARE
 else if(isLoglaw) then
  !$OMP WORKSHARE
  Te(Ib1:Ib2,1)=ustar**3/(kapa*Yp)
  !$OMP END WORKSHARE
 end if
end if
!$OMP END PARALLEL
end Subroutine Wallfunc