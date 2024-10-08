Subroutine Turvis
use Aero2DCOM
implicit none
integer i,j
real(8) Xi,fnu1,Dwplus,phi1,F1,betai,St,phi2,F2,Ret,alphastar,Ymax,Ym
logical(1) sstlowre
 sstlowre=.false.
 Ym=11.225
 if(Turmod=='ke') then
  Ymax=maxval(Ystar)
 else
  Ymax=maxval(Yplus)
 end if
 DO j=1,Jc
  DO i=1,Ic
   if(Turmod=='sa') then
    if(Walltreat=='lr') then
     Xi=rho(i,j)*Tn(i,j)/mu(i,j)+Cks*ksi/d(i,j)
    else
     Xi=rho(i,j)*Tn(i,j)/mu(i,j)
    end if
    fnu1=Xi**3/(Xi**3+Cnu1**3)
    !if(Ymax>Ym) fnu1=1.0
    mut(i,j)=rho(i,j)*Tn(i,j)*fnu1
   else if(Turmod=='ke') then
    mut(i,j)=(rho(i,j)*Cu*Tk(i,j)**2)/Te(i,j)
   else if(Turmod=='sst') then
    Dwplus=max(2*rho(i,j)*(Tkx(i,j)*Twx(i,j)+Tky(i,j)*Twy(i,j))/(sigmaw2*Tw(i,j)),1e-10)
    phi1=min(max(sqrt(Tk(i,j))/(betastarf*Tw(i,j)*d(i,j)),500*mu(i,j)/(rho(i,j)*d(i,j)**2*Tw(i,j))),&
    4*rho(i,j)*Tk(i,j)/(sigmaw2*Dwplus*d(i,j)**2))
    F1=tanh(phi1**4)
    betai=F1*betai1+(1-F1)*betai2
    St=sqrt(2*(Ux(i,j)**2+Vy(i,j)**2)+(Uy(i,j)+Vx(i,j))**2)
    phi2=max(2*sqrt(Tk(i,j))/(betastarf*Tw(i,j)*d(i,j)),500*mu(i,j)/(rho(i,j)*d(i,j)**2*Tw(i,j)))
    F2=tanh(phi2**2)
    if(sstlowre) then
      Ret=rho(i,j)*Tk(i,j)/(mu(i,j)*Tw(i,j))
      alphastar=alphastarf*(0.024+Ret/Rk)/(1+Ret/Rk)
    else
      alphastar=alphastarf
    end if
    mut(i,j)=(rho(i,j)*Tk(i,j))/(Tw(i,j)*max(1./alphastar,St*F2/(alpha1*Tw(i,j))))
   else if(Turmod=='lam'.or.Turmod=='inv') then
    mut(i,j)=0
   end if
  end DO
 end DO
end Subroutine Turvis
