Subroutine Densityinterp
use Aero2DCOM
implicit none
integer i,j
real(8) rhoc,rhow,rhoe,rhos,rhon,rhoww,rhoss,dkc,dkcc,dkcw,dkce,dkcww,dac,dacc,dacs,dacn,dacss
real(8) wk,wa,rwp,rwm,rsp,rsm,Psiwp,Psiwm,Psisp,Psism,Xfc,Yfc,Xc0,Yc0,Xc1,Yc1,Gxupc0,Gxupc1,Gyupc0,Gyupc1
real(8),external:: interpl
logical(1) isCom,isIncom,isCenter,isUp,is2ndUp,isQuick,isTvd,isFromm

isCom = Proctrl=='com'
isIncom = Proctrl=='incom'
isCenter = denface=='center'
isUp = denface=='1upwind'
is2ndUp = denface=='2upwind'
isQuick = denface=='Quick'
isTvd = denface=='tvd'
isFromm = denface=='Fromm'

!$OMP PARALLEL
  !$OMP DO PRIVATE(wk,rhoc,rhow,rhoe,rhoww,dkc,dkcc,dkcw,dkce,dkcww,rwp,rwm,Psiwp,Psiwm,Xfc,Yfc,Xc0,Yc0,Xc1,Yc1,Gxupc0,Gxupc1,Gyupc0,Gyupc1,i)
  DO j=1,Jc-1
   DO i=Is,Ie+1
    wk=sign(0.5,Unk(i,j))
    dkc=dkw(i,j)
    if(i==Ic+1) then
     rhoc=rho(1,j)
     dkcc=dk(1,j)
     Gxupc1=rhox(1,j)
     Gyupc1=rhoy(1,j)
     Xc1=Xc(i,j)
     Yc1=Yc(i,j)
    else
     rhoc=rho(i,j)
     dkcc=dk(i,j)
     Gxupc1=rhox(i,j)
     Gyupc1=rhoy(i,j)
     Xc1=Xc(i,j)
     Yc1=Yc(i,j)
    end if
    if(i==1) then
     rhow=rho(Ic,j)
     dkcw=dk(Ic,j)
     Gxupc0=rhox(Ic,j)
     Gyupc0=rhoy(Ic,j)
     Xc0=Xc(Ic,j)
     Yc0=Yc(Ic,j)
    else
     rhow=rho(i-1,j)
     dkcw=dk(i-1,j)
     Gxupc0=rhox(i-1,j)
     Gyupc0=rhoy(i-1,j)
     Xc0=Xc(i-1,j)
     Yc0=Yc(i-1,j)
    end if
    if(i==Ic+1) then
     rhoe=rho(2,j)
     dkce=dk(2,j)
    else if(i==Ic) then
     if(Is>1) then
      rhoe=rho(i,j)
      dkce=dk(i,j)
     else
      rhoe=rho(1,j)
      dkce=dk(1,j)
     end if
    else
     rhoe=rho(i+1,j)
     dkce=dk(i+1,j)
    end if
    if(i==1) then
     rhoww=rho(Ic-1,j)
     dkcww=dk(Ic-1,j)
    else if(i==2) then
     if(Is>1) then
      rhoww=rho(i-1,j)
      dkcww=dk(i-1,j)
     else
      rhoww=rho(Ic,j)
      dkcww=dk(Ic,j)
     end if
    else
     rhoww=rho(i-2,j)
     dkcww=dk(i-2,j)
    end if
    if(isIncom) then
     rhok(i,j)=rhoc
    else if(isCom) then
     if(isCenter) then
      rhok(i,j)=interpl(rhow,rhoc,dkc)
     else if(isUp) then
      rhok(i,j)=(0.5-wk)*rhoc+(0.5+wk)*rhow
     else if(is2ndUp) then
      rhok(i,j)=(0.5-wk)*((2*dkcc+dkce)*rhoc-dkcc*rhoe)/(dkcc+dkce)+(0.5+wk)*((2*dkcw+dkcww)*rhow-dkcw*rhoww)/(dkcw+dkcww)
     else if(isFromm) then
      Xfc=0.5*(Xg(i,j)+Xg(i,j+1))
      Yfc=0.5*(Yg(i,j)+Yg(i,j+1))
      rhok(i,j)=(0.5-wk)*(rhoc+Gxupc1*(Xfc-Xc1)+Gyupc1*(Yfc-Yc1))+(0.5+wk)*(rhow+Gxupc0*(Xfc-Xc0)+Gyupc0*(Yfc-Yc0))
     else if(isQuick) then
      rhok(i,j)=(0.5-wk)*(3*(dkcw*rhoc+dkcc*rhow)/(dkcc+dkcw)+((2*dkcc+dkce)*rhoc-dkcc*rhoe)/(dkcc+dkce))/4+&
      (0.5+wk)*(3*(dkcw*rhoc+dkcc*rhow)/(dkcc+dkcw)+((2*dkcw+dkcww)*rhow-dkcw*rhoww)/(dkcw+dkcww))/4
     else if(isTvd) then
      if(abs(rhoc-rhow)>0) then
       rwp=(rhow-rhoww)/(rhoc-rhow)
      else
       rwp=0
      end if
      if(abs(rhoc-rhow)>0) then
       rwm=(rhoe-rhoc)/(rhoc-rhow)
      else
       rwm=0
      end if
      Psiwp=(rwp+rwp**2)/(1+rwp**2)
      Psiwm=(rwm+rwm**2)/(1+rwm**2)
      rhok(i,j)=(0.5-wk)*(rhoc+Psiwm*(rhow-rhoc)/2)+(0.5+wk)*(rhow+Psiwp*(rhoc-rhow)/2)
     end if
    end if
   end DO
  end DO
  !$OMP END DO
  !$OMP DO PRIVATE(wa,rhoc,rhos,rhon,rhoss,dac,dacc,dacs,dacn,dacss,rsp,rsm,Psisp,Psism,Xfc,Yfc,Xc0,Yc0,Xc1,Yc1,Gxupc0,Gxupc1,Gyupc0,Gyupc1,i)
  DO j=1,Jc
   DO i=Is,Ie
    wa=sign(0.5,Vna(i,j))
    dac=daw(i,j)
    rhoc=rho(i,j)
    dacc=da(i,j)
    Xc1=Xc(i,j)
    Yc1=Yc(i,j)
    Gxupc1=rhox(i,j)
    Gyupc1=rhoy(i,j)
    if(j==Jc) then
     rhon=rho(i,j)
     dacn=da(i,j)
    else
     rhon=rho(i,j+1)
     dacn=da(i,j+1)
    end if
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     rhos=rho(i,j)
     dacs=da(i,j)
     Xc0=Xc(i,j)
     Yc0=Yc(i,j)
     Gxupc0=0
     Gyupc0=0
    else if(j==1) then
     rhos=rho(Ic+1-i,j)
     dacs=da(Ic+1-i,j)
     Xc0=Xc(Ic+1-i,j)
     Yc0=Yc(Ic+1-i,j)
     Gxupc0=rhox(Ic+1-i,j)
     Gyupc0=rhoy(Ic+1-i,j)
    else
     rhos=rho(i,j-1)
     dacs=da(i,j-1)
     Xc0=Xc(i,j-1)
     Yc0=Yc(i,j-1)
     Gxupc0=rhox(i,j-1)
     Gyupc0=rhoy(i,j-1)
    end if
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     rhoss=rho(i,j)
     dacss=da(i,j)
    else if(j==1) then
     rhoss=rho(Ic+1-i,j+1)
     dacss=da(Ic+1-i,j+1)
    else if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
     rhoss=rho(i,j-1)
     dacss=da(i,j-1)
    else if(j==2) then
     rhoss=rho(Ic+1-i,j-1)
     dacss=da(Ic+1-i,j-1)
    else
     rhoss=rho(i,j-2)
     dacss=da(i,j-2)
    end if
    if(isIncom) then
     rhoa(i,j)=rhoc
    else if(isCom) then
     if(isCenter) then
      rhoa(i,j)=interpl(rhos,rhoc,dac)
     else if(isUp) then
      rhoa(i,j)=(0.5-wa)*rhoc+(0.5+wa)*rhos
     else if(is2ndUp) then
      rhoa(i,j)=(0.5-wa)*((2*dacc+dacn)*rhoc-dacc*rhon)/(dacc+dacn)+(0.5+wa)*((2*dacs+dacss)*rhos-dacs*rhoss)/(dacs+dacss)
     else if(isFromm) then
      Xfc=0.5*(Xg(i,j)+Xg(i+1,j))
      Yfc=0.5*(Yg(i,j)+Yg(i+1,j))
      rhoa(i,j)=(0.5-wa)*(rhoc+Gxupc1*(Xfc-Xc1)+Gyupc1*(Yfc-Yc1))+(0.5+wa)*(rhos+Gxupc0*(Xfc-Xc0)+Gyupc0*(Yfc-Yc0))
     else if(isQuick) then
      rhoa(i,j)=(0.5-wa)*(3*(dacs*rhoc+dacc*rhos)/(dacc+dacs)+((2*dacc+dacn)*rhoc-dacc*rhon)/(dacc+dacn))/4+&
      (0.5+wa)*(3*(dacs*rhoc+dacc*rhos)/(dacc+dacs)+((2*dacs+dacss)*rhos-dacs*rhoss)/(dacs+dacss))/4
     else if(isTvd) then
      if(abs(rhoc-rhos)>0) then
       rsp=(rhos-rhoss)/(rhoc-rhos)
      else
       rsp=0
      end if
      if(abs(rhoc-rhos)>0) then
       rsm=(rhon-rhoc)/(rhoc-rhos)
      else
       rsm=0
      end if
      Psisp=(rsp+rsp**2)/(1+rsp**2)
      Psism=(rsm+rsm**2)/(1+rsm**2)
      rhoa(i,j)=(0.5-wa)*(rhoc+Psism*(rhos-rhoc)/2)+(0.5+wa)*(rhos+Psisp*(rhoc-rhos)/2)
     end if
    end if
   end DO
  end DO
  !$OMP END DO
!$OMP END PARALLEL
end Subroutine Densityinterp