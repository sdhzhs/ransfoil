Subroutine Densityinterp
use Aero2DCOM
implicit none
integer i,j
real(8) rhoc,rhow,rhoe,rhos,rhon,rhoww,rhoss,dkc,dkcc,dkcw,dkce,dac,dacc,dacs,dacn
real(8) wk,wa,rwp,rwm,rsp,rsm,Psiwp,Psiwm,Psisp,Psism
real(8),external:: interpl
logical(1) isCom,isIncom,isCenter,isUp,is2ndUp,isQuick,isTvd

isCom = Proctrl=='com'
isIncom = Proctrl=='incom'
isCenter = denface=='center'
isUp = denface=='1upwind'
is2ndUp = denface=='2upwind'
isQuick = denface=='Quick'
isTvd = denface=='tvd'

!$OMP PARALLEL
  !$OMP DO PRIVATE(wk,rhoc,rhow,rhoe,rhoww,dkc,dkcc,dkcw,dkce,rwp,rwm,Psiwp,Psiwm,i)
  DO j=1,Jc-1
   DO i=Is,Ie+1
    wk=sign(0.5,Unk(i,j))
    dkc=dkw(i,j)
    dkcc=dkd(i,j)
    if(i==Ic+1) then
     rhoc=rho(1,j)
     dkce=dkd(2,j)
    else
     rhoc=rho(i,j)
     dkce=dkd(i+1,j)
    end if
    if(i==1) then
     rhow=rho(Ic,j)
     dkcw=dkd(Ic,j)
    else
     rhow=rho(i-1,j)
     dkcw=dkd(i-1,j)
    end if
    if(i==Ic+1) then
     rhoe=rho(2,j)
    else if(i==Ic) then
     if(Is>1) then
      rhoe=rho(i,j)
     else
      rhoe=rho(1,j)
     end if
    else
     rhoe=rho(i+1,j)
    end if
    if(i==1) then
     rhoww=rho(Ic-1,j)
    else if(i==2) then
     if(Is>1) then
      rhoww=rho(i-1,j)
     else
      rhoww=rho(Ic,j)
     end if
    else
     rhoww=rho(i-2,j)
    end if
    if(isIncom) then
     rhok(i,j)=rhoc
    else if(isCom) then
     if(isCenter) then
      rhok(i,j)=interpl(rhow,rhoc,dkc)
     else if(isUp) then
      rhok(i,j)=(0.5-wk)*rhoc+(0.5+wk)*rhow
     else if(is2ndUp) then
      rhok(i,j)=(0.5-wk)*((2*dkcc+dkce)*rhoc-dkcc*rhoe)/(dkcc+dkce)+(0.5+wk)*((2*dkcc+dkcw)*rhow-dkcc*rhoww)/(dkcw+dkcc)
     else if(isQuick) then
      rhok(i,j)=(0.5-wk)*(3*interpl(rhow,rhoc,dkc)+((2*dkcc+dkce)*rhoc-dkcc*rhoe)/(dkcc+dkce))/4+&
      (0.5+wk)*(3*interpl(rhow,rhoc,dkc)+((2*dkcc+dkcw)*rhow-dkcc*rhoww)/(dkcc+dkcw))/4
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
  !$OMP DO PRIVATE(wa,rhoc,rhos,rhon,rhoss,dac,dacc,dacs,dacn,rsp,rsm,Psisp,Psism,i)
  DO j=1,Jc
   DO i=Is,Ie
    wa=sign(0.5,Vna(i,j))
    rhoc=rho(i,j)
    dac=daw(i,j)
    dacc=dad(i,j)
    if(j==Jc) then
     rhon=rho(i,j)
     dacn=2*dad(i,j+1)
    else
     rhon=rho(i,j+1)
     dacn=dad(i,j+1)
    end if
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     rhos=rho(i,j)
     dacs=dad(i,j)
    else if(j==1) then
     rhos=rho(Ic+1-i,j)
     dacs=dad(Ic+1-i,j+1)
    else
     rhos=rho(i,j-1)
     dacs=dad(i,j-1)
    end if
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
     rhoss=rho(i,j)
    else if(j==1) then
     rhoss=rho(Ic+1-i,j+1)
    else if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
     rhoss=rho(i,j-1)
    else if(j==2) then
     rhoss=rho(Ic+1-i,j-1)
    else
     rhoss=rho(i,j-2)
    end if
    if(isIncom) then
     rhoa(i,j)=rhoc
    else if(isCom) then
     if(isCenter) then
      rhoa(i,j)=interpl(rhos,rhoc,dac)
     else if(isUp) then
      rhoa(i,j)=(0.5-wa)*rhoc+(0.5+wa)*rhos
     else if(is2ndUp) then
      rhoa(i,j)=(0.5-wa)*((2*dacc+dacn)*rhoc-dacc*rhon)/(dacc+dacn)+(0.5+wa)*((2*dacc+dacs)*rhos-dacc*rhoss)/(dacc+dacs)
     else if(isQuick) then
      rhoa(i,j)=(0.5-wa)*(3*interpl(rhos,rhoc,dac)+((2*dacc+dacn)*rhoc-dacc*rhon)/(dacc+dacn))/4+&
      (0.5+wa)*(3*interpl(rhos,rhoc,dac)+((2*dacc+dacs)*rhos-dacc*rhoss)/(dacc+dacs))/4
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