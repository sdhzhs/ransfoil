Subroutine Densityinterp
use Aero2DCOM
implicit none
integer i,j
real(8) rhoc,rhow,rhoe,rhos,rhon,rhoww,rhoss,dkc,dkw,dke,dkww,dac,das,dan,dass
real(8) wk,wa,rwp,rwm,rsp,rsm,Psiwp,Psiwm,Psisp,Psism
real(8),external:: interpl
  DO j=1,Jc-1
   DO i=2,Ic
   wk=sign(0.5,Unk(i,j))
   rhoc=rho(i,j)
   rhow=rho(i-1,j)
   dkc=dk(i,j)
   dkw=dk(i-1,j)
   if(i==Ic) then
    rhoe=rho(i,j)
    dke=dk(i,j)
   else
    rhoe=rho(i+1,j)
    dke=dk(i+1,j)
   end if
   if(i==2) then
    rhoww=rho(i-1,j)
    dkww=dk(i-1,j)
   else
    rhoww=rho(i-2,j)
    dkww=dk(i-2,j)
   end if
   if(Proctrl=='incom') then
   rhok(i,j)=rhoc
   else if(Proctrl=='com') then
   if(denface=='center') then
   rhok(i,j)=interpl(rhow,rhoc,dkw,dkc)
   else if(denface=='1upwind') then
   rhok(i,j)=(0.5-wk)*rhoc+(0.5+wk)*rhow
   else if(denface=='2upwind') then
   rhok(i,j)=(0.5-wk)*((2*dkc+dke)*rhoc-dkc*rhoe)/(dkc+dke)+(0.5+wk)*((2*dkw+dkww)*rhow-dkw*rhoww)/(dkw+dkww)
   else if(denface=='Quick') then
   rhok(i,j)=(0.5-wk)*(3*(dkw*rhoc+dkc*rhow)/(dkc+dkw)+((2*dkc+dke)*rhoc-dkc*rhoe)/(dkc+dke))/4+&
   (0.5+wk)*(3*(dkw*rhoc+dkc*rhow)/(dkc+dkw)+((2*dkw+dkww)*rhow-dkw*rhoww)/(dkw+dkww))/4
   else if(denface=='tvd') then
   if(abs(rhoc-rhow)<1e-30) then
   rwp=0
   else
   rwp=(rhow-rhoww)/(rhoc-rhow)
   end if
   if(abs(rhoc-rhow)<1e-30) then
   rwm=0
   else
   rwm=(rhoe-rhoc)/(rhoc-rhow)
   end if
   Psiwp=(rwp+rwp**2)/(1+rwp**2)
   Psiwm=(rwm+rwm**2)/(1+rwm**2)
   rhok(i,j)=(0.5-wk)*(rhoc+Psiwm*(rhow-rhoc)/2)+(0.5+wk)*(rhow+Psiwp*(rhoc-rhow)/2)
   end if
   end if
   end DO
  end DO
  DO j=1,Jc
   DO i=2,Ic-1
   wa=sign(0.5,Vna(i,j))
   rhoc=rho(i,j)
   dac=da(i,j)
   if(j==Jc) then
    rhon=rho(i,j)
    dan=da(i,j)
   else
    rhon=rho(i,j+1)
    dan=da(i,j+1)
   end if
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    rhos=rho(i,j)
    das=da(i,j)
   else if(j==1) then
    rhos=rho(Ic+1-i,j)
    das=da(Ic+1-i,j)
   else
    rhos=rho(i,j-1)
    das=da(i,j-1)
   end if
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    rhoss=rho(i,j)
    dass=da(i,j)
   else if(j==1) then
    rhoss=rho(Ic+1-i,j+1)
    dass=da(Ic+1-i,j+1)
   else if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
    rhoss=rho(i,j-1)
    dass=da(i,j-1)
   else if(j==2) then
    rhoss=rho(Ic+1-i,j-1)
    dass=da(Ic+1-i,j-1)
   else
    rhoss=rho(i,j-2)
    dass=da(i,j-2)
   end if
   if(Proctrl=='incom') then
   rhoa(i,j)=rhoc
   else if(Proctrl=='com') then
   if(denface=='center') then
   rhoa(i,j)=interpl(rhos,rhoc,das,dac)
   else if(denface=='1upwind') then
   rhoa(i,j)=(0.5-wa)*rhoc+(0.5+wa)*rhos
   else if(denface=='2upwind') then
   rhoa(i,j)=(0.5-wa)*((2*dac+dan)*rhoc-dac*rhon)/(dac+dan)+(0.5+wa)*((2*das+dass)*rhos-das*rhoss)/(das+dass)
   else if(denface=='Quick') then
   rhoa(i,j)=(0.5-wa)*(3*(das*rhoc+dac*rhos)/(dac+das)+((2*dac+dan)*rhoc-dac*rhon)/(dac+dan))/4+&
   (0.5+wa)*(3*(das*rhoc+dac*rhos)/(dac+das)+((2*das+dass)*rhos-das*rhoss)/(das+dass))/4
   else if(denface=='tvd') then
   if(abs(rhoc-rhos)<1e-30) then
   rsp=0
   else
   rsp=(rhos-rhoss)/(rhoc-rhos)
   end if
   if(abs(rhoc-rhos)<1e-30) then
   rsm=0
   else
   rsm=(rhon-rhoc)/(rhoc-rhos)
   end if
   Psisp=(rsp+rsp**2)/(1+rsp**2)
   Psism=(rsm+rsm**2)/(1+rsm**2)
   rhoa(i,j)=(0.5-wa)*(rhoc+Psism*(rhos-rhoc)/2)+(0.5+wa)*(rhos+Psisp*(rhoc-rhos)/2)
   end if
   end if
   end DO
  end DO
end Subroutine Densityinterp
