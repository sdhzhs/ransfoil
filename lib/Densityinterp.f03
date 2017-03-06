Subroutine Densityinterp
use Aero2DCOM
implicit none
integer i,j
real(8) rouc,rouw,roue,rous,roun,rouww,rouss,dkc,dkw,dke,dkww,dac,das,dan,dass
real(8) wk,wa,rwp,rwm,rsp,rsm,Psiwp,Psiwm,Psisp,Psism

  DO j=1,Jc-1
   DO i=2,Ic
   wk=sign(0.5,Unk(i,j))
   rouc=rou(i,j)
   rouw=rou(i-1,j)
   dkc=dk(i,j)
   dkw=dk(i-1,j)
   if(i==Ic) then
    roue=rou(i,j)
    dke=dk(i,j)
   else
    roue=rou(i+1,j)
    dke=dk(i+1,j)
   end if
   if(i==2) then
    rouww=rou(i-1,j)
    dkww=dk(i-1,j)
   else
    rouww=rou(i-2,j)
    dkww=dk(i-2,j)
   end if
   if(Proctrl=='incom') then
   rouk(i,j)=rouc
   else if(Proctrl=='com') then
   if(denface=='center') then
   rouk(i,j)=interpl(rouw,rouc,dkw,dkc)
   else if(denface=='1upwind') then
   rouk(i,j)=(0.5-wk)*rouc+(0.5+wk)*rouw
   else if(denface=='2upwind') then
   rouk(i,j)=(0.5-wk)*((2*dkc+dke)*rouc-dkc*roue)/(dkc+dke)+(0.5+wk)*((2*dkw+dkww)*rouw-dkw*rouww)/(dkw+dkww)
   else if(denface=='Quick') then
   rouk(i,j)=(0.5-wk)*(3*(dkw*rouc+dkc*rouw)/(dkc+dkw)+((2*dkc+dke)*rouc-dkc*roue)/(dkc+dke))/4+&
   (0.5+wk)*(3*(dkw*rouc+dkc*rouw)/(dkc+dkw)+((2*dkw+dkww)*rouw-dkw*rouww)/(dkw+dkww))/4
   else if(denface=='tvd') then
   if(abs(rouc-rouw)<1e-30) then
   rwp=0
   else
   rwp=(rouw-rouww)/(rouc-rouw)
   end if
   if(abs(rouc-rouw)<1e-30) then
   rwm=0
   else
   rwm=(roue-rouc)/(rouc-rouw)
   end if
   Psiwp=(rwp+rwp**2)/(1+rwp**2)
   Psiwm=(rwm+rwm**2)/(1+rwm**2)
   rouk(i,j)=(0.5-wk)*(rouc+Psiwm*(rouw-rouc)/2)+(0.5+wk)*(rouw+Psiwp*(rouc-rouw)/2)
   end if
   end if
   end DO
  end DO
  DO j=1,Jc
   DO i=2,Ic-1
   wa=sign(0.5,Vna(i,j))
   rouc=rou(i,j)
   dac=da(i,j)
   if(j==Jc) then
    roun=rou(i,j)
    dan=da(i,j)
   else
    roun=rou(i,j+1)
    dan=da(i,j+1)
   end if
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    rous=rou(i,j)
    das=da(i,j)
   else if(j==1) then
    rous=rou(Ic+1-i,j)
    das=da(Ic+1-i,j)
   else
    rous=rou(i,j-1)
    das=da(i,j-1)
   end if
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    rouss=rou(i,j)
    dass=da(i,j)
   else if(j==1) then
    rouss=rou(Ic+1-i,j+1)
    dass=da(Ic+1-i,j+1)
   else if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
    rouss=rou(i,j-1)
    dass=da(i,j-1)
   else if(j==2) then
    rouss=rou(Ic+1-i,j-1)
    dass=da(Ic+1-i,j-1)
   else
    rouss=rou(i,j-2)
    dass=da(i,j-2)
   end if
   if(Proctrl=='incom') then
   roua(i,j)=rouc
   else if(Proctrl=='com') then
   if(denface=='center') then
   roua(i,j)=interpl(rous,rouc,das,dac)
   else if(denface=='1upwind') then
   roua(i,j)=(0.5-wa)*rouc+(0.5+wa)*rous
   else if(denface=='2upwind') then
   roua(i,j)=(0.5-wa)*((2*dac+dan)*rouc-dac*roun)/(dac+dan)+(0.5+wa)*((2*das+dass)*rous-das*rouss)/(das+dass)
   else if(denface=='Quick') then
   roua(i,j)=(0.5-wa)*(3*(das*rouc+dac*rous)/(dac+das)+((2*dac+dan)*rouc-dac*roun)/(dac+dan))/4+&
   (0.5+wa)*(3*(das*rouc+dac*rous)/(dac+das)+((2*das+dass)*rous-das*rouss)/(das+dass))/4
   else if(denface=='tvd') then
   if(abs(rouc-rous)<1e-30) then
   rsp=0
   else
   rsp=(rous-rouss)/(rouc-rous)
   end if
   if(abs(rouc-rous)<1e-30) then
   rsm=0
   else
   rsm=(roun-rouc)/(rouc-rous)
   end if
   Psisp=(rsp+rsp**2)/(1+rsp**2)
   Psism=(rsm+rsm**2)/(1+rsm**2)
   roua(i,j)=(0.5-wa)*(rouc+Psism*(rous-rouc)/2)+(0.5+wa)*(rous+Psisp*(rouc-rous)/2)
   end if
   end if
   end DO
  end DO
end Subroutine Densityinterp
