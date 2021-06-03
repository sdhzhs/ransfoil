Subroutine Defercorrect(F,Fwall,cor,Fw,Fe,Fs,Fn)
use Aero2DCOM
implicit none
integer i,j
real(8) Fcc,Fcw,Fce,Fcs,Fcn,Fcww,Fcee,Fcss,Fcnn,dkc,dkw,dkww,dke,dkee,dac,das,dass,dan,dann
real(8) rp,rm,Psip,Psim
real(8) F(Ic,Jc),Fwall(Ib1:Ib2),cor(Ic,Jc),Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc)
!$OMP DO
DO j=1,Jc-1
 DO i=2,Ic-1
 Fcc=F(i,j)
 Fcw=F(i-1,j)
 Fce=F(i+1,j)
 Fcn=F(i,j+1)
 dkc=dk(i,j)
 dac=da(i,j)
 dkw=dk(i-1,j)
 dke=dk(i+1,j)
 dan=da(i,j+1)
 if(i==2) then
  Fcww=F(i-1,j)
  dkww=dk(i-1,j)
 else
  Fcww=F(i-2,j)
  dkww=dk(i-2,j)
 end if
 if(i==Ic-1) then
  Fcee=F(i+1,j)
  dkee=dk(i+1,j)
 else
  Fcee=F(i+2,j)
  dkee=dk(i+2,j)
 end if
 if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
  if(Walltreat=='lr') then
   if(Discret=='Quick') then
    Fcs=(F(i,j+1)+10*F(i,j)-8*Fwall(i))/3
    Fcss=6*Fcs+3*F(i,j)-8*Fwall(i)
   else
    Fcs=2*Fwall(i)-F(i,j)
    Fcss=2*Fwall(i)-F(i,j+1)
   end if
  else
   Fcs=F(i,j)
   Fcss=F(i,j)
  end if
  das=da(i,j)
  dass=da(i,j)
 else if(j==1) then
  Fcs=F(Ic+1-i,j)
  Fcss=F(Ic+1-i,j+1)
  das=da(Ic+1-i,j)
  dass=da(Ic+1-i,j+1)
 else if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
  Fcs=F(i,j-1)
  if(Walltreat=='lr') then
   if(Discret=='Quick') then
    Fcss=(F(i,j)+10*F(i,j-1)-8*Fwall(i))/3
   else
    Fcss=2*Fwall(i)-F(i,j-1)
   end if
  else
   Fcss=F(i,j-1)
  end if
  das=da(i,j-1)
  dass=da(i,j-1)
 else if(j==2) then
  Fcs=F(i,j-1)
  Fcss=F(Ic+1-i,j-1)
  das=da(i,j-1)
  dass=da(Ic+1-i,j-1)
 else
  Fcs=F(i,j-1)
  Fcss=F(i,j-2)
  das=da(i,j-1)
  dass=da(i,j-2)
 end if
 if(j==Jc-1) then
  Fcnn=F(i,j+1)
  dann=da(i,j+1)
 else
  Fcnn=F(i,j+2)
  dann=da(i,j+2)
 end if
 if(Discret=='1upwind') then
  cor(i,j)=0
 else if(Discret=='2upwind') then
  Psip=dkw*(Fcw-Fcww)/(dkw+dkww)
  Psim=dkc*(Fce-Fcc)/(dkc+dke)
  cor(i,j)=max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim
  Psip=dkc*(Fcw-Fcc)/(dkc+dkw)
  Psim=dke*(Fce-Fcee)/(dke+dkee)
  cor(i,j)=cor(i,j)+max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim
  Psip=das*(Fcs-Fcss)/(das+dass)
  Psim=dac*(Fcn-Fcc)/(dac+dan)
  cor(i,j)=cor(i,j)+max(Fs(i,j),0.0)*Psip+max(-Fs(i,j),0.0)*Psim
  Psip=dac*(Fcs-Fcc)/(dac+das)
  Psim=dan*(Fcn-Fcnn)/(dan+dann)
  cor(i,j)=cor(i,j)+max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim
 else if(Discret=='Quick') then
  if(i==2) then
   Psip=0
  else
   Psip=dkw*(3*(Fcc-Fcw)/(dkc+dkw)+(Fcw-Fcww)/(dkw+dkww))/4
   !Psip=(3*Fcc-2*Fcw-Fcww)/8
  end if
  Psim=dkc*(3*(Fcw-Fcc)/(dkc+dkw)+(Fcc-Fce)/(dkc+dke))/4
  !Psim=(3*Fcw-2*Fcc-Fce)/8
  cor(i,j)=max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim
  if(i==Ic-1) then
   Psim=0
  else
   Psim=dke*(3*(Fce-Fcc)/(dkc+dke)+(Fcee-Fce)/(dke+dkee))/4
   !Psim=(2*Fce+Fcee-3*Fcc)/8
  end if
  Psip=dkc*(3*(Fcc-Fce)/(dkc+dke)+(Fcw-Fcc)/(dkc+dkw))/4
  !Psip=(2*Fcc+Fcw-3*Fce)/8
  cor(i,j)=cor(i,j)+max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim
  if(j==Jc-1) then
   Psim=0
  else
   Psim=dan*(3*(Fcn-Fcc)/(dac+dan)+(Fcnn-Fcn)/(dan+dann))/4
   !Psim=(2*Fcn+Fcnn-3*Fcc)/8
  end if
  Psip=dac*(3*(Fcc-Fcn)/(dac+dan)+(Fcs-Fcc)/(dac+das))/4
  !Psip=(2*Fcc+Fcs-3*Fcn)/8
  cor(i,j)=cor(i,j)+max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim
  Psip=das*(3*(Fcc-Fcs)/(dac+das)+(Fcs-Fcss)/(das+dass))/4
  !Psip=(3*Fcc-2*Fcs-Fcss)/8
  Psim=dac*(3*(Fcs-Fcc)/(dac+das)+(Fcc-Fcn)/(dac+dan))/4
  !Psim=(3*Fcs-2*Fcc-Fcn)/8
  cor(i,j)=cor(i,j)+max(Fs(i,j),0.0)*Psip+max(-Fs(i,j),0.0)*Psim
 else if(Discret=='tvd') then
  if(abs(Fcc-Fcw)>0) then
   rp=(Fcw-Fcww)/(Fcc-Fcw)
   rm=(Fce-Fcc)/(Fcc-Fcw)
  else
   rp=0
   rm=0
  end if
  Psip=(rp+rp**2)/(1+rp**2)
  Psim=(rm+rm**2)/(1+rm**2)
  cor(i,j)=(max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim)*(Fcc-Fcw)/2
  if(abs(Fce-Fcc)>0) then
   rp=(Fcc-Fcw)/(Fce-Fcc)
   rm=(Fcee-Fce)/(Fce-Fcc)
  else
   rp=0
   rm=0
  end if
  Psip=(rp+rp**2)/(1+rp**2)
  Psim=(rm+rm**2)/(1+rm**2)
  cor(i,j)=cor(i,j)+(max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim)*(Fcc-Fce)/2
  if(abs(Fcc-Fcs)>0) then
   rp=(Fcs-Fcss)/(Fcc-Fcs)
   rm=(Fcn-Fcc)/(Fcc-Fcs)
  else
   rp=0
   rm=0
  end if
  Psip=(rp+rp**2)/(1+rp**2)
  Psim=(rm+rm**2)/(1+rm**2)
  cor(i,j)=cor(i,j)+(max(Fs(i,j),0.0)*Psip+max(-Fs(i,j),0.0)*Psim)*(Fcc-Fcs)/2
  if(abs(Fcn-Fcc)>0) then
   rp=(Fcc-Fcs)/(Fcn-Fcc)
   rm=(Fcnn-Fcn)/(Fcn-Fcc)
  else
   rp=0
   rm=0
  end if
  Psip=(rp+rp**2)/(1+rp**2)
  Psim=(rm+rm**2)/(1+rm**2)
  cor(i,j)=cor(i,j)+(max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim)*(Fcc-Fcn)/2
 end if
 end DO
end DO
!$OMP END DO
end Subroutine Defercorrect