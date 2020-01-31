Subroutine Defercorrect(F,cor,Fw,Fe,Fs,Fn)
use Aero2DCOM
implicit none
integer i,j
real(8) Fcc,Fcw,Fce,Fcs,Fcn,Fcww,Fcee,Fcss,Fcnn,dkc,dkw,dkww,dke,dkee,dac,das,dass,dan,dann
real(8) rp,rm,Psip,Psim
real(8) F(Ic,Jc),cor(Ic,Jc),Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc)
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
  Fcs=F(i,j)
  Fcss=F(i,j)
  das=da(i,j)
  dass=da(i,j)
 else if(j==1) then
  Fcs=F(Ic+1-i,j)
  Fcss=F(Ic+1-i,j+1)
  das=da(Ic+1-i,j)
  dass=da(Ic+1-i,j+1)
 else if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
  Fcs=F(i,j-1)
  Fcss=F(i,j-1)
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
  cor(i,j)=dkw*(Fcw-Fcww)*max(Fw(i,j),0.0)/(dkw+dkww)+dkc*(Fce-Fcc)*max(-Fw(i,j),0.0)/(dkc+dke)+dkc*(Fcw-Fcc)*max(Fe(i,j),0.0)/(dkc+dkw)+&
  dke*(Fce-Fcee)*max(-Fe(i,j),0.0)/(dke+dkee)+das*(Fcs-Fcss)*max(Fs(i,j),0.0)/(das+dass)+dac*(Fcn-Fcc)*max(-Fs(i,j),0.0)/(dac+dan)+&
  dac*(Fcs-Fcc)*max(Fn(i,j),0.0)/(dac+das)+dan*(Fcn-Fcnn)*max(-Fn(i,j),0.0)/(dan+dann)
 else if(Discret=='Quick') then
  cor(i,j)=dkw*(3*(Fcc-Fcw)/(dkc+dkw)+(Fcw-Fcww)/(dkw+dkww))*max(Fw(i,j),0.0)/4+dkc*(3*(Fcc-Fce)/(dkc+dke)+(Fcw-Fcc)/(dkc+dkw))*max(Fe(i,j),0.0)/4+&
  dkc*(3*(Fcw-Fcc)/(dkc+dkw)+(Fcc-Fce)/(dkc+dke))*max(-Fw(i,j),0.0)/4+dke*(3*(Fce-Fcc)/(dkc+dke)+(Fcee-Fce)/(dke+dkee))*max(-Fe(i,j),0.0)/4+&
  das*(3*(Fcc-Fcs)/(dac+das)+(Fcs-Fcss)/(das+dass))*max(Fs(i,j),0.0)/4+dac*(3*(Fcc-Fcn)/(dac+dan)+(Fcs-Fcc)/(dac+das))*max(Fn(i,j),0.0)/4+&
  dac*(3*(Fcs-Fcc)/(dac+das)+(Fcc-Fcn)/(dac+dan))*max(-Fs(i,j),0.0)/4+dan*(3*(Fcn-Fcc)/(dac+dan)+(Fcnn-Fcn)/(dan+dann))*max(-Fn(i,j),0.0)/4
 else if(Discret=='tvd') then
  if(abs(Fcc-Fcw)<1e-30) then
   rp=0
   rm=0
  else
   rp=(Fcw-Fcww)/(Fcc-Fcw)
   rm=(Fce-Fcc)/(Fcc-Fcw)
  end if
  Psip=(rp+rp**2)/(1+rp**2)
  Psim=(rm+rm**2)/(1+rm**2)
  cor(i,j)=(max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim)*(Fcc-Fcw)/2
  if(abs(Fce-Fcc)<1e-30) then
   rp=0
   rm=0
  else
   rp=(Fcc-Fcw)/(Fce-Fcc)
   rm=(Fcee-Fce)/(Fce-Fcc)
  end if
  Psip=(rp+rp**2)/(1+rp**2)
  Psim=(rm+rm**2)/(1+rm**2)
  cor(i,j)=cor(i,j)+(max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim)*(Fcc-Fce)/2
  if(abs(Fcc-Fcs)<1e-30) then
   rp=0
   rm=0
  else
   rp=(Fcs-Fcss)/(Fcc-Fcs)
   rm=(Fcn-Fcc)/(Fcc-Fcs)
  end if
  Psip=(rp+rp**2)/(1+rp**2)
  Psim=(rm+rm**2)/(1+rm**2)
  cor(i,j)=cor(i,j)+(max(Fs(i,j),0.0)*Psip+max(-Fs(i,j),0.0)*Psim)*(Fcc-Fcs)/2
  if(abs(Fcn-Fcc)<1e-30) then
   rp=0
   rm=0
  else
   rp=(Fcc-Fcs)/(Fcn-Fcc)
   rm=(Fcnn-Fcn)/(Fcn-Fcc)
  end if
  Psip=(rp+rp**2)/(1+rp**2)
  Psim=(rm+rm**2)/(1+rm**2)
  cor(i,j)=cor(i,j)+(max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim)*(Fcc-Fcn)/2
 end if
 end DO
end DO
end Subroutine Defercorrect
