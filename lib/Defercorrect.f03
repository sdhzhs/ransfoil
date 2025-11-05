Subroutine Defercorrect(F,Fwall,cor,Fw,Fe,Fs,Fn)
use Aero2DCOM
implicit none
integer i,j
real(8) Fcc,Fcw,Fce,Fcs,Fcn,Fcww,Fcee,Fcss,Fcnn,dkcc,dkcw,dkce,dkcee,dacc,dacs,dacn,dacnn
real(8) rp,rm,Psip,Psim
real(8) F(Ic,Jc),Fwall(Ib1:Ib2),cor(Ic,Jc),Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc)
logical(1) isLr,isUp,is2ndUp,isQuick,isTvd

isLr = Walltreat=='lr'
isUp = Discret=='1upwind'
is2ndUp = Discret=='2upwind'
isQuick = Discret=='Quick'
isTvd = Discret=='tvd'

!$OMP DO
DO j=1,Jc-1
 DO i=Is,Ie
 Fcc=F(i,j)
 dkcc=dkd(i,j)
 dkce=dkd(i+1,j)
 if(i==1) then
  Fcw=F(Ic,j)
  dkcw=dkd(Ic,j)
 else
  Fcw=F(i-1,j)
  dkcw=dkd(i-1,j)
 end if
 if(i==Ic) then
  Fce=F(1,j)
  dkcee=dkd(2,j)
 else
  Fce=F(i+1,j)
  dkcee=dkd(i+2,j)
 end if
 if(i==1) then
  Fcww=F(Ic-1,j)
 else if(i==2) then
  if(Is>1) then
   Fcww=F(i-1,j)
  else
   Fcww=F(Ic,j)
  end if
 else
  Fcww=F(i-2,j)
 end if
 if(i==Ic) then
  Fcee=F(2,j)
 else if(i==Ic-1) then
  if(Is>1) then
   Fcee=F(i+1,j)
  else
   Fcee=F(1,j)
  end if
 else
  Fcee=F(i+2,j)
 end if
 Fcn=F(i,j+1)
 dacc=dad(i,j)
 dacn=dad(i,j+1)
 dacnn=dad(i,j+2)
 if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
  if(isLr) then
   if(isQuick) then
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
  dacs=dad(i,j)
 else if(j==1) then
  Fcs=F(Ic+1-i,j)
  Fcss=F(Ic+1-i,j+1)
  dacs=dad(Ic+1-i,j+1)
 else if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
  Fcs=F(i,j-1)
  if(isLr) then
   if(isQuick) then
    Fcss=(F(i,j)+10*F(i,j-1)-8*Fwall(i))/3
   else
    Fcss=2*Fwall(i)-F(i,j-1)
   end if
  else
   Fcss=F(i,j-1)
  end if
  dacs=dad(i,j-1)
 else if(j==2) then
  Fcs=F(i,j-1)
  Fcss=F(Ic+1-i,j-1)
  dacs=dad(i,j-1)
 else
  Fcs=F(i,j-1)
  Fcss=F(i,j-2)
  dacs=dad(i,j-1)
 end if
 if(j==Jc-1) then
  Fcnn=F(i,j+1)
 else
  Fcnn=F(i,j+2)
 end if
 if(isUp) then
  cor(i,j)=0
 else if(is2ndUp) then
  Psip=dkcc*(Fcw-Fcww)/(dkcc+dkcw)
  Psim=dkcc*(Fce-Fcc)/(dkcc+dkce)
  cor(i,j)=max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim
  Psip=(Fcw-Fcc)/2
  Psim=dkce*(Fce-Fcee)/(dkce+dkcee)
  cor(i,j)=cor(i,j)+max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim
  Psip=dacc*(Fcs-Fcss)/(dacc+dacs)
  Psim=dacc*(Fcn-Fcc)/(dacc+dacn)
  cor(i,j)=cor(i,j)+max(Fs(i,j),0.0)*Psip+max(-Fs(i,j),0.0)*Psim
  Psip=(Fcs-Fcc)/2
  Psim=dacn*(Fcn-Fcnn)/(dacn+dacnn)
  cor(i,j)=cor(i,j)+max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim
 else if(isQuick) then
  if(i==2) then
   Psip=0
  else
   Psip=dkcc*(3*(Fcc-Fcw)/(2*dkcc)+(Fcw-Fcww)/(dkcc+dkcw))/4
   !Psip=(3*Fcc-2*Fcw-Fcww)/8
  end if
  Psim=dkcc*(3*(Fcw-Fcc)/(2*dkcc)+(Fcc-Fce)/(dkcc+dkce))/4
  !Psim=(3*Fcw-2*Fcc-Fce)/8
  cor(i,j)=max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim
  if(i==Ic-1) then
   Psim=0
  else
   Psim=dkce*(3*(Fce-Fcc)/(dkcc+dkce)+(Fcee-Fce)/(dkce+dkcee))/4
   !Psim=(2*Fce+Fcee-3*Fcc)/8
  end if
  Psip=dkcc*(3*(Fcc-Fce)/(dkcc+dkce)+(Fcw-Fcc)/(2*dkcc))/4
  !Psip=(2*Fcc+Fcw-3*Fce)/8
  cor(i,j)=cor(i,j)+max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim
  if(j==Jc-1) then
   Psim=0
  else
   Psim=dacn*(3*(Fcn-Fcc)/(dacc+dacn)+(Fcnn-Fcn)/(dacn+dacnn))/4
   !Psim=(2*Fcn+Fcnn-3*Fcc)/8
  end if
  Psip=dacc*(3*(Fcc-Fcn)/(dacc+dacn)+(Fcs-Fcc)/(2*dacc))/4
  !Psip=(2*Fcc+Fcs-3*Fcn)/8
  cor(i,j)=cor(i,j)+max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim
  Psip=dacc*(3*(Fcc-Fcs)/(2*dacc)+(Fcs-Fcss)/(dacc+dacs))/4
  !Psip=(3*Fcc-2*Fcs-Fcss)/8
  Psim=dacc*(3*(Fcs-Fcc)/(2*dacc)+(Fcc-Fcn)/(dacc+dacn))/4
  !Psim=(3*Fcs-2*Fcc-Fcn)/8
  cor(i,j)=cor(i,j)+max(Fs(i,j),0.0)*Psip+max(-Fs(i,j),0.0)*Psim
 else if(isTvd) then
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