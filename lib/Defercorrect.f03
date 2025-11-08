Subroutine Defercorrect(F,Fx,Fy,Fwall,cor,Fw,Fe,Fs,Fn)
use Aero2DCOM
implicit none
integer i,j
real(8) Fcc,Fcw,Fce,Fcs,Fcn,Fcww,Fcee,Fcss,Fcnn,dkcc,dkcw,dkcww,dkce,dkcee,dacc,dacs,dacss,dacn,dacnn
real(8) rp,rm,Psip,Psim,Xfc,Yfc,Xcc(5),Ycc(5),Gxup(5),Gyup(5)
real(8) F(Ic,Jc),Fx(Ic,Jc),Fy(Ic,Jc),Fwall(Ib1:Ib2),cor(Ic,Jc),Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc)
logical(1) isLr,isUp,is2ndUp,isQuick,isTvd,isFromm

isLr = Walltreat=='lr'
isUp = Discret=='1upwind'
is2ndUp = Discret=='2upwind'
isQuick = Discret=='Quick'
isTvd = Discret=='tvd'
isFromm = Discret=='Fromm'

!$OMP DO
DO j=1,Jc-1
 DO i=Is,Ie
 Fcc=F(i,j)
 dkcc=dk(i,j)
 dacc=da(i,j)
 Xcc(1)=Xc(i,j)
 Ycc(1)=Yc(i,j)
 Gxup(1)=Fx(i,j)
 Gyup(1)=Fy(i,j)
 if(i==1) then
  Fcw=F(Ic,j)
  dkcw=dk(Ic,j)
  Xcc(2)=Xc(Ic,j)
  Ycc(2)=Yc(Ic,j)
  Gxup(2)=Fx(Ic,j)
  Gyup(2)=Fy(Ic,j)
 else
  Fcw=F(i-1,j)
  dkcw=dk(i-1,j)
  Xcc(2)=Xc(i-1,j)
  Ycc(2)=Yc(i-1,j)
  Gxup(2)=Fx(i-1,j)
  Gyup(2)=Fy(i-1,j)
 end if
 if(i==Ic) then
  Fce=F(1,j)
  dkce=dk(1,j)
  Xcc(3)=Xc(1,j)
  Ycc(3)=Yc(1,j)
  Gxup(3)=Fx(1,j)
  Gyup(3)=Fy(1,j)
 else
  Fce=F(i+1,j)
  dkce=dk(i+1,j)
  Xcc(3)=Xc(i+1,j)
  Ycc(3)=Yc(i+1,j)
  Gxup(3)=Fx(i+1,j)
  Gyup(3)=Fy(i+1,j)
 end if
 Fcn=F(i,j+1)
 dacn=da(i,j+1)
 Xcc(5)=Xc(i,j+1)
 Ycc(5)=Yc(i,j+1)
 Gxup(5)=Fx(i,j+1)
 Gyup(5)=Fy(i,j+1)
 if(i==1) then
  Fcww=F(Ic-1,j)
  dkcww=dk(Ic-1,j)
 else if(i==2) then
  if(Is>1) then
   Fcww=F(i-1,j)
   dkcww=dk(i-1,j)
  else
   Fcww=F(Ic,j)
   dkcww=dk(Ic,j)
  end if
 else
  Fcww=F(i-2,j)
  dkcww=dk(i-2,j)
 end if
 if(i==Ic) then
  Fcee=F(2,j)
  dkcee=dk(2,j)
 else if(i==Ic-1) then
  if(Is>1) then
   Fcee=F(i+1,j)
   dkcee=dk(i+1,j)
  else
   Fcee=F(1,j)
   dkcee=dk(1,j)
  end if
 else
  Fcee=F(i+2,j)
  dkcee=dk(i+2,j)
 end if
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
  dacs=da(i,j)
  dacss=da(i,j)
  Xcc(4)=Xc(i,j)
  Ycc(4)=Yc(i,j)
  Gxup(4)=Fx(i,j)
  Gyup(4)=Fy(i,j)
 else if(j==1) then
  Fcs=F(Ic+1-i,j)
  Fcss=F(Ic+1-i,j+1)
  dacs=da(Ic+1-i,j)
  dacss=da(Ic+1-i,j+1)
  Xcc(4)=Xc(Ic+1-i,j)
  Ycc(4)=Yc(Ic+1-i,j)
  Gxup(4)=Fx(Ic+1-i,j)
  Gyup(4)=Fy(Ic+1-i,j)
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
  dacs=da(i,j-1)
  dacss=da(i,j-1)
  Xcc(4)=Xc(i,j-1)
  Ycc(4)=Yc(i,j-1)
  Gxup(4)=Fx(i,j-1)
  Gyup(4)=Fy(i,j-1)
 else if(j==2) then
  Fcs=F(i,j-1)
  Fcss=F(Ic+1-i,j-1)
  dacs=da(i,j-1)
  dacss=da(Ic+1-i,j-1)
  Xcc(4)=Xc(i,j-1)
  Ycc(4)=Yc(i,j-1)
  Gxup(4)=Fx(i,j-1)
  Gyup(4)=Fy(i,j-1)
 else
  Fcs=F(i,j-1)
  Fcss=F(i,j-2)
  dacs=da(i,j-1)
  dacss=da(i,j-2)
  Xcc(4)=Xc(i,j-1)
  Ycc(4)=Yc(i,j-1)
  Gxup(4)=Fx(i,j-1)
  Gyup(4)=Fy(i,j-1)
 end if
 if(j==Jc-1) then
  Fcnn=F(i,j+1)
  dacnn=da(i,j+1)
 else
  Fcnn=F(i,j+2)
  dacnn=da(i,j+2)
 end if
 if(isUp) then
  cor(i,j)=0
 else if(is2ndUp) then
  Psip=dkcw*(Fcw-Fcww)/(dkcw+dkcww)
  Psim=dkcc*(Fce-Fcc)/(dkcc+dkce)
  cor(i,j)=max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim
  Psip=dkcc*(Fcw-Fcc)/(dkcc+dkcw)
  Psim=dkce*(Fce-Fcee)/(dkce+dkcee)
  cor(i,j)=cor(i,j)+max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim
  Psip=dacs*(Fcs-Fcss)/(dacs+dacss)
  Psim=dacc*(Fcn-Fcc)/(dacc+dacn)
  cor(i,j)=cor(i,j)+max(Fs(i,j),0.0)*Psip+max(-Fs(i,j),0.0)*Psim
  Psip=dacc*(Fcs-Fcc)/(dacc+dacs)
  Psim=dacn*(Fcn-Fcnn)/(dacn+dacnn)
  cor(i,j)=cor(i,j)+max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim
 else if(isFromm) then
  Xfc=0.5*(Xg(i,j)+Xg(i,j+1))
  Yfc=0.5*(Yg(i,j)+Yg(i,j+1))
  Psip=Gxup(2)*(Xfc-Xcc(2))+Gyup(2)*(Yfc-Ycc(2))
  Psim=-Gxup(1)*(Xfc-Xcc(1))-Gyup(1)*(Yfc-Ycc(1))
  cor(i,j)=max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim
  Xfc=0.5*(Xg(i+1,j)+Xg(i+1,j+1))
  Yfc=0.5*(Yg(i+1,j)+Yg(i+1,j+1))
  Psip=-Gxup(1)*(Xfc-Xcc(1))-Gyup(1)*(Yfc-Ycc(1))
  Psim=Gxup(3)*(Xfc-Xcc(3))+Gyup(3)*(Yfc-Ycc(3))
  cor(i,j)=cor(i,j)+max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim
  Xfc=0.5*(Xg(i,j)+Xg(i+1,j))
  Yfc=0.5*(Yg(i,j)+Yg(i+1,j))
  Psip=Gxup(4)*(Xfc-Xcc(4))+Gyup(4)*(Yfc-Ycc(4))
  Psim=-Gxup(1)*(Xfc-Xcc(1))-Gyup(1)*(Yfc-Ycc(1))
  cor(i,j)=cor(i,j)+max(Fs(i,j),0.0)*Psip+max(-Fs(i,j),0.0)*Psim
  Xfc=0.5*(Xg(i,j+1)+Xg(i+1,j+1))
  Yfc=0.5*(Yg(i,j+1)+Yg(i+1,j+1))
  Psip=-Gxup(1)*(Xfc-Xcc(1))-Gyup(1)*(Yfc-Ycc(1))
  Psim=Gxup(5)*(Xfc-Xcc(5))+Gyup(5)*(Yfc-Ycc(5))
  cor(i,j)=cor(i,j)+max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim
 else if(isQuick) then
  if(i==2) then
   Psip=0
  else
   Psip=dkcw*(3*(Fcc-Fcw)/(dkcc+dkcw)+(Fcw-Fcww)/(dkcw+dkcww))/4
   !Psip=(3*Fcc-2*Fcw-Fcww)/8
  end if
  Psim=dkcc*(3*(Fcw-Fcc)/(dkcc+dkcw)+(Fcc-Fce)/(dkcc+dkce))/4
  !Psim=(3*Fcw-2*Fcc-Fce)/8
  cor(i,j)=max(Fw(i,j),0.0)*Psip+max(-Fw(i,j),0.0)*Psim
  if(i==Ic-1) then
   Psim=0
  else
   Psim=dkce*(3*(Fce-Fcc)/(dkcc+dkce)+(Fcee-Fce)/(dkce+dkcee))/4
   !Psim=(2*Fce+Fcee-3*Fcc)/8
  end if
  Psip=dkcc*(3*(Fcc-Fce)/(dkcc+dkce)+(Fcw-Fcc)/(dkcc+dkcw))/4
  !Psip=(2*Fcc+Fcw-3*Fce)/8
  cor(i,j)=cor(i,j)+max(Fe(i,j),0.0)*Psip+max(-Fe(i,j),0.0)*Psim
  if(j==Jc-1) then
   Psim=0
  else
   Psim=dacn*(3*(Fcn-Fcc)/(dacc+dacn)+(Fcnn-Fcn)/(dacn+dacnn))/4
   !Psim=(2*Fcn+Fcnn-3*Fcc)/8
  end if
  Psip=dacc*(3*(Fcc-Fcn)/(dacc+dacn)+(Fcs-Fcc)/(dacc+dacs))/4
  !Psip=(2*Fcc+Fcs-3*Fcn)/8
  cor(i,j)=cor(i,j)+max(Fn(i,j),0.0)*Psip+max(-Fn(i,j),0.0)*Psim
  Psip=dacs*(3*(Fcc-Fcs)/(dacc+dacs)+(Fcs-Fcss)/(dacs+dacss))/4
  !Psip=(3*Fcc-2*Fcs-Fcss)/8
  Psim=dacc*(3*(Fcs-Fcc)/(dacc+dacs)+(Fcc-Fcn)/(dacc+dacn))/4
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