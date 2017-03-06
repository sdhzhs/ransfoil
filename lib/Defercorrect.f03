Subroutine Defercorrect(F,cor,Fw,Fe,Fs,Fn)
use Aero2DCOM
implicit none
integer i,j
real*8 Fcc,Fcw,Fce,Fcs,Fcn,Fcww,Fcee,Fcss,Fcnn
real*8 rwp,rwm,rep,rem,rsp,rsm,rnp,rnm,Psiwp,Psiwm,Psiep,Psiem,Psisp,Psism,Psinp,Psinm
real*8 F(Ic,Jc),cor(Ic,Jc),Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc)
DO j=1,Jc-1
 DO i=2,Ic-1
 Fcc=F(i,j)
 Fcw=F(i-1,j)
 Fce=F(i+1,j)
 Fcn=F(i,j+1)
 if(i==2) then
  Fcww=F(i-1,j)
 else
  Fcww=F(i-2,j)
 end if
 if(i==Ic-1) then
  Fcee=F(i+1,j)
 else
  Fcee=F(i+2,j)
 end if
 if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
  Fcs=F(i,j)
  Fcss=F(i,j)
 else if(j==1) then
  Fcs=F(Ic+1-i,j)
  Fcss=F(Ic+1-i,j+1)
 else if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
  Fcs=F(i,j-1)
  Fcss=F(i,j-1)
 else if(j==2) then
  Fcs=F(i,j-1)
  Fcss=F(Ic+1-i,j-1)
 else
  Fcs=F(i,j-1)
  Fcss=F(i,j-2)
 end if
 if(j==Jc-1) then
  Fcnn=F(i,j+1)
 else
  Fcnn=F(i,j+2)
 end if
 if(Discret=='1upwind') then
  cor(i,j)=0
 else if(Discret=='2upwind') then
  cor(i,j)=(Fcw-Fcww)*max(Fw(i,j),0.0)/2+(Fce-Fcc)*max(-Fw(i,j),0.0)/2+(Fcw-Fcc)*max(Fe(i,j),0.0)/2+&
  (Fce-Fcee)*max(-Fe(i,j),0.0)/2+(Fcs-Fcss)*max(Fs(i,j),0.0)/2+(Fcn-Fcc)*max(-Fs(i,j),0.0)/2+&
  (Fcs-Fcc)*max(Fn(i,j),0.0)/2+(Fcn-Fcnn)*max(-Fn(i,j),0.0)/2
 else if(Discret=='Quick') then
  if(i==2) then
   Fcww=3*F(i,j)-2*F(i-1,j)
  end if
  if(i==Ic-1) then
   Fcee=3*F(i,j)-2*F(i+1,j)
  end if
  if(j==2.and.(i>=Ib1.and.i<=Ib2)) then
   Fcss=3*F(i,j)-2*F(i,j-1)
  end if
  if(j==Jc-1) then
   Fcnn=3*F(i,j)-2*F(i,j+1)
  end if
  cor(i,j)=(3*Fcc-2*Fcw-Fcww)*max(Fw(i,j),0.0)/8+(Fcw+2*Fcc-3*Fce)*max(Fe(i,j),0.0)/8+&
  (3*Fcw-2*Fcc-Fce)*max(-Fw(i,j),0.0)/8+(Fcee+2*Fce-3*Fcc)*max(-Fe(i,j),0.0)/8+&
  (3*Fcc-2*Fcs-Fcss)*max(Fs(i,j),0.0)/8+(Fcs+2*Fcc-3*Fcn)*max(Fn(i,j),0.0)/8+&
  (3*Fcs-2*Fcc-Fcn)*max(-Fs(i,j),0.0)/8+(Fcnn+2*Fcn-3*Fcc)*max(-Fn(i,j),0.0)/8
 else if(Discret=='tvd') then
  if(abs(Fcc-Fcw)<1e-30) then
  rwp=0
  else
  rwp=(Fcw-Fcww)/(Fcc-Fcw)
  end if
  if(abs(Fcc-Fcw)<1e-30) then
  rwm=0
  else
  rwm=(Fce-Fcc)/(Fcc-Fcw)
  end if
  if(abs(Fce-Fcc)<1e-30) then
  rep=0
  else
  rep=(Fcc-Fcw)/(Fce-Fcc)
  end if
  if(abs(Fce-Fcc)<1e-30) then
  rem=0
  else
  rem=(Fcee-Fce)/(Fce-Fcc)
  end if
  if(abs(Fcc-Fcs)<1e-30) then
  rsp=0
  else
  rsp=(Fcs-Fcss)/(Fcc-Fcs)
  end if
  if(abs(Fcc-Fcs)<1e-30) then
  rsm=0
  else
  rsm=(Fcn-Fcc)/(Fcc-Fcs)
  end if
  if(abs(Fcn-Fcc)<1e-30) then
  rnp=0
  else
  rnp=(Fcc-Fcs)/(Fcn-Fcc)
  end if
  if(abs(Fcn-Fcc)<1e-30) then
  rnm=0
  else
  rnm=(Fcnn-Fcn)/(Fcn-Fcc)
  end if
  Psiwp=(rwp+rwp**2)/(1+rwp**2)
  Psiwm=(rwm+rwm**2)/(1+rwm**2)
  Psiep=(rep+rep**2)/(1+rep**2)
  Psiem=(rem+rem**2)/(1+rem**2)
  Psisp=(rsp+rsp**2)/(1+rsp**2)
  Psism=(rsm+rsm**2)/(1+rsm**2)
  Psinp=(rnp+rnp**2)/(1+rnp**2)
  Psinm=(rnm+rnm**2)/(1+rnm**2)
  cor(i,j)=-(max(-Fe(i,j),0.0)*Psiem+max(Fe(i,j),0.0)*Psiep)*(Fce-Fcc)/2+(max(Fw(i,j),0.0)*Psiwp+&
  max(-Fw(i,j),0.0)*Psiwm)*(Fcc-Fcw)/2-(max(-Fn(i,j),0.0)*Psinm+max(Fn(i,j),0.0)*Psinp)*(Fcn-Fcc)/2+&
  (max(Fs(i,j),0.0)*Psisp+max(-Fs(i,j),0.0)*Psism)*(Fcc-Fcs)/2
 end if
 end DO
end DO
end Subroutine Defercorrect
