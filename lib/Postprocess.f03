Subroutine Postprocess(iter)
use Aero2DCOM
implicit none
integer i,n,iter
real(8) Cx,Cy,Ctx,Cty,Cpm,Ctm
n=3
if(mod(iter,100)==0) then
 DO i=1,n-1
  Clrec(i)=Clrec(i+1)
  Cdrec(i)=Cdrec(i+1)
  Cfrec(i)=Cfrec(i+1)
  Cmrec(i)=Cmrec(i+1)
  Xpcrec(i)=Xpcrec(i+1)
  Ypcrec(i)=Ypcrec(i+1)
  Pnw(:,i)=Pnw(:,i+1)
  Unw(:,i)=Unw(:,i+1)
  Vnw(:,i)=Vnw(:,i+1)
  Tnw(:,i)=Tnw(:,i+1)
  mutnw(:,i)=mutnw(:,i+1)
  hcvnw(:,i)=hcvnw(:,i+1)
  Axnw(:,i)=Axnw(:,i+1)
  Aynw(:,i)=Aynw(:,i+1)
  Ypnw(:,i)=Ypnw(:,i+1)
 end DO
 Clrec(n)=Cl
 Cdrec(n)=Cd
 Cfrec(n)=Cf
 Cmrec(n)=Cm
 Xpcrec(n)=Xpc
 Ypcrec(n)=Ypc
 Pnw(:,n)=P(Ib1:Ib2,1)
 Unw(:,n)=U(Ib1:Ib2,1)
 Vnw(:,n)=V(Ib1:Ib2,1)
 Tnw(:,n)=T(Ib1:Ib2,1)
 mutnw(:,n)=mut(Ib1:Ib2,1)
 hcvnw(:,n)=hcv
 Axnw(:,n)=Ax
 Aynw(:,n)=Ay
 if(Turmod=='sa'.or.Turmod=='sst') then
  Ypnw(:,n)=Yplus
 else if(Turmod=='ke') then
  Ypnw(:,n)=Ystar
 end if
end if
DO i=Ib1,Ib2
 if(Turmod=='sa'.and.Walltreat=='lr'.or.(Turmod=='sst'.and.Walltreat=='lr').or.Turmod=='lam'.or.Turmod=='inv') then
  if(Ta-Tf/=0) then
   hcv(i)=ca*(ka/ca+mut(i,1)/Prt)*(T(i,1)-Tf)/Yp(i)/(Ta-Tf)
  else
   hcv(i)=ca*(ka/ca+mut(i,1)/Prt)/Yp(i)
  end if
  Ax(i)=(mu(i,1)+mut(i,1))*U(i,1)/Yp(i)
  Ay(i)=(mu(i,1)+mut(i,1))*V(i,1)/Yp(i)
 else if(Turmod=='ke'.or.(Turmod=='sa'.and.Walltreat=='wf').or.(Turmod=='sst'.and.Walltreat=='wf')) then
  if(Ta-Tf/=0) then
   hcv(i)=Q(i)/(Tf-Ta)
  else
   hcv(i)=Q(i)/(Tf-T(i,1))
  end if
  Ax(i)=rho(i,1)*ustar(i)*Un(i,1)*Yga(i,1)/Uplus(i)/da(i,1)**2
  Ay(i)=-rho(i,1)*ustar(i)*Un(i,1)*Xga(i,1)/Uplus(i)/da(i,1)**2
 end if
end DO
Cx=0
Cy=0
Ctx=0
Cty=0
Cpm=0
Ctm=0
DO i=Ib1,Ib2
 Cx=Cx+P(i,1)*(Yg(i+1,1)-Yg(i,1))
 Cy=Cy+P(i,1)*(-Xg(i+1,1)+Xg(i,1))
 Ctx=Ctx+Ax(i)*DR(i)
 Cty=Cty+Ay(i)*DR(i)
 Cpm=Cpm+(Xw(i)-0.25*c)*P(i,1)*(-Xg(i+1,1)+Xg(i,1))-Yw(i)*P(i,1)*(Yg(i+1,1)-Yg(i,1))
 Ctm=Ctm+(Xw(i)-0.25*c)*Ay(i)*DR(i)-Yw(i)*Ax(i)*DR(i)
end DO
Xpc=0.25+Cpm/(c*Cy)
Ypc=0
Cl=((Cty+Cy)*cos(AoA*Pi/180)-(Ctx+Cx)*sin(AoA*Pi/180))/(0.5*rhoi*Vfar**2*c)
Cd=(Cy*sin(AoA*Pi/180)+Cx*cos(AoA*Pi/180))/(0.5*rhoi*Vfar**2*c)
Cf=(Cty*sin(AoA*Pi/180)+Ctx*cos(AoA*Pi/180))/(0.5*rhoi*Vfar**2*c)
Cm=-(Cpm+Ctm)/(0.5*rhoi*Vfar**2*c**2)
end Subroutine Postprocess
