Subroutine Postprocess
use Aero2DCOM
implicit none
integer i
real(8) Cx,Cy,Ctx,Cty,Cpm
DO i=Ib1,Ib2
if(Turmod=='sa'.and.Walltreat=='lr'.or.(Turmod=='sst'.and.Walltreat=='lr').or.Turmod=='lam'.or.Turmod=='inv') then
if(Ta-Tf/=0) then
hcv(i)=ca*(ka/ca+miut(i,1)/Prt)*(T(i,1)-Tf)/Dyp(i)/(Ta-Tf)
else
hcv(i)=ca*(ka/ca+miut(i,1)/Prt)/Dyp(i)
end if
Ax(i)=(miu(i,1)+miut(i,1))*U(i,1)/Dyp(i)
Ay(i)=(miu(i,1)+miut(i,1))*V(i,1)/Dyp(i)
else if(Turmod=='ke'.or.(Turmod=='sa'.and.Walltreat=='wf').or.(Turmod=='sst'.and.Walltreat=='wf')) then
if(Ta-Tf/=0) then
hcv(i)=Q(i)/(Tf-Ta)
else
hcv(i)=Q(i)/(Tf-T(i,1))
end if
Ax(i)=rou(i,1)*ustar(i)**2*U(i,1)/sqrt(U(i,1)**2+V(i,1)**2)
Ay(i)=rou(i,1)*ustar(i)**2*V(i,1)/sqrt(U(i,1)**2+V(i,1)**2)
end if
end DO
Cx=0
Cy=0
Ctx=0
Cty=0
Cpm=0
DO i=Ib1,Ib2
Cx=Cx+P(i,1)*(Yg(i+1,1)-Yg(i,1))
Cy=Cy+P(i,1)*(-Xg(i+1,1)+Xg(i,1))
Ctx=Ctx+Ax(i)*DR(i)
Cty=Cty+Ay(i)*DR(i)
Cpm=Cpm+(Xw(i)-0.25*c)*P(i,1)*(-Xg(i+1,1)+Xg(i,1))-Yw(i)*P(i,1)*(Yg(i+1,1)-Yg(i,1))
end DO
Xpc=0.25+Cpm/(c*Cy)
Ypc=0
Cl=(Cy*cos(AoA*Pi/180)-Cx*sin(AoA*Pi/180))/(0.5*roui*Vfar**2*c)
Cd=(Cy*sin(AoA*Pi/180)+Cx*cos(AoA*Pi/180))/(0.5*roui*Vfar**2*c)
Cf=(Cty*sin(AoA*Pi/180)+Ctx*cos(AoA*Pi/180))/(0.5*roui*Vfar**2*c)
Cm=-Cpm/(0.5*roui*Vfar**2*c**2)
end Subroutine Postprocess