Subroutine Jacobi
use Aero2DCOM
implicit none
integer i,j,k,Imin
real(8) Xmin,Tricx1,Tricy1,Triarea1,Tricx2,Tricy2,Triarea2,d0,d1,d01x,d01y,Sf,sintheta
real(8) theta(Ic,Jc),dw(Ib1:Ib2)

Xg=Xg*c
Yg=Yg*c
DO j=1,Jc
  DO i=1,Ic
   Tricx1=(Xg(i,j)+Xg(i+1,j)+Xg(i,j+1))/3
   Tricy1=(Yg(i,j)+Yg(i+1,j)+Yg(i,j+1))/3
   Triarea1=abs((Xg(i+1,j)-Xg(i,j))*(Yg(i,j+1)-Yg(i,j))-(Yg(i+1,j)-Yg(i,j))*(Xg(i,j+1)-Xg(i,j)))
   Tricx2=(Xg(i+1,j)+Xg(i,j+1)+Xg(i+1,j+1))/3
   Tricy2=(Yg(i+1,j)+Yg(i,j+1)+Yg(i+1,j+1))/3
   Triarea2=abs((Xg(i+1,j)-Xg(i+1,j+1))*(Yg(i,j+1)-Yg(i+1,j+1))-(Yg(i+1,j)-Yg(i+1,j+1))*(Xg(i,j+1)-Xg(i+1,j+1)))
   Xc(i,j)=(Tricx1*Triarea1+Tricx2*Triarea2)/(Triarea1+Triarea2)
   Yc(i,j)=(Tricy1*Triarea1+Tricy2*Triarea2)/(Triarea1+Triarea2)
   Vol(i,j)=0.5*(Triarea1+Triarea2)
   sintheta=Triarea1/sqrt(((Xg(i+1,j)-Xg(i,j))**2+(Yg(i+1,j)-Yg(i,j))**2)*((Xg(i,j+1)-Xg(i,j))**2+(Yg(i,j+1)-Yg(i,j))**2))
   if(sintheta>1.0) sintheta=1.0
   theta(i,j)=180*asin(sintheta)/Pi
  end DO
end DO
DO j=1,Jc
  DO i=1,Ip
   Xfk(i,j)=Yg(i,j+1)-Yg(i,j)
   Yfk(i,j)=-(Xg(i,j+1)-Xg(i,j))
   Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
   if(i==1) then
    if(Is>1) then
     dkw(i,j)=0.0
     dkd(i,j)=2*sqrt((Xc(i,j)-0.5*(Xg(i,j)+Xg(i,j+1)))**2+(Yc(i,j)-0.5*(Yg(i,j)+Yg(i,j+1)))**2)
    else
     d0=sqrt((Xc(Ic,j)-0.5*(Xg(i,j)+Xg(i,j+1)))**2+(Yc(Ic,j)-0.5*(Yg(i,j)+Yg(i,j+1)))**2)
     d1=sqrt((Xc(i,j)-0.5*(Xg(i,j)+Xg(i,j+1)))**2+(Yc(i,j)-0.5*(Yg(i,j)+Yg(i,j+1)))**2)
     dkw(i,j)=d1/(d0+d1)
     dkd(i,j)=sqrt((Xc(i,j)-Xc(Ic,j))**2+(Yc(i,j)-Yc(Ic,j))**2)
     d01x=Xc(i,j)-Xc(Ic,j)
     d01y=Yc(i,j)-Yc(Ic,j)
     dkd(i,j)=max((d01x*Xfk(i,j)+d01y*Yfk(i,j))/Sf,0.05*dkd(i,j))
    end if
   else if(i==Ip) then
    if(Ie<Ic) then
     dkw(i,j)=1.0
     dkd(i,j)=2*sqrt((Xc(i-1,j)-0.5*(Xg(i,j)+Xg(i,j+1)))**2+(Yc(i-1,j)-0.5*(Yg(i,j)+Yg(i,j+1)))**2)
    else
     d0=sqrt((Xc(i-1,j)-0.5*(Xg(i,j)+Xg(i,j+1)))**2+(Yc(i-1,j)-0.5*(Yg(i,j)+Yg(i,j+1)))**2)
     d1=sqrt((Xc(1,j)-0.5*(Xg(i,j)+Xg(i,j+1)))**2+(Yc(1,j)-0.5*(Yg(i,j)+Yg(i,j+1)))**2)
     dkw(i,j)=d1/(d0+d1)
     dkd(i,j)=sqrt((Xc(1,j)-Xc(i-1,j))**2+(Yc(1,j)-Yc(i-1,j))**2)
     d01x=Xc(1,j)-Xc(i-1,j)
     d01y=Yc(1,j)-Yc(i-1,j)
     dkd(i,j)=max((d01x*Xfk(i,j)+d01y*Yfk(i,j))/Sf,0.05*dkd(i,j))
    end if
   else
    d0=sqrt((Xc(i-1,j)-0.5*(Xg(i,j)+Xg(i,j+1)))**2+(Yc(i-1,j)-0.5*(Yg(i,j)+Yg(i,j+1)))**2)
    d1=sqrt((Xc(i,j)-0.5*(Xg(i,j)+Xg(i,j+1)))**2+(Yc(i,j)-0.5*(Yg(i,j)+Yg(i,j+1)))**2)
    dkw(i,j)=d1/(d0+d1)
    dkd(i,j)=sqrt((Xc(i,j)-Xc(i-1,j))**2+(Yc(i,j)-Yc(i-1,j))**2)
    d01x=Xc(i,j)-Xc(i-1,j)
    d01y=Yc(i,j)-Yc(i-1,j)
    dkd(i,j)=max((d01x*Xfk(i,j)+d01y*Yfk(i,j))/Sf,0.05*dkd(i,j))
   end if
  end DO
end DO
DO j=1,Jp
  DO i=1,Ic
   Xfa(i,j)=-(Yg(i+1,j)-Yg(i,j))
   Yfa(i,j)=Xg(i+1,j)-Xg(i,j)
   Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
   if(j==1) then
    if(i>=Ib1.and.i<=Ib2) then
     daw(i,j)=0.0
     dad(i,j)=2*sqrt((Xc(i,j)-0.5*(Xg(i+1,j)+Xg(i,j)))**2+(Yc(i,j)-0.5*(Yg(i+1,j)+Yg(i,j)))**2)
    else
     d0=sqrt((Xc(Ic+1-i,j)-0.5*(Xg(i+1,j)+Xg(i,j)))**2+(Yc(Ic+1-i,j)-0.5*(Yg(i+1,j)+Yg(i,j)))**2)
     d1=sqrt((Xc(i,j)-0.5*(Xg(i+1,j)+Xg(i,j)))**2+(Yc(i,j)-0.5*(Yg(i+1,j)+Yg(i,j)))**2)
     daw(i,j)=d1/(d0+d1)
     dad(i,j)=sqrt((Xc(i,j)-Xc(Ic+1-i,j))**2+(Yc(i,j)-Yc(Ic+1-i,j))**2)
     d01x=Xc(i,j)-Xc(Ic+1-i,j)
     d01y=Yc(i,j)-Yc(Ic+1-i,j)
     dad(i,j)=max((d01x*Xfa(i,j)+d01y*Yfa(i,j))/Sf,0.05*dad(i,j))
    end if
   else if(j==Jp) then
    daw(i,j)=1.0
    dad(i,j)=2*sqrt((Xc(i,j-1)-0.5*(Xg(i+1,j)+Xg(i,j)))**2+(Yc(i,j-1)-0.5*(Yg(i+1,j)+Yg(i,j)))**2)
   else
    d0=sqrt((Xc(i,j-1)-0.5*(Xg(i+1,j)+Xg(i,j)))**2+(Yc(i,j-1)-0.5*(Yg(i+1,j)+Yg(i,j)))**2)
    d1=sqrt((Xc(i,j)-0.5*(Xg(i+1,j)+Xg(i,j)))**2+(Yc(i,j)-0.5*(Yg(i+1,j)+Yg(i,j)))**2)
    daw(i,j)=d1/(d0+d1)
    dad(i,j)=sqrt((Xc(i,j)-Xc(i,j-1))**2+(Yc(i,j)-Yc(i,j-1))**2)
    d01x=Xc(i,j)-Xc(i,j-1)
    d01y=Yc(i,j)-Yc(i,j-1)
    dad(i,j)=max((d01x*Xfa(i,j)+d01y*Yfa(i,j))/Sf,0.05*dad(i,j))
   end if
  end DO
end DO
Xmin=1e+30
Imin=(Ib1+Ib2+1)/2
DO i=Ib1,Ib2
 Xw(i)=0.5*(Xg(i+1,1)+Xg(i,1))
 Yw(i)=0.5*(Yg(i+1,1)+Yg(i,1))
 Yp(i)=sqrt((Xc(i,1)-Xw(i))**2+(Yc(i,1)-Yw(i))**2)
 DR(i)=sqrt((Xg(i+1,1)-Xg(i,1))**2+(Yg(i+1,1)-Yg(i,1))**2)
 if(Xg(i,1)<Xmin) then
  Xmin=Xg(i,1)
  if(Is>1) Imin=i
 end if
end DO
Sw(Imin)=DR(Imin)/2
Sw(Imin-1)=-DR(Imin-1)/2
DO i=Imin+1,Ib2
 Sw(i)=Sw(i-1)+DR(i-1)/2+DR(i)/2
end DO
DO i=Imin-2,Ib1,-1
 Sw(i)=Sw(i+1)-DR(i+1)/2-DR(i)/2
end DO
DO j=1,Jc
 DO i=1,Ic
  DO k=Ib1,Ib2
   dw(k)=sqrt((Xc(i,j)-Xw(k))**2+(Yc(i,j)-Yw(k))**2)
  end DO
  d(i,j)=minval(dw)
 end DO
end DO
print *,'The minimum and maximum angle of grid cells are:'
print *,minval(theta),maxval(theta)
!print *,minval(Xc),maxval(Xc)
!print *,minval(Yc),maxval(Yc)
!print *,minval(Vol),maxval(Vol)
!print *,minval(dkw(2:Ic,:)),maxval(dkw(2:Ic,:))
!print *,minval(dkd(2:Ic,:)),maxval(dkd(2:Ic,:))
!print *,minval(daw(:,2:Jc)),maxval(daw(:,2:Jc))
!print *,minval(dad(:,2:Jc)),maxval(dad(:,2:Jc))
end Subroutine Jacobi
