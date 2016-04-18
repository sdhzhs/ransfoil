Subroutine Jacobi
use Aero2DCOM
implicit none
integer i,j,k
real(8) theta(Ic,Jc),dw(Ib1:Ib2)
dx=1
dy=1
Xg=Xg*c
Yg=Yg*c
DO j=1,Jc
  DO i=1,Ic
  Xc(i,j)=(Xg(i,j)+Xg(i+1,j)+Xg(i,j+1)+Xg(i+1,j+1))/4
  Yc(i,j)=(Yg(i,j)+Yg(i+1,j)+Yg(i,j+1)+Yg(i+1,j+1))/4
  end DO
end DO
DO j=1,Jc
   DO i=1,Ic
   Xgk(i,j)=((Xg(i+1,j)+Xg(i+1,j+1))-(Xg(i,j)+Xg(i,j+1)))/(2*dx)
   Ygk(i,j)=((Yg(i+1,j)+Yg(i+1,j+1))-(Yg(i,j)+Yg(i,j+1)))/(2*dx)
   Xga(i,j)=((Xg(i,j+1)+Xg(i+1,j+1))-(Xg(i,j)+Xg(i+1,j)))/(2*dy)
   Yga(i,j)=((Yg(i,j+1)+Yg(i+1,j+1))-(Yg(i,j)+Yg(i+1,j)))/(2*dy)
   end DO
end DO
Jg=Xgk*Yga-Xga*Ygk
a1=(Xga**2+Yga**2)/Jg
y1=(Xgk**2+Ygk**2)/Jg
b1=(Xgk*Xga+Ygk*Yga)/Jg
theta=180*acos(b1/sqrt(a1*y1))/Pi
DO i=Ib1,Ib2
Xw(i)=0.5*(Xg(i+1,1)+Xg(i,1))
Yw(i)=0.5*(Yg(i+1,1)+Yg(i,1))
Dyp(i)=sqrt((Xc(i,1)-Xw(i))**2+(Yc(i,1)-Yw(i))**2)
DR(i)=sqrt((Xg(i+1,1)-Xg(i,1))**2+(Yg(i+1,1)-Yg(i,1))**2)
end DO
Sw((Ib1+Ib2+1)/2)=sqrt((Xw((Ib1+Ib2+1)/2)-Xw((Ib1+Ib2-1)/2))**2+(Yw((Ib1+Ib2+1)/2)-Yw((Ib1+Ib2-1)/2))**2)/2
Sw((Ib1+Ib2-1)/2)=-Sw((Ib1+Ib2+1)/2)
DO i=(Ib1+Ib2+3)/2,Ib2
Sw(i)=Sw(i-1)+sqrt((Xw(i)-Xw(i-1))**2+(Yw(i)-Yw(i-1))**2)
end DO
DO i=(Ib1+Ib2-3)/2,Ib1,-1
Sw(i)=Sw(i+1)-sqrt((Xw(i+1)-Xw(i))**2+(Yw(i+1)-Yw(i))**2)
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
end Subroutine Jacobi
