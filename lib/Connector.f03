Subroutine Connector(X,Y,X0,Y0,fd,ed,Ip,Ip0)
implicit none
integer i,Ip,Ip0
real(8) fd,ed,lenth,ds,s(Ip),s0(Ip0),X(Ip),Y(Ip),X0(Ip0),Y0(Ip0)
s0(1)=0
DO i=1,Ip0-1
 ds=sqrt((X0(i+1)-X0(i))**2+(Y0(i+1)-Y0(i))**2)
 s0(i+1)=s0(i)+ds
end DO
lenth=s0(Ip0)
Call tanhgridline2(lenth,fd,ed,s,Ip)
Call interp1(s,X,s0,X0,Ip0,Ip,'spline')
Call interp1(s,Y,s0,Y0,Ip0,Ip,'spline')
end Subroutine Connector

Subroutine Connector2(X,Y,X0,Y0,fd,ed,Ip,Ip0)
implicit none
integer i,Ip,Ip0
real(8) fd,ed,lenth,ds
real(8) X(Ip),Y(Ip),X0(Ip0),Y0(Ip0)
real(8) s(Ip),s0(Ip0),sd((Ip+1)/2),su((Ip+1)/2)
s0(1)=0
DO i=1,Ip0-1
 ds=sqrt((X0(i+1)-X0(i))**2+(Y0(i+1)-Y0(i))**2)
 s0(i+1)=s0(i)+ds
end DO
lenth=s0((Ip0+1)/2)
Call tanhgridline2(lenth,fd,ed,sd,(Ip+1)/2)
lenth=s0(Ip0)-s0((Ip0+1)/2)
Call tanhgridline2(lenth,fd,ed,su,(Ip+1)/2)
s(1)=0
DO i=1,(Ip-1)/2
 ds=sd((Ip+3)/2-i)-sd((Ip+1)/2-i)
 s(i+1)=s(i)+ds
end DO
DO i=1,(Ip-1)/2
 s((Ip+1)/2+i)=su(i+1)+s((Ip+1)/2)
end DO
Call interp1(s,X,s0,X0,Ip0,Ip,'spline')
Call interp1(s,Y,s0,Y0,Ip0,Ip,'spline')
end Subroutine Connector2

Subroutine tanhgridline(lenth,fb,s,M)
implicit none
integer i,l,maxl,M
real(8) err,lenth,fb,alpha,alpha0,B,s(M)
maxl=100
err=1d-10
B=lenth/((M-1)*fb)
alpha=5
DO l=1,maxl
 alpha0=alpha
 alpha=alpha-(sinh(alpha)/alpha-B)/((cosh(alpha)*alpha-sinh(alpha))/alpha**2)
 if(abs(alpha-alpha0)<err) exit
end DO
DO i=1,M
 s(i)=lenth*(1-tanh(alpha*(M-i)/(M-1)/2)/tanh(alpha/2))
end DO
end Subroutine tanhgridline

Subroutine tanhgridline2(lenth,fb,eb,s,M)
implicit none
integer i,l,maxl,M
real(8) err,lenth,fb,eb,alpha,alpha0,A,B,s(M),u(M)
maxl=100
err=1d-10
A=sqrt(eb/fb)
B=lenth/((M-1)*sqrt(eb*fb))
alpha=5
DO l=1,maxl
 alpha0=alpha
 alpha=alpha-(sinh(alpha)/alpha-B)/((cosh(alpha)*alpha-sinh(alpha))/alpha**2)
 if(abs(alpha-alpha0)<err) exit
end DO
DO i=1,M
 u(i)=0.5*(1-tanh(alpha/2-alpha*(i-1)/(M-1))/tanh(alpha/2))
 s(i)=lenth*u(i)/(A+(1-A)*u(i))
end DO
end Subroutine tanhgridline2