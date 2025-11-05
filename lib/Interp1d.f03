Subroutine interp1(x,y,x0,y0,M,N,order,bctype)
implicit none
integer i,k,M,N
integer::In=1
integer Ierr
real(8) x(N),y(N),x0(M),y0(M),b(M),c(M),d(M)
real(8),external::seval
character(*) order,bctype
if(order=='linear') then
 DO i=1,N
  if(x(i)>=x0(1).and.x(i)<=x0(M)) then
   DO k=1,M-1
    if(x(i)>=x0(k).and.x(i)<=x0(k+1)) then
     In=k
     exit
    end if
   end DO
  else if(x(i)<x0(1)) then
   In=1
  else if(x(i)>x0(M)) then
   In=M-1
  end if
  y(i)=(y0(In)*(x0(In+1)-x(i))+y0(In+1)*(x(i)-x0(In)))/(x0(In+1)-x0(In))
 end DO
else if(order=='spline') then
 if(bctype=='nature') then
  Call CSFIT(M,x0,y0,2,0d+0,2,0d+0,b,c,d,Ierr)
  Call CSEVAL(M,x0,y0,N,x,b,c,d,y)
 else if(bctype=='cubic4') then
  Call CSFIT(M,x0,y0,0,0d+0,0,0d+0,b,c,d,Ierr)
  Call CSEVAL(M,x0,y0,N,x,b,c,d,y)
 else if(bctype=='cyclic') then
  Call CSFIT(M,x0,y0,4,0d+0,4,0d+0,b,c,d,Ierr)
  Call CSEVAL(-M,x0,y0,N,x,b,c,d,y)
 end if
end if
end Subroutine interp1

real(8) function interpl(phi0,phi1,weight)
implicit none
real(8) phi0,phi1,weight

interpl=weight*phi0+(1-weight)*phi1

end function interpl
