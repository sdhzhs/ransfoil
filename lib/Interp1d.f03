Subroutine interp1(x,y,x0,y0,M,N,order)
implicit none
integer i,k,M,N
integer::In=1
real(8) x(N),y(N),x0(M),y0(M),b(M),c(M),d(M)
real(8),external::seval
character(*) order
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
 Call spline(M,x0,y0,b,c,d)
 DO i=1,N
  y(i)=seval(M,x(i),x0,y0,b,c,d)
 end DO
end if
end Subroutine interp1
