subroutine spline(n, x, y, b, c, d)
!Boundary conditions: The 3rd derivative of the spline at the
!                     left and right endpoints is to match the 3rd deriv-
!                     ative of the cubic passing through the
!                     first or last 4 data points
implicit none
integer n
real(8) x(n), y(n), b(n), c(n), d(n)
integer nm1, ib, i
real(8) t

nm1 = n-1
if(n<2) return
if(n<3) then
  b(1) = (y(2)-y(1))/(x(2)-x(1))
  c(1) = 0.
  d(1) = 0.
  b(2) = b(1)
  c(2) = 0.
  d(2) = 0.
  return
end if

d(1) = x(2) - x(1)
c(2) = (y(2) - y(1))/d(1)
do i = 2, nm1
  d(i) = x(i+1) - x(i)
  b(i) = 2.*(d(i-1) + d(i))
  c(i+1) = (y(i+1) - y(i))/d(i)
  c(i) = c(i+1) - c(i)
end do

b(1) = -d(1)
b(n) = -d(n-1)
c(1) = 0.
c(n) = 0.
if(n>3) then
  c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
  c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
  c(1) = c(1)*d(1)**2/(x(4)-x(1))
  c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
end if

do i = 2, n
  t = d(i-1)/b(i-1)
  b(i) = b(i) - t*d(i-1)
  c(i) = c(i) - t*c(i-1)
end do

c(n) = c(n)/b(n)
do ib = 1, nm1
  i = n-ib
  c(i) = (c(i) - d(i)*c(i+1))/b(i)
end do

b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
do i = 1, nm1
  b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
  d(i) = (c(i+1) - c(i))/d(i)
  c(i) = 3.*c(i)
end do
c(n) = 3.*c(n)
d(n) = d(n-1)
return

end subroutine spline

real(8) function seval(n, u, x, y, b, c, d)
implicit none
integer n
real(8) u, x(n), y(n), b(n), c(n), d(n)
integer i, j, k
real(8) dx

i = 1
if(u>x(2)) then
  j = n+1
  do while(j>i+1)
    k = (i+j)/2
    if (u<x(k)) j = k
    if (u>=x(k)) i = k
  end do
end if

dx = u - x(i)
seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
return

end function seval
