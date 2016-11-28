Subroutine tanhgridline(lenth,fb,s,M)
implicit none
integer i,l,maxl,M
real(8) err,lenth,fb,alpha,alpha0,B,s(M)
maxl=100
err=1e-15
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
