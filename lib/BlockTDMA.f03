Subroutine BlockTDMA(Ap,Aw,Ae,B,X,M,N)
implicit none
integer i,k,M,N
real(8) Ap(M,M,N),Aw(M,M,N),Ae(M,M,N),B(M,N),X(M,N),c(M,M,N),d(M,N),Dp(M,M,N),invAp(M,M,N),invDp(M,M,N),Awc(M,M,N),Awd(M,N),&
cX(M,N)
c=0
d=0
DO k=1,N
Call invmatrix(invAp(:,:,k),Ap(:,:,k),M)
end DO
Call Matrixmultiply(invAp(:,:,1),Ae(:,:,1),c(:,:,1),M,M)
DO k=2,N
Call Matrixmultiply(Aw(:,:,k),c(:,:,k-1),Awc(:,:,k),M,M)
Dp(:,:,k)=Ap(:,:,k)-Awc(:,:,k)
Call invmatrix(invDp(:,:,k),Dp(:,:,k),M)
Call Matrixmultiply(invDp(:,:,k),Ae(:,:,k),c(:,:,k),M,M)
end DO
Call Matrixmultiply(invAp(:,:,1),B(:,1),d(:,1),M,1)
DO k=2,N
Call Matrixmultiply(Aw(:,:,k),d(:,k-1),Awd(:,k),M,1)
Call Matrixmultiply(invDp(:,:,k),B(:,k)+Awd(:,k),d(:,k),M,1)
end DO
DO i=1,M
X(i,N)=d(i,N)
end DO
DO k=N-1,1,-1
Call Matrixmultiply(c(:,:,k),X(:,k+1),cX(:,k),M,1)
X(:,k)=d(:,k)+cX(:,k)
end DO
end Subroutine BlockTDMA

Subroutine invMatrix(invA,A,N)
implicit none
integer i,j,k,N,kp,iq
real(8) invA(N,N),A(N,N),B(N,N),R
B=0
DO i=1,N
B(i,i)=1
end DO
DO k=1,N-1
kp=k+1
DO i=kp,N
R=A(i,k)/A(k,k)
DO j=kp,N
A(i,j)=A(i,j)-R*A(k,j)
end DO
DO j=1,N
B(i,j)=B(i,j)-R*B(k,j)
end DO
end DO
end DO
DO k=1,N
invA(N,k)=B(N,k)/A(N,N)
DO i=N-1,1,-1
iq=i+1
DO j=iq,N
invA(i,k)=B(i,k)-A(i,j)*invA(j,k)
end DO
invA(i,k)=invA(i,k)/A(i,i)
end DO
end DO
end subroutine invmatrix

Subroutine Matrixmultiply(A,B,C,N,M)
implicit none
integer i,j,k,N,M
real(8) A(N,N),B(N,M),C(N,M)
C=0
DO i=1,N
 DO k=1,N
  DO j=1,M
  C(i,j)=C(i,j)+A(i,k)*B(k,j)
  end DO
 end DO
end DO
end subroutine Matrixmultiply
