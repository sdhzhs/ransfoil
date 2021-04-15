Subroutine BlockTDMA(Ap,Aw,Ae,B,X,M,N)
implicit none
integer k,M,N
real(8) Ap(M,M,N),Aw(M,M,N),Ae(M,M,N),B(M,N),X(M,N),c(M,M,N),d(M,N),invAp(M,M),Dp(M,M),invDp(M,M),Awc(M,M),Awd(M),cX(M)

c=0
d=0
Call invmatrix(invAp,Ap(:,:,1),M)
Call Matrixmultiply(invAp,Ae(:,:,1),c(:,:,1),M,M)
Call Matrixmultiply(invAp,B(:,1),d(:,1),M,1)
DO k=2,N
 Call Matrixmultiply(Aw(:,:,k),c(:,:,k-1),Awc,M,M)
 Call Matrixmultiply(Aw(:,:,k),d(:,k-1),Awd,M,1)
 Dp=Ap(:,:,k)-Awc
 Call invmatrix(invDp,Dp,M)
 Call Matrixmultiply(invDp,Ae(:,:,k),c(:,:,k),M,M)
 Call Matrixmultiply(invDp,B(:,k)+Awd,d(:,k),M,1)
end DO
X(:,N)=d(:,N)
DO k=N-1,1,-1
 Call Matrixmultiply(c(:,:,k),X(:,k+1),cX,M,1)
 X(:,k)=d(:,k)+cX
end DO

end Subroutine BlockTDMA

Subroutine invMatrix(invA,A,N)
implicit none
integer i,j,k,N,kp,iq
real(8) R,maxa,swap
real(8) invA(N,N),A(N,N),B(N,N)

B=0
DO i=1,N
 B(i,i)=1
end DO
DO k=1,N-1
 maxa=0
 kp=k
 DO i=k,N
  if(abs(A(i,k))>maxa) then
   maxa=abs(A(i,k))
   kp=i
  end if
 end DO
 if(kp/=k) then
  DO j=k,N
   swap=A(k,j)
   A(k,j)=A(kp,j)
   A(kp,j)=swap
  end DO
  DO j=1,N
   swap=B(k,j)
   B(k,j)=B(kp,j)
   B(kp,j)=swap
  end DO
 end if
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
  invA(i,k)=B(i,k)
  iq=i+1
  DO j=iq,N
   invA(i,k)=invA(i,k)-A(i,j)*invA(j,k)
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
DO j=1,M
 DO k=1,N
  DO i=1,N
   C(i,j)=C(i,j)+A(i,k)*B(k,j)
  end DO
 end DO
end DO

end subroutine Matrixmultiply
