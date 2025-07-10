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

Subroutine BlockCTDMA(Ap,Aw,Ae,B,X,M,N)
implicit none
integer k,M,N
real(8) Ap(M,M,N),Aw(M,M,N),Ae(M,M,N),B(M,N),X(M,N),c(M,M,N),d(M,N),e(M,M,N),P(M,M,N),Q(M,N),IDA(M,M),invAp(M,M),Dp(M,M),invDp(M,M),&
invEp(M,M),Awc(M,M),Awd(M),Awe(M,M),cc(M,M),cd(M),ce(M,M),invce(M,M),Pc(M,M),Pd(M),Pe(M,M),invPe(M,M),invP(M,M),cX(M),eX(M)
c=0
d=0
e=0
IDA=0
DO k=1,M
IDA(k,k)=1
end DO
Call invmatrix(invAp,Ap(:,:,1),M)
Call Matrixmultiply(invAp,Ae(:,:,1),c(:,:,1),M,M)
Call Matrixmultiply(invAp,B(:,1),d(:,1),M,1)
Call Matrixmultiply(invAp,Aw(:,:,1),e(:,:,1),M,M)
DO k=2,N-1
Call Matrixmultiply(Aw(:,:,k),c(:,:,k-1),Awc,M,M)
Call Matrixmultiply(Aw(:,:,k),d(:,k-1),Awd,M,1)
Call Matrixmultiply(Aw(:,:,k),e(:,:,k-1),Awe,M,M)
Dp=Ap(:,:,k)-Awc
Call invmatrix(invDp,Dp,M)
Call Matrixmultiply(invDp,Ae(:,:,k),c(:,:,k),M,M)
Call Matrixmultiply(invDp,B(:,k)+Awd,d(:,k),M,1)
Call Matrixmultiply(invDp,Awe,e(:,:,k),M,M)
end DO
Call Matrixmultiply(Aw(:,:,N),c(:,:,N-1),Awc,M,M)
Call Matrixmultiply(Aw(:,:,N),d(:,N-1),Awd,M,1)
Call Matrixmultiply(Aw(:,:,N),e(:,:,N-1),Awe,M,M)
Call invmatrix(invEp,Ap(:,:,N)-Awc-Awe,M)
Call Matrixmultiply(invEp,Ae(:,:,N),c(:,:,N),M,M)
Call Matrixmultiply(invEp,B(:,N)+Awd,d(:,N),M,1)
Call Matrixmultiply(c(:,:,N),c(:,:,1),cc,M,M)
Call Matrixmultiply(c(:,:,N),d(:,1),cd,M,1)
Call Matrixmultiply(c(:,:,N),e(:,:,1),ce,M,M)
Call invmatrix(invce,IDA-ce,M)
Call Matrixmultiply(invce,cc,P(:,:,1),M,M)
Call Matrixmultiply(invce,cd+d(:,N),Q(:,1),M,1)
DO k=2,N-1
Call Matrixmultiply(P(:,:,k-1),c(:,:,k-1),Pc,M,M)
Call Matrixmultiply(P(:,:,k-1),d(:,k-1),Pd,M,1)
Call Matrixmultiply(P(:,:,k-1),e(:,:,k-1),Pe,M,M)
Call invmatrix(invPe,IDA-Pe,M)
Call Matrixmultiply(invPe,Pc,P(:,:,k),M,M)
Call Matrixmultiply(invPe,Pd+Q(:,k-1),Q(:,k),M,1)
end DO
Call invmatrix(invP,IDA-P(:,:,N-1),M)
Call Matrixmultiply(invP,Q(:,N-1),X(:,N),M,1)
DO k=N-1,1,-1
Call Matrixmultiply(c(:,:,k),X(:,k+1),cX,M,1)
Call Matrixmultiply(e(:,:,k),X(:,N),eX,M,1)
X(:,k)=d(:,k)+cX+eX
end DO
end Subroutine BlockCTDMA

Subroutine invMatrix(invA,A,N)
implicit none
integer i,j,k,N,kp,iq
real(8) R,maxa,swap
real(8) invA(N,N),A(N,N)

invA=0
DO i=1,N
 invA(i,i)=1
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
   swap=invA(k,j)
   invA(k,j)=invA(kp,j)
   invA(kp,j)=swap
  end DO
 end if
 kp=k+1
 DO i=kp,N
  R=A(i,k)/A(k,k)
  DO j=kp,N
   A(i,j)=A(i,j)-R*A(k,j)
  end DO
  DO j=1,N
   invA(i,j)=invA(i,j)-R*invA(k,j)
  end DO
 end DO
end DO
DO k=1,N
 invA(N,k)=invA(N,k)/A(N,N)
 DO i=N-1,1,-1
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
