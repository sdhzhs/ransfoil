Subroutine sor(aM,b,F,F0,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer maxl,i,j,k,Ic,Jc,Ib1,Ib2
real(8) err,omega,a,rms
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real(8) Fo(Ic,Jc)
character(*) scalar
!$OMP PARALLEL PRIVATE(k,maxl,err,omega)
maxl=1000
if(scalar=='dP') then
 err=1e-4
 omega=1.9
else
 err=1e-6
 omega=1
end if
DO k=1,maxl
  !$OMP WORKSHARE
  Fo=F
  !$OMP END WORKSHARE
  !$OMP DO
  DO j=1,Jc-1
    DO i=2,Ic-1
      if(j>1) then
        F(i,j)=omega*(a*(aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j-1)+aM(5,i,j)*F(i,j+1)+b(i,j))/aM(1,i,j)+(1-a)*F0(i,j))+&
        (1-omega)*F(i,j)
      else if(j==1.and.(i>Ib2.or.i<Ib1)) then
        F(i,j)=omega*(a*(aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(Ic+1-i,j)+aM(5,i,j)*F(i,j+1)+b(i,j))/aM(1,i,j)+(1-a)*F0(i,j))+&
        (1-omega)*F(i,j)
      else
        if(scalar=='Te'.or.scalar=='Tw') then
          cycle
        else
          F(i,j)=omega*(a*(aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j)+aM(5,i,j)*F(i,j+1)+b(i,j))/aM(1,i,j)+(1-a)*F0(i,j))+&
          (1-omega)*F(i,j)
        end if
      end if
    end DO
  end DO
  !$OMP END DO
  !$OMP SINGLE
  rms=0
  !$OMP END SINGLE
  !$OMP DO REDUCTION(+:rms)
  DO j=1,Jc
   DO i=1,Ic
    if(abs(Fo(i,j))>0) then
     rms=rms+abs((F(i,j)-Fo(i,j))/Fo(i,j))
    else
     rms=rms+abs(F(i,j)-Fo(i,j))
    end if
   end DO
  end DO
  !$OMP END DO
  if(rms/(Ic*Jc)<err) exit
end DO
!$OMP END PARALLEL
!print *,rms/(Ic*Jc),k
end Subroutine sor

Subroutine CGSTAB(aM,b,F,F0,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer maxl,k,Ic,Jc,Ib1,Ib2
real(8) err,a
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real(8) alpha,beta,rho,omega,rho0,sumrmsi,sumvrms,sumptz,sumpts
real(8) rmsi(Ic,Jc),rms(Ic,Jc),p(Ic,Jc),v(Ic,Jc),s(Ic,Jc),t(Ic,Jc),pt(Ic,Jc),y(Ic,Jc),z(Ic,Jc),aD(Ic,Jc),vec0(Ic,Jc)
character(*) scalar
character(4) pretype

!$OMP PARALLEL PRIVATE(alpha,beta,rho,omega,rho0,k,maxl,err,pretype)
maxl=1000
if(scalar=='dP') then
 err=1e-8
else
 err=1e-10
end if
pretype='ILU'
Call Residual(aM,b,F,F0,rms,a,Ic,Jc,Ib1,Ib2,scalar)
!$OMP WORKSHARE
rmsi=rms
p=0
v=p
aD=aM(1,:,:)
!$OMP END WORKSHARE
if(pretype=='ILU') then
 Call DILU(aM,aD,Ic,Jc,Ib2)
end if
alpha=1
rho=1
omega=1
DO k=1,maxl
 rho0=rho
 !$OMP WORKSHARE
 sumrmsi=sum(rmsi*rms)
 !$OMP END WORKSHARE
 rho=sumrmsi
 beta=(rho*alpha)/(rho0*omega)
 !$OMP WORKSHARE
 p=rms+beta*(p-omega*v)
 !$OMP END WORKSHARE
 if(pretype=='JAC') then
  !$OMP WORKSHARE
  y=a*p/aM(1,:,:)
  !$OMP END WORKSHARE
 else
  !$OMP WORKSHARE
  y=p
  !$OMP END WORKSHARE
 end if
 if(pretype=='ILU'.or.pretype=='SSOR') then
  Call Sorprecond(aM,aD,y,vec0,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
 end if
 Call Multmatrixvector(aM,y,v,a,Ic,Jc,Ib1,Ib2,scalar)
 !$OMP WORKSHARE
 sumvrms=sum(v*rmsi)
 !$OMP END WORKSHARE
 alpha=rho/sumvrms
 !$OMP WORKSHARE
 s=rms-alpha*v
 !$OMP END WORKSHARE
 if(pretype=='JAC') then
  !$OMP WORKSHARE
  z=a*s/aM(1,:,:)
  !$OMP END WORKSHARE
 else
  !$OMP WORKSHARE
  z=s
  !$OMP END WORKSHARE
 end if
 if(pretype=='ILU'.or.pretype=='SSOR') then
  Call Sorprecond(aM,aD,z,vec0,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
 end if
 Call Multmatrixvector(aM,z,t,a,Ic,Jc,Ib1,Ib2,scalar)
 if(pretype=='JAC') then
  !$OMP WORKSHARE
  pt=a*t/aM(1,:,:)
  !$OMP END WORKSHARE
 else
  !$OMP WORKSHARE
  pt=t
  !$OMP END WORKSHARE
 end if
 if(pretype=='ILU'.or.pretype=='SSOR') then
  Call Sorprecond(aM,aD,pt,vec0,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
 end if
 !$OMP WORKSHARE
 sumptz=sum(pt*z)
 sumpts=sum(pt**2)
 !$OMP END WORKSHARE
 omega=sumptz/sumpts
 !$OMP WORKSHARE
 F=F+alpha*y+omega*z
 rms=s-omega*t
 !$OMP END WORKSHARE
 if(sum(abs(rms))/(Ic*Jc)<err) exit
end DO
!$OMP END PARALLEL
!print *,sum(abs(rms))/(Ic*Jc),k
end Subroutine CGSTAB

Subroutine Residual(aM,b,F,F0,rms,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer i,j,Ic,Jc,Ib1,Ib2
real(8) a
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc),rms(Ic,Jc)
character(*) scalar

!$OMP WORKSHARE
rms=0
!$OMP END WORKSHARE
!$OMP DO
DO j=1,Jc-1
 DO i=2,Ic-1
  if(j>1) then
   rms(i,j)=aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j-1)+aM(5,i,j)*F(i,j+1)+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-aM(1,i,j)*F(i,j)/a
  else if(j==1.and.(i>Ib2.or.i<Ib1)) then
   rms(i,j)=aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(Ic+1-i,j)+aM(5,i,j)*F(i,j+1)+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-aM(1,i,j)*F(i,j)/a
  else if(scalar/='Te'.and.scalar/='Tw') then
   rms(i,j)=aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j)+aM(5,i,j)*F(i,j+1)+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-aM(1,i,j)*F(i,j)/a
  end if
 end DO
end DO
!$OMP END DO
end Subroutine Residual

Subroutine DILU(aM,aD,Ic,Jc,Ib2)
implicit none
integer i,j,Ic,Jc,Ib2
real(8) aM(5,Ic,Jc),aD(Ic,Jc)

!$OMP DO
DO j=1,Jc-1
 DO i=2,Ic-1
  if(j==1) then
   if(i>Ib2) then
    aD(i,j)=aD(i,j)-(aM(2,i,j)*aM(3,i-1,j)/aD(i-1,j)+aM(4,i,j)*aM(4,Ic+1-i,j)/aD(Ic+1-i,j))
   else
    aD(i,j)=aD(i,j)-aM(2,i,j)*aM(3,i-1,j)/aD(i-1,j)
   end if
  else
   aD(i,j)=aD(i,j)-(aM(2,i,j)*aM(3,i-1,j)/aD(i-1,j)+aM(4,i,j)*aM(5,i,j-1)/aD(i,j-1))
  end if
 end DO
end DO
!$OMP END DO
end Subroutine DILU

Subroutine Multmatrixvector(aM,u,v,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer i,j,Ic,Jc,Ib1,Ib2
real(8) a
real(8) aM(5,Ic,Jc),u(Ic,Jc),v(Ic,Jc)
character(*) scalar

!$OMP WORKSHARE
v=u
!$OMP END WORKSHARE
!$OMP DO
DO j=1,Jc-1
  DO i=2,Ic-1
   if(j>1) then
    v(i,j)=aM(1,i,j)*u(i,j)/a-aM(3,i,j)*u(i+1,j)-aM(2,i,j)*u(i-1,j)-aM(4,i,j)*u(i,j-1)-aM(5,i,j)*u(i,j+1)
   else if(j==1.and.(i>Ib2.or.i<Ib1)) then
    v(i,j)=aM(1,i,j)*u(i,j)/a-aM(3,i,j)*u(i+1,j)-aM(2,i,j)*u(i-1,j)-aM(4,i,j)*u(Ic+1-i,j)-aM(5,i,j)*u(i,j+1)
   else if(scalar/='Te'.and.scalar/='Tw') then
    v(i,j)=aM(1,i,j)*u(i,j)/a-aM(3,i,j)*u(i+1,j)-aM(2,i,j)*u(i-1,j)-aM(4,i,j)*u(i,j)-aM(5,i,j)*u(i,j+1)
   end if
  end DO
end DO
!$OMP END DO
end Subroutine Multmatrixvector

Subroutine Sorprecond(aM,aD,b,b0,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
!$ use omp_lib
implicit none
integer i,j,Ic,Jc,Ib1,Ib2,nt,tid,npt,rpt
character(*) scalar,pretype
real(8) aM(5,Ic,Jc),aD(Ic,Jc),b(Ic,Jc),b0(Ic,Jc)
real(8) a,omega
logical(1) cond(4)

if(pretype=='SSOR'.and.scalar=='dP') then
 omega=1.5
else
 omega=1.0
end if

nt=1
tid=0
!$ nt=OMP_GET_NUM_THREADS()
!$ tid=OMP_GET_THREAD_NUM()
npt=(Jc-1)/nt
rpt=mod(Jc-1,nt)

!$OMP WORKSHARE
b(1,:)=omega*b(1,:)
b(Ic,:)=omega*b(Ic,:)
b(2:Ic-1,Jc)=omega*b(2:Ic-1,Jc)
!$ b0=b
!$OMP END WORKSHARE

cond(1)=(npt==0)
cond(2)=(nt>1.and.tid>0)
!$OMP DO
!!$OMP SINGLE
DO j=1,Jc-1
 DO i=2,Ic-1
  if(j>1) then
!$   cond(3)=(tid<rpt.and.mod(j,npt+1)==1)
!$   cond(4)=(tid>=rpt.and.(npt==1.or.mod(j-rpt*(npt+1),npt)==1))
!$   if(cond(1).or.cond(2).and.(cond(3).or.cond(4))) then
!$    b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(4,i,j)*b0(i,j-1))/aD(i,j)
!$   else
    b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(4,i,j)*b(i,j-1))/aD(i,j)
!$   end if
  else if(j==1.and.i>Ib2) then
   b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(4,i,j)*b(Ic+1-i,j))/aD(i,j)
  else if(i>=Ib1.and.(scalar=='Te'.or.scalar=='Tw')) then
   b(i,j)=omega*b(i,j)
  else
   b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j))/aD(i,j)
  end if
 end DO
end DO
!!$OMP END SINGLE
!$OMP END DO

!$OMP WORKSHARE
b=(2-omega)*b/omega
!$OMP END WORKSHARE

!$OMP DO
DO j=1,Jc-1
 DO i=2,Ic-1
  if(.not.(j==1.and.i>=Ib1.and.i<=Ib2.and.(scalar=='Te'.or.scalar=='Tw'))) then
   b(i,j)=aD(i,j)*b(i,j)/a
  end if
 end DO
end DO
!$OMP END DO

!$OMP WORKSHARE
b(1,:)=omega*b(1,:)
b(Ic,:)=omega*b(Ic,:)
b(2:Ic-1,Jc)=omega*b(2:Ic-1,Jc)
!$ b0=b
!$OMP END WORKSHARE

!$OMP DO
!!$OMP SINGLE
DO j=Jc-1,1,-1
 DO i=Ic-1,2,-1
  if(j==1.and.i<Ib1) then
!$   if(npt==0.or.npt==1) then
!$    b(i,j)=a*omega*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(5,i,j)*b0(i,j+1)+aM(4,i,j)*b(Ic+1-i,j))/aD(i,j)
!$   else
    b(i,j)=a*omega*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(5,i,j)*b(i,j+1)+aM(4,i,j)*b(Ic+1-i,j))/aD(i,j)
!$   end if
  else if(j==1.and.i>=Ib1.and.i<=Ib2.and.(scalar=='Te'.or.scalar=='Tw')) then
   b(i,j)=omega*b(i,j)
  else
!$   cond(3)=(tid<rpt.and.mod(Jc-j,npt+1)==1)
!$   cond(4)=(tid>=rpt.and.(npt==1.or.mod(Jc-j-rpt*(npt+1),npt)==1))
!$   if(cond(1).or.cond(2).and.(cond(3).or.cond(4))) then
!$    b(i,j)=a*omega*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(5,i,j)*b0(i,j+1))/aD(i,j)
!$   else
    b(i,j)=a*omega*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(5,i,j)*b(i,j+1))/aD(i,j)
!$   end if
  end if
 end DO
end DO
!!$OMP END SINGLE
!$OMP END DO
end Subroutine Sorprecond
