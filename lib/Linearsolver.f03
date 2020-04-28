Subroutine sor(aM,b,F,F0,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer maxl,i,j,k,Ic,Jc,Ib1,Ib2
real(8) err,omiga,a
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc),Fo(Ic,Jc),rms(Ic,Jc)
character(*) scalar
!$OMP PARALLEL
!$OMP SINGLE
maxl=1000
if(scalar=='dP') then
 err=1e-4
 omiga=1.9
else
 err=1e-8
 omiga=1
end if
!$OMP END SINGLE
DO k=1,maxl
  !$OMP WORKSHARE
  Fo=F
  !$OMP END WORKSHARE
  !$OMP DO
  DO j=1,Jc-1
    DO i=2,Ic-1
      if(j>1) then
        F(i,j)=omiga*(a*(aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j-1)+aM(5,i,j)*F(i,j+1)+b(i,j))/aM(1,i,j)+(1-a)*F0(i,j))+&
        (1-omiga)*F(i,j)
      else if(j==1.and.(i>Ib2.or.i<Ib1)) then
        F(i,j)=omiga*(a*(aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(Ic+1-i,j)+aM(5,i,j)*F(i,j+1)+b(i,j))/aM(1,i,j)+(1-a)*F0(i,j))+&
        (1-omiga)*F(i,j)
      else
        if(scalar=='Te'.or.scalar=='Tw') then
          cycle
        else
          F(i,j)=omiga*(a*(aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j)+aM(5,i,j)*F(i,j+1)+b(i,j))/aM(1,i,j)+(1-a)*F0(i,j))+&
          (1-omiga)*F(i,j)
        end if
      end if
    end DO
  end DO
  !$OMP END DO
  !$OMP DO
  DO j=1,Jc
   DO i=1,Ic
    if(abs(Fo(i,j))>1e-15) then
     rms(i,j)=abs(F(i,j)-Fo(i,j))/abs(Fo(i,j))
    else
     rms(i,j)=abs(F(i,j)-Fo(i,j))/(abs(Fo(i,j))+1)
    end if
   end DO
  end DO
  !$OMP END DO
  if(sum(rms)/(Ic*Jc)<err) exit
end DO
!$OMP END PARALLEL
end Subroutine sor

Subroutine CGSTAB(aM,b,F,F0,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer maxl,i,j,k,Ic,Jc,Ib1,Ib2
real(8) err,a
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real(8) alpha,beta,rho,omiga,rho0,sumrmsi,sumvrms,sumptz,sumpts
real(8) rmsi(Ic,Jc),rms(Ic,Jc),p(Ic,Jc),v(Ic,Jc),s(Ic,Jc),t(Ic,Jc),pt(Ic,Jc),y(Ic,Jc),z(Ic,Jc)
character(*) scalar
!$OMP PARALLEL PRIVATE(alpha,beta,rho,omiga,rho0)
!$OMP SINGLE
maxl=1000
if(scalar=='dP') then
 err=1e-8
else
 err=1e-10
end if
!$OMP END SINGLE
!$OMP WORKSHARE
rms=0
!$OMP END WORKSHARE
!$OMP DO
DO j=1,Jc-1
  DO i=2,Ic-1
   if(j>1) then
    rms(i,j)=aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j-1)+aM(5,i,j)*F(i,j+1)+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-aM(1,i,j)*F(i,j)/a
   else if(j==1.and.(i>Ib2.or.i<Ib1)) then
    rms(i,j)=aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(Ic+1-i,j)+aM(5,i,j)*F(i,j+1)+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-&
    aM(1,i,j)*F(i,j)/a
   else
    if(scalar=='Te'.or.scalar=='Tw') then
     cycle
    else
     rms(i,j)=aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j)+aM(5,i,j)*F(i,j+1)+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-aM(1,i,j)*F(i,j)/a
    end if
   end if
  end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
rmsi=rms
p=0
v=p
!$OMP END WORKSHARE
alpha=1
rho=1
omiga=1
DO k=1,maxl
rho0=rho
!$OMP WORKSHARE
sumrmsi=sum(rmsi*rms)
!$OMP END WORKSHARE
rho=sumrmsi
beta=(rho*alpha)/(rho0*omiga)
!$OMP WORKSHARE
p=rms+beta*(p-omiga*v)
!y=a*p/aM(1,:,:)
!v=y
!$OMP END WORKSHARE
Call Sorprecond(aM,p,y,a,Ic,Jc,Ib1,Ib2)
!$OMP WORKSHARE
v=y
!$OMP END WORKSHARE
!$OMP DO
DO j=1,Jc-1
 DO i=2,Ic-1
  if(j>1) then
   v(i,j)=aM(1,i,j)*y(i,j)/a-aM(3,i,j)*y(i+1,j)-aM(2,i,j)*y(i-1,j)-aM(4,i,j)*y(i,j-1)-aM(5,i,j)*y(i,j+1)
  else if(j==1.and.(i>Ib2.or.i<Ib1)) then
   v(i,j)=aM(1,i,j)*y(i,j)/a-aM(3,i,j)*y(i+1,j)-aM(2,i,j)*y(i-1,j)-aM(4,i,j)*y(Ic+1-i,j)-aM(5,i,j)*y(i,j+1)
  else
   v(i,j)=aM(1,i,j)*y(i,j)/a-aM(3,i,j)*y(i+1,j)-aM(2,i,j)*y(i-1,j)-aM(4,i,j)*y(i,j)-aM(5,i,j)*y(i,j+1)
  end if
 end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
sumvrms=sum(v*rmsi)
!$OMP END WORKSHARE
alpha=rho/sumvrms
!$OMP WORKSHARE
s=rms-alpha*v
!z=a*s/aM(1,:,:)
!t=z
!$OMP END WORKSHARE
Call Sorprecond(aM,s,z,a,Ic,Jc,Ib1,Ib2)
!$OMP WORKSHARE
t=z
!$OMP END WORKSHARE
!$OMP DO
DO j=1,Jc-1
  DO i=2,Ic-1
   if(j>1) then
    t(i,j)=aM(1,i,j)*z(i,j)/a-aM(3,i,j)*z(i+1,j)-aM(2,i,j)*z(i-1,j)-aM(4,i,j)*z(i,j-1)-aM(5,i,j)*z(i,j+1)
   else if(j==1.and.(i>Ib2.or.i<Ib1)) then
    t(i,j)=aM(1,i,j)*z(i,j)/a-aM(3,i,j)*z(i+1,j)-aM(2,i,j)*z(i-1,j)-aM(4,i,j)*z(Ic+1-i,j)-aM(5,i,j)*z(i,j+1)
   else
    t(i,j)=aM(1,i,j)*z(i,j)/a-aM(3,i,j)*z(i+1,j)-aM(2,i,j)*z(i-1,j)-aM(4,i,j)*z(i,j)-aM(5,i,j)*z(i,j+1)
   end if
  end DO
end DO
!$OMP END DO
!!$OMP WORKSHARE
!pt=a*t/aM(1,:,:)
!!$OMP END WORKSHARE
Call Sorprecond(aM,t,pt,a,Ic,Jc,Ib1,Ib2)
!$OMP WORKSHARE
sumptz=sum(pt*z)
sumpts=sum(pt**2)
!$OMP END WORKSHARE
omiga=sumptz/sumpts
!$OMP WORKSHARE
F=F+alpha*y+omiga*z
rms=s-omiga*t
!$OMP END WORKSHARE
if(sum(abs(rms))/(Ic*Jc)<err) exit
end DO
!$OMP END PARALLEL
end Subroutine CGSTAB

Subroutine Sorprecond(aM,b,F,a,Ic,Jc,Ib1,Ib2)
implicit none
integer i,j,Ic,Jc,Ib1,Ib2
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc)
real(8) a,omiga
real(8),allocatable,save,dimension(:,:)::u,v
omiga=1.5
!$OMP WORKSHARE
u=b
!$OMP END WORKSHARE
!$OMP DO
DO j=1,Jc-1
  DO i=2,Ic-1
  if(j>1) then
    u(i,j)=a*omiga*(b(i,j)+aM(2,i,j)*u(i-1,j)+aM(4,i,j)*u(i,j-1))/aM(1,i,j)
  else if(j==1.and.i<Ib1) then
    u(i,j)=a*omiga*(b(i,j)+aM(2,i,j)*u(i-1,j))/aM(1,i,j)
  else if(j==1.and.i>Ib2) then
    u(i,j)=a*omiga*(b(i,j)+aM(2,i,j)*u(i-1,j)+aM(4,i,j)*u(Ic+1-i,j))/aM(1,i,j)
  else
    u(i,j)=a*omiga*(b(i,j)+aM(2,i,j)*u(i-1,j)+aM(4,i,j)*u(i,j))/aM(1,i,j)
  end if
  end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
v=u
!$OMP END WORKSHARE
!$OMP DO
DO j=1,Jc-1
  DO i=2,Ic-1
  v(i,j)=(2-omiga)*aM(1,i,j)*u(i,j)/(a*omiga)
  end DO
end DO
!$OMP END DO
!$OMP WORKSHARE
F=v
!$OMP END WORKSHARE
!$OMP DO
DO j=Jc-1,1,-1
  DO i=Ic-1,2,-1
    if(j>1.or.(j==1.and.i>=Ib1)) then
      F(i,j)=a*omiga*(v(i,j)+aM(3,i,j)*F(i+1,j)+aM(5,i,j)*F(i,j+1))/aM(1,i,j)
    else
      F(i,j)=a*omiga*(v(i,j)+aM(3,i,j)*F(i+1,j)+aM(5,i,j)*F(i,j+1)+aM(4,i,j)*F(Ic+1-i,j))/aM(1,i,j)
    end if
  end DO
end DO
!$OMP END DO
end Subroutine Sorprecond
