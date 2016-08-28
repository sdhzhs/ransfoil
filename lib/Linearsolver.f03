Subroutine sor(aP,aW,aE,aS,aN,b,F,F0,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer maxl,i,j,k,Ic,Jc,Ib1,Ib2
real(8) err,omiga,a,aP(Ic,Jc),aW(Ic,Jc),aE(Ic,Jc),aS(Ic,Jc),aN(Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc),Fo(Ic,Jc),rms(Ic,Jc)
character(*) scalar
maxl=1000
if(scalar=='dP') then
err=1d-4
omiga=1.9
else
err=1d-8
omiga=1
end if
!$OMP PARALLEL
DO k=1,maxl
   !$OMP SINGLE
   Fo=F
   !$OMP END SINGLE
   !$OMP DO
   DO j=1,Jc-1
     DO i=2,Ic-1
         if(j>1) then
         F(i,j)=omiga*(a*(aE(i,j)*F(i+1,j)+aW(i,j)*F(i-1,j)+aS(i,j)*F(i,j-1)+aN(i,j)*F(i,j+1)+b(i,j))/aP(i,j)+(1-a)*F0(i,j))+&
         (1-omiga)*F(i,j)
         else if(j==1.and.(i>Ib2.or.i<Ib1)) then
         F(i,j)=omiga*(a*(aE(i,j)*F(i+1,j)+aW(i,j)*F(i-1,j)+aS(i,j)*F(Ic+1-i,j)+aN(i,j)*F(i,j+1)+b(i,j))/aP(i,j)+(1-a)*F0(i,j))+&
         (1-omiga)*F(i,j)
         else
         if(scalar=='Te'.or.scalar=='Tw') then
         cycle
         else
         F(i,j)=omiga*(a*(aE(i,j)*F(i+1,j)+aW(i,j)*F(i-1,j)+aS(i,j)*F(i,j)+aN(i,j)*F(i,j+1)+b(i,j))/aP(i,j)+(1-a)*F0(i,j))+&
         (1-omiga)*F(i,j)
         end if
         end if
     end DO
   end DO
   !$OMP END DO
   !$OMP DO
   DO j=1,Jc
    DO i=1,Ic
     if(abs(Fo(i,j))>1d-15) then
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

Subroutine CGSTAB(aP,aW,aE,aS,aN,b,F,F0,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer maxl,i,j,k,Ic,Jc,Ib1,Ib2
real(8) err,a,aP(Ic,Jc),aW(Ic,Jc),aE(Ic,Jc),aS(Ic,Jc),aN(Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real(8) alpha,beta,rou,omiga,rou0,sumrmsi,sumvrms,sumptz,sumpts
real(8) rmsi(Ic,Jc),rms(Ic,Jc),p(Ic,Jc),v(Ic,Jc),s(Ic,Jc),t(Ic,Jc),pt(Ic,Jc),y(Ic,Jc),z(Ic,Jc)
character(*) scalar
maxl=1000
if(scalar=='dP') then
err=1d-8
else
err=1d-12
end if
rms=0
   DO j=1,Jc-1
     DO i=2,Ic-1
      if(j>1) then
        rms(i,j)=aE(i,j)*F(i+1,j)+aW(i,j)*F(i-1,j)+aS(i,j)*F(i,j-1)+aN(i,j)*F(i,j+1)+b(i,j)+(1-a)*aP(i,j)*F0(i,j)/a-aP(i,j)*F(i,j)/a
      else if(j==1.and.(i>Ib2.or.i<Ib1)) then
        rms(i,j)=aE(i,j)*F(i+1,j)+aW(i,j)*F(i-1,j)+aS(i,j)*F(Ic+1-i,j)+aN(i,j)*F(i,j+1)+b(i,j)+(1-a)*aP(i,j)*F0(i,j)/a-&
        aP(i,j)*F(i,j)/a
      else
        if(scalar=='Te'.or.scalar=='Tw') then
        cycle
        else
        rms(i,j)=aE(i,j)*F(i+1,j)+aW(i,j)*F(i-1,j)+aS(i,j)*F(i,j)+aN(i,j)*F(i,j+1)+b(i,j)+(1-a)*aP(i,j)*F0(i,j)/a-aP(i,j)*F(i,j)/a
        end if
      end if
     end DO
   end DO
rmsi=rms
p=0
v=p
!$OMP PARALLEL PRIVATE(alpha,beta,rou,omiga,rou0)
alpha=1
rou=1
omiga=1
DO k=1,maxl
rou0=rou
!$OMP WORKSHARE
sumrmsi=sum(rmsi*rms)
!$OMP END WORKSHARE
rou=sumrmsi
beta=(rou*alpha)/(rou0*omiga)
!$OMP WORKSHARE
p=rms+beta*(p-omiga*v)
!y=a*p/aP
!v=y
!$OMP END WORKSHARE
Call Sorprecond(aP/a,aW,aE,aS,aN,p,y,Ic,Jc,Ib1,Ib2,scalar)
!$OMP WORKSHARE
v=y
!$OMP END WORKSHARE
   !$OMP DO
   DO j=1,Jc-1
     DO i=2,Ic-1
      if(j>1) then
        v(i,j)=aP(i,j)*y(i,j)/a-aE(i,j)*y(i+1,j)-aW(i,j)*y(i-1,j)-aS(i,j)*y(i,j-1)-aN(i,j)*y(i,j+1)
      else if(j==1.and.(i>Ib2.or.i<Ib1)) then
        v(i,j)=aP(i,j)*y(i,j)/a-aE(i,j)*y(i+1,j)-aW(i,j)*y(i-1,j)-aS(i,j)*y(Ic+1-i,j)-aN(i,j)*y(i,j+1)
      else
        if(scalar=='Te'.or.scalar=='Tw') then
        cycle
        else
        v(i,j)=aP(i,j)*y(i,j)/a-aE(i,j)*y(i+1,j)-aW(i,j)*y(i-1,j)-aS(i,j)*y(i,j)-aN(i,j)*y(i,j+1)
        end if
      end if
     end DO
   end DO
   !$OMP END DO
!$OMP WORKSHARE
sumvrms=sum(v*rmsi)
!$OMP END WORKSHARE
alpha=rou/sumvrms
!$OMP WORKSHARE
s=rms-alpha*v
!z=a*s/aP
!t=z
!$OMP END WORKSHARE
Call Sorprecond(aP/a,aW,aE,aS,aN,s,z,Ic,Jc,Ib1,Ib2,scalar)
!$OMP WORKSHARE
t=z
!$OMP END WORKSHARE
   !$OMP DO
   DO j=1,Jc-1
     DO i=2,Ic-1
      if(j>1) then
        t(i,j)=aP(i,j)*z(i,j)/a-aE(i,j)*z(i+1,j)-aW(i,j)*z(i-1,j)-aS(i,j)*z(i,j-1)-aN(i,j)*z(i,j+1)
      else if(j==1.and.(i>Ib2.or.i<Ib1)) then
        t(i,j)=aP(i,j)*z(i,j)/a-aE(i,j)*z(i+1,j)-aW(i,j)*z(i-1,j)-aS(i,j)*z(Ic+1-i,j)-aN(i,j)*z(i,j+1)
      else
        if(scalar=='Te'.or.scalar=='Tw') then
        cycle
        else
        t(i,j)=aP(i,j)*z(i,j)/a-aE(i,j)*z(i+1,j)-aW(i,j)*z(i-1,j)-aS(i,j)*z(i,j)-aN(i,j)*z(i,j+1)
        end if
      end if
     end DO
   end DO
   !$OMP END DO
!!$OMP WORKSHARE
!pt=a*t/aP
!!$OMP END WORKSHARE
Call Sorprecond(aP/a,aW,aE,aS,aN,t,pt,Ic,Jc,Ib1,Ib2,scalar)
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

Subroutine Sorprecond(aP,aW,aE,aS,aN,b,F,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer i,j,Ic,Jc,Ib1,Ib2
real(8) aP(Ic,Jc),aW(Ic,Jc),aE(Ic,Jc),aS(Ic,Jc),aN(Ic,Jc),b(Ic,Jc),F(Ic,Jc)
real(8) omiga
real(8),allocatable,save,dimension(:,:)::u,v
character(*) scalar
omiga=1.5
!$OMP WORKSHARE
u=b
!$OMP END WORKSHARE
!$OMP DO
DO j=1,Jc-1
  DO i=2,Ic-1
  if(j>1) then
    u(i,j)=omiga*(b(i,j)+aW(i,j)*u(i-1,j)+aS(i,j)*u(i,j-1))/aP(i,j)
  else if(j==1.and.i<Ib1) then
    u(i,j)=omiga*(b(i,j)+aW(i,j)*u(i-1,j))/aP(i,j)
  else if(j==1.and.i>Ib2) then
    u(i,j)=omiga*(b(i,j)+aW(i,j)*u(i-1,j)+aS(i,j)*u(Ic+1-i,j))/aP(i,j)
  else
    if(scalar=='Te'.or.scalar=='Tw') then
    cycle
    else
    u(i,j)=omiga*(b(i,j)+aW(i,j)*u(i-1,j)+aS(i,j)*u(i,j))/aP(i,j)
    end if
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
  if(j==1.and.i>=Ib1.and.i<=Ib2.and.(scalar=='Te'.or.scalar=='Tw')) then
  cycle
  else
  v(i,j)=(2-omiga)*aP(i,j)*u(i,j)/omiga
  end if
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
    if(j==1.and.i<=Ib2.and.(scalar=='Te'.or.scalar=='Tw')) then
    cycle
    else
    F(i,j)=omiga*(v(i,j)+aE(i,j)*F(i+1,j)+aN(i,j)*F(i,j+1))/aP(i,j)
    end if
  else
    F(i,j)=omiga*(v(i,j)+aE(i,j)*F(i+1,j)+aN(i,j)*F(i,j+1)+aS(i,j)*F(Ic+1-i,j))/aP(i,j)
  end if
  end DO
end DO
!$OMP END DO
end Subroutine Sorprecond
