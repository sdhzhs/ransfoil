Subroutine sor(aM,b,F,F0,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer maxl,i,j,k,Ic,Jc,Ib1,Ib2
real(8) err,omiga,a,rms
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real(8) Fo(Ic,Jc)
character(*) scalar
!$OMP PARALLEL PRIVATE(k)
!$OMP SINGLE
maxl=1000
if(scalar=='dP') then
 err=1e-4
 omiga=1.9
else
 err=1e-6
 omiga=1
end if
!$OMP END SINGLE
DO k=1,maxl
  !$OMP WORKSHARE
  Fo=F
  !$OMP END WORKSHARE
  !$OMP DO PRIVATE(i)
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
  !$OMP SINGLE
  rms=0
  !$OMP END SINGLE
  !$OMP DO REDUCTION(+:rms) PRIVATE(i)
  DO j=1,Jc
   DO i=1,Ic
    if(abs(Fo(i,j))>0) then
     rms=rms+abs(F(i,j)-Fo(i,j))/abs(Fo(i,j))
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
integer maxl,i,j,k,Ic,Jc,Ib1,Ib2
real(8) err,a
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real(8) alpha,beta,rho,omiga,rho0,sumrmsi,sumvrms,sumptz,sumpts
real(8) rmsi(Ic,Jc),rms(Ic,Jc),p(Ic,Jc),v(Ic,Jc),s(Ic,Jc),t(Ic,Jc),pt(Ic,Jc),y(Ic,Jc),z(Ic,Jc),aD(Ic,Jc)
character(*) scalar
character(4) pretype

!$OMP PARALLEL PRIVATE(alpha,beta,rho,omiga,rho0,k)
!$OMP SINGLE
maxl=1000
if(scalar=='dP') then
 err=1e-8
else
 err=1e-10
end if
pretype='ILU'
!$OMP END SINGLE
!$OMP WORKSHARE
rms=0
!$OMP END WORKSHARE
!$OMP DO PRIVATE(i)
DO j=1,Jc-1
 DO i=2,Ic-1
  if(j>1) then
   rms(i,j)=aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(i,j-1)+aM(5,i,j)*F(i,j+1)+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-aM(1,i,j)*F(i,j)/a
  else if(j==1.and.(i>Ib2.or.i<Ib1)) then
   rms(i,j)=aM(3,i,j)*F(i+1,j)+aM(2,i,j)*F(i-1,j)+aM(4,i,j)*F(Ic+1-i,j)+aM(5,i,j)*F(i,j+1)+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-&
   aM(1,i,j)*F(i,j)/a
  else
   if(scalar/='Te'.and.scalar/='Tw') then
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
aD=aM(1,:,:)
!$OMP END WORKSHARE
if(pretype=='ILU') then
 !$OMP DO PRIVATE(i)
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
end if
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
 if(pretype=='JAC') then
  y=a*p/aM(1,:,:)
 else
  y=p
 end if
 !$OMP END WORKSHARE
 if(pretype=='ILU'.or.pretype=='SSOR') then
  Call Sorprecond(aM,aD,y,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
 end if
 !$OMP WORKSHARE
 v=y
 !$OMP END WORKSHARE
 !$OMP DO PRIVATE(i)
 DO j=1,Jc-1
  DO i=2,Ic-1
   if(j>1) then
    v(i,j)=aM(1,i,j)*y(i,j)/a-aM(3,i,j)*y(i+1,j)-aM(2,i,j)*y(i-1,j)-aM(4,i,j)*y(i,j-1)-aM(5,i,j)*y(i,j+1)
   else if(j==1.and.(i>Ib2.or.i<Ib1)) then
    v(i,j)=aM(1,i,j)*y(i,j)/a-aM(3,i,j)*y(i+1,j)-aM(2,i,j)*y(i-1,j)-aM(4,i,j)*y(Ic+1-i,j)-aM(5,i,j)*y(i,j+1)
   else
    if(scalar/='Te'.and.scalar/='Tw') then
     v(i,j)=aM(1,i,j)*y(i,j)/a-aM(3,i,j)*y(i+1,j)-aM(2,i,j)*y(i-1,j)-aM(4,i,j)*y(i,j)-aM(5,i,j)*y(i,j+1)
    end if
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
 if(pretype=='JAC') then
  z=a*s/aM(1,:,:)
 else
  z=s
 end if
 !$OMP END WORKSHARE
 if(pretype=='ILU'.or.pretype=='SSOR') then
  Call Sorprecond(aM,aD,z,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
 end if
 !$OMP WORKSHARE
 t=z
 !$OMP END WORKSHARE
 !$OMP DO PRIVATE(i)
 DO j=1,Jc-1
   DO i=2,Ic-1
    if(j>1) then
     t(i,j)=aM(1,i,j)*z(i,j)/a-aM(3,i,j)*z(i+1,j)-aM(2,i,j)*z(i-1,j)-aM(4,i,j)*z(i,j-1)-aM(5,i,j)*z(i,j+1)
    else if(j==1.and.(i>Ib2.or.i<Ib1)) then
     t(i,j)=aM(1,i,j)*z(i,j)/a-aM(3,i,j)*z(i+1,j)-aM(2,i,j)*z(i-1,j)-aM(4,i,j)*z(Ic+1-i,j)-aM(5,i,j)*z(i,j+1)
    else
     if(scalar/='Te'.and.scalar/='Tw') then
      t(i,j)=aM(1,i,j)*z(i,j)/a-aM(3,i,j)*z(i+1,j)-aM(2,i,j)*z(i-1,j)-aM(4,i,j)*z(i,j)-aM(5,i,j)*z(i,j+1)
     end if
    end if
   end DO
 end DO
 !$OMP END DO
 !$OMP WORKSHARE
 if(pretype=='JAC') then
  pt=a*t/aM(1,:,:)
 else
  pt=t
 end if
 !$OMP END WORKSHARE
 if(pretype=='ILU'.or.pretype=='SSOR') then
  Call Sorprecond(aM,aD,pt,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
 end if
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
!print *,sum(abs(rms))/(Ic*Jc),k
end Subroutine CGSTAB

Subroutine Sorprecond(aM,aD,b,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
implicit none
integer i,j,Ic,Jc,Ib1,Ib2
character(*) scalar,pretype
real(8) aM(5,Ic,Jc),aD(Ic,Jc),b(Ic,Jc)
real(8) a,omiga

if(pretype=='SSOR'.and.scalar=='dP') then
 omiga=1.5
else
 omiga=1.0
end if

!$OMP WORKSHARE
b(1,:)=omiga*b(1,:)
b(Ic,:)=omiga*b(Ic,:)
b(2:Ic-1,Jc)=omiga*b(2:Ic-1,Jc)
!$OMP END WORKSHARE

!$OMP DO
DO j=1,Jc-1
 DO i=2,Ic-1
  if(j>1) then
   b(i,j)=a*omiga*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(4,i,j)*b(i,j-1))/aD(i,j)
  else if(j==1.and.i>Ib2) then
   b(i,j)=a*omiga*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(4,i,j)*b(Ic+1-i,j))/aD(i,j)
  else if(i>=Ib1.and.(scalar=='Te'.or.scalar=='Tw')) then
   b(i,j)=omiga*b(i,j)
  else
   b(i,j)=a*omiga*(b(i,j)+aM(2,i,j)*b(i-1,j))/aD(i,j)
  end if
 end DO
end DO
!$OMP END DO

!$OMP WORKSHARE
b=(2-omiga)*b/omiga
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
b(1,:)=omiga*b(1,:)
b(Ic,:)=omiga*b(Ic,:)
b(2:Ic-1,Jc)=omiga*b(2:Ic-1,Jc)
!$OMP END WORKSHARE

!$OMP DO
DO j=Jc-1,1,-1
 DO i=Ic-1,2,-1
  if(j==1.and.i<Ib1) then
   b(i,j)=a*omiga*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(5,i,j)*b(i,j+1)+aM(4,i,j)*b(Ic+1-i,j))/aD(i,j)
  else if(j==1.and.i>=Ib1.and.i<=Ib2.and.(scalar=='Te'.or.scalar=='Tw')) then
   b(i,j)=omiga*b(i,j)
  else
   b(i,j)=a*omiga*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(5,i,j)*b(i,j+1))/aD(i,j)
  end if
 end DO
end DO
!$OMP END DO
end Subroutine Sorprecond
