Subroutine sor(aM,b,F,F0,a,Ic,Jc,Ib1,Ib2,scalar)
implicit none
integer maxl,i,j,k,Ic,Jc,Ib1,Ib2,Is,Ie
real(8) err,omega,a,rms,FW,FE,FS,FN
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real(8) Fo(Ic,Jc)
character(*) scalar
logical(1) isP,isTe,isTw

isP = scalar=='dP'
isTe = scalar=='Te'
isTw = scalar=='Tw'

!$OMP PARALLEL PRIVATE(k,maxl,err,omega,Is,Ie)
maxl=1000
if(isP) then
 err=1e-4
 omega=1.9
else
 err=1e-6
 omega=1
end if
if(Ib1>1.and.Ib2<Ic) then
 Is=2
 Ie=Ic-1
else
 Is=1
 Ie=Ic
end if
DO k=1,maxl
  !$OMP WORKSHARE
  Fo=F
  !$OMP END WORKSHARE
  if(isTe.or.isTw) then
  !$OMP DO
  DO j=1,Jc-1
    DO i=Is,Ie
      if(i==1) then
        FW=F(Ic,j)
      else
        FW=F(i-1,j)
      end if
      if(i==Ic) then
        FE=F(1,j)
      else
        FE=F(i+1,j)
      end if
      if(j==1.and.(i>Ib2.or.i<Ib1)) then
        FS=F(Ic+1-i,j)
      else if(j==1) then
        FS=F(i,j)
      else
        FS=F(i,j-1)
      end if
      FN=F(i,j+1)
      if(j==1.and.(i>Ib2.or.i<Ib1).or.j>1) then
        F(i,j)=omega*(a*(aM(3,i,j)*FE+aM(2,i,j)*FW+aM(4,i,j)*FS+aM(5,i,j)*FN+b(i,j))/aM(1,i,j)+(1-a)*F0(i,j))+(1-omega)*F(i,j)
      else
        cycle
      end if
    end DO
  end DO
  !$OMP END DO
  else
  !$OMP DO
  DO j=1,Jc-1
    DO i=Is,Ie
      if(i==1) then
        FW=F(Ic,j)
      else
        FW=F(i-1,j)
      end if
      if(i==Ic) then
        FE=F(1,j)
      else
        FE=F(i+1,j)
      end if
      if(j==1.and.(i>Ib2.or.i<Ib1)) then
        FS=F(Ic+1-i,j)
      else if(j==1) then
        FS=F(i,j)
      else
        FS=F(i,j-1)
      end if
      FN=F(i,j+1)
      F(i,j)=omega*(a*(aM(3,i,j)*FE+aM(2,i,j)*FW+aM(4,i,j)*FS+aM(5,i,j)*FN+b(i,j))/aM(1,i,j)+(1-a)*F0(i,j))+(1-omega)*F(i,j)
    end DO
  end DO
  !$OMP END DO
  end if
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
integer maxl,k,Ic,Jc,Ib1,Ib2,iter,preid
real(8) err,a,normb
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real(8) alpha,beta,rho,omega,rho0,sumrmsi,sumvrms,sumptz,sumpts,sumrms
real(8) rmsi(Ic,Jc),rms(Ic,Jc),p(Ic,Jc),v(Ic,Jc),s(Ic,Jc),t(Ic,Jc),pt(Ic,Jc),y(Ic,Jc),z(Ic,Jc),aD(Ic,Jc),vec0(Ic,Jc)
character(*) scalar
character(4) pretype
logical(1) isP,isTe,isTw

isP = scalar=='dP'
isTe = scalar=='Te'
isTw = scalar=='Tw'
pretype = 'ILU'

if(pretype=='JAC') then
 preid=0
else if(pretype=='SSOR') then
 preid=1
else if(pretype=='ILU') then
 preid=2
else
 preid=-1
end if

!$OMP PARALLEL PRIVATE(alpha,beta,rho,omega,rho0,k,maxl,err)
maxl=1000
if(isP) then
 err=1e-8
else
 err=1e-10
end if

Call Residual(aM,b,F,F0,rms,a,Ic,Jc,Ib1,Ib2,isTe,isTw)
!$OMP WORKSHARE
rmsi=rms
p=0
v=p
aD=aM(1,:,:)
normb=sum(abs(b))
!$OMP END WORKSHARE
if(preid==2) then
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
 if(preid==0) then
  !$OMP WORKSHARE
  y=a*p/aM(1,:,:)
  !$OMP END WORKSHARE
 else
  !$OMP WORKSHARE
  y=p
  !$OMP END WORKSHARE
 end if
 if(preid==1.or.preid==2) then
!  Call Sorprecond(aM,aD,y,vec0,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
  Call SSorpcond(y,vec0)
 end if
 Call Multmatrixvector(aM,y,v,a,Ic,Jc,Ib1,Ib2,isTe,isTw)
 !$OMP WORKSHARE
 sumvrms=sum(v*rmsi)
 !$OMP END WORKSHARE
 alpha=rho/sumvrms
 !$OMP WORKSHARE
 s=rms-alpha*v
 !$OMP END WORKSHARE
 if(preid==0) then
  !$OMP WORKSHARE
  z=a*s/aM(1,:,:)
  !$OMP END WORKSHARE
 else
  !$OMP WORKSHARE
  z=s
  !$OMP END WORKSHARE
 end if
 if(preid==1.or.preid==2) then
!  Call Sorprecond(aM,aD,z,vec0,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
  Call SSorpcond(z,vec0)
 end if
 Call Multmatrixvector(aM,z,t,a,Ic,Jc,Ib1,Ib2,isTe,isTw)
 if(preid==0) then
  !$OMP WORKSHARE
  pt=a*t/aM(1,:,:)
  !$OMP END WORKSHARE
 else
  !$OMP WORKSHARE
  pt=t
  !$OMP END WORKSHARE
 end if
 if(preid==1.or.preid==2) then
!  Call Sorprecond(aM,aD,pt,vec0,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
  Call SSorpcond(pt,vec0)
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
 !$OMP WORKSHARE
 sumrms=sum(abs(rms))
 !$OMP END WORKSHARE
 !!$OMP SINGLE
 !iter=k
 !!$OMP END SINGLE
 !if(normb==0) then
  if(sumrms/(Ic*Jc)<err) exit
 !else
 ! if(sumrms/normb<err) exit
 !end if
end DO
!$OMP END PARALLEL
!print *,sumrms/(Ic*Jc),iter

contains

Subroutine SSorpcond(vt,vt0)
!$ use omp_lib
implicit none
integer i,j,nt,tid,npt,rpt,Is,Ie
real(8) vt(Ic,Jc),vt0(Ic,Jc)
real(8) omega
logical(1) cond(4)

if(preid==1.and.isP) then
 omega=1.5
else
 omega=1.0
end if
if(Ib1>1.and.Ib2<Ic) then
 Is=2
 Ie=Ic-1
else
 Is=1
 Ie=Ic
end if

nt=1
tid=0
!$ nt=OMP_GET_NUM_THREADS()
!$ tid=OMP_GET_THREAD_NUM()
npt=(Jc-1)/nt
rpt=mod(Jc-1,nt)

if(Ib1>1.and.Ib2<Ic) then
!$OMP WORKSHARE
 vt(1,:)=omega*vt(1,:)
 vt(Ic,:)=omega*vt(Ic,:)
!$OMP END WORKSHARE
end if
!$OMP WORKSHARE
vt(:,Jc)=omega*vt(:,Jc)
!$ vt0=vt
!$OMP END WORKSHARE

cond(1)=(npt==0)
cond(2)=(nt>1.and.tid>0)
if(isTe.or.isTw) then
!$OMP DO
!!$OMP SINGLE
DO j=1,Jc-1
 DO i=Is,Ie
  if(j>1) then
!$   cond(3)=(tid<rpt.and.mod(j,npt+1)==1)
!$   cond(4)=(tid>=rpt.and.(npt==1.or.mod(j-rpt*(npt+1),npt)==1))
!$   if(cond(1).or.cond(2).and.(cond(3).or.cond(4))) then
!$    if(i==1) then
!$     vt(i,j)=a*omega*(vt(i,j)+aM(4,i,j)*vt0(i,j-1))/aD(i,j)
!$    else if(i==Ic) then
!$     vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(3,i,j)*vt(1,j)+aM(4,i,j)*vt0(i,j-1))/aD(i,j)
!$    else
!$     vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(4,i,j)*vt0(i,j-1))/aD(i,j)
!$    end if
!$   else
    if(i==1) then
     vt(i,j)=a*omega*(vt(i,j)+aM(4,i,j)*vt(i,j-1))/aD(i,j)
    else if(i==Ic) then
     vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(3,i,j)*vt(1,j)+aM(4,i,j)*vt(i,j-1))/aD(i,j)
    else
     vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(4,i,j)*vt(i,j-1))/aD(i,j)
    end if
!$   end if
  else if(j==1.and.i>Ib2) then
   vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(4,i,j)*vt(Ic+1-i,j))/aD(i,j)
  else if(i>=Ib1) then
   vt(i,j)=omega*vt(i,j)
  else
   if(i==1) then
    vt(i,j)=a*omega*vt(i,j)/aD(i,j)
   else if(i==Ic) then
    vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(3,i,j)*vt(1,j))/aD(i,j)
   else
    vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j))/aD(i,j)
   end if
  end if
 end DO
end DO
!!$OMP END SINGLE
!$OMP END DO
else
!$OMP DO
!!$OMP SINGLE
DO j=1,Jc-1
 DO i=Is,Ie
  if(j>1) then
!$   cond(3)=(tid<rpt.and.mod(j,npt+1)==1)
!$   cond(4)=(tid>=rpt.and.(npt==1.or.mod(j-rpt*(npt+1),npt)==1))
!$   if(cond(1).or.cond(2).and.(cond(3).or.cond(4))) then
!$    if(i==1) then
!$     vt(i,j)=a*omega*(vt(i,j)+aM(4,i,j)*vt0(i,j-1))/aD(i,j)
!$    else if(i==Ic) then
!$     vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(3,i,j)*vt(1,j)+aM(4,i,j)*vt0(i,j-1))/aD(i,j)
!$    else
!$     vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(4,i,j)*vt0(i,j-1))/aD(i,j)
!$    end if
!$   else
    if(i==1) then
     vt(i,j)=a*omega*(vt(i,j)+aM(4,i,j)*vt(i,j-1))/aD(i,j)
    else if(i==Ic) then
     vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(3,i,j)*vt(1,j)+aM(4,i,j)*vt(i,j-1))/aD(i,j)
    else
     vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(4,i,j)*vt(i,j-1))/aD(i,j)
    end if
!$   end if
  else if(j==1.and.i>Ib2) then
   vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(4,i,j)*vt(Ic+1-i,j))/aD(i,j)
  else
   if(i==1) then
    vt(i,j)=a*omega*vt(i,j)/aD(i,j)
   else if(i==Ic) then
    vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j)+aM(3,i,j)*vt(1,j))/aD(i,j)
   else
    vt(i,j)=a*omega*(vt(i,j)+aM(2,i,j)*vt(i-1,j))/aD(i,j)
   end if
  end if
 end DO
end DO
!!$OMP END SINGLE
!$OMP END DO
end if

!$OMP WORKSHARE
vt=(2-omega)*vt/omega
!$OMP END WORKSHARE

if(isTe.or.isTw) then
!$OMP DO
DO j=1,Jc-1
 DO i=Is,Ie
  if(.not.(j==1.and.i>=Ib1.and.i<=Ib2)) then
   vt(i,j)=aD(i,j)*vt(i,j)/a
  end if
 end DO
end DO
!$OMP END DO
else
!$OMP DO
DO j=1,Jc-1
 DO i=Is,Ie
  vt(i,j)=aD(i,j)*vt(i,j)/a
 end DO
end DO
!$OMP END DO
end if

if(Ib1>1.and.Ib2<Ic) then
!$OMP WORKSHARE
 vt(1,:)=omega*vt(1,:)
 vt(Ic,:)=omega*vt(Ic,:)
!$OMP END WORKSHARE
end if
!$OMP WORKSHARE
vt(:,Jc)=omega*vt(:,Jc)
!$ vt0=vt
!$OMP END WORKSHARE

if(isTe.or.isTw) then
!$OMP DO
!!$OMP SINGLE
DO j=Jc-1,1,-1
 DO i=Ie,Is,-1
  if(j==1.and.i<Ib1) then
!$   if(npt==0.or.npt==1) then
!$    vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(5,i,j)*vt0(i,j+1)+aM(4,i,j)*vt(Ic+1-i,j))/aD(i,j)
!$   else
    vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(5,i,j)*vt(i,j+1)+aM(4,i,j)*vt(Ic+1-i,j))/aD(i,j)
!$   end if
  else if(j==1.and.i>=Ib1.and.i<=Ib2) then
   vt(i,j)=omega*vt(i,j)
  else
!$   cond(3)=(tid<rpt.and.mod(Jc-j,npt+1)==1)
!$   cond(4)=(tid>=rpt.and.(npt==1.or.mod(Jc-j-rpt*(npt+1),npt)==1))
!$   if(cond(1).or.cond(2).and.(cond(3).or.cond(4))) then
!$    if(i==1) then
!$     vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(2,i,j)*vt(Ic,j)+aM(5,i,j)*vt0(i,j+1))/aD(i,j)
!$    else if(i==Ic) then
!$     vt(i,j)=a*omega*(vt(i,j)+aM(5,i,j)*vt0(i,j+1))/aD(i,j)
!$    else
!$     vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(5,i,j)*vt0(i,j+1))/aD(i,j)
!$    end if
!$   else
    if(i==1) then
     vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(2,i,j)*vt(Ic,j)+aM(5,i,j)*vt(i,j+1))/aD(i,j)
    else if(i==Ic) then
     vt(i,j)=a*omega*(vt(i,j)+aM(5,i,j)*vt(i,j+1))/aD(i,j)
    else
     vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(5,i,j)*vt(i,j+1))/aD(i,j)
    end if
!$   end if
  end if
 end DO
end DO
!!$OMP END SINGLE
!$OMP END DO
else
!$OMP DO
!!$OMP SINGLE
DO j=Jc-1,1,-1
 DO i=Ie,Is,-1
  if(j==1.and.i<Ib1) then
!$   if(npt==0.or.npt==1) then
!$    vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(5,i,j)*vt0(i,j+1)+aM(4,i,j)*vt(Ic+1-i,j))/aD(i,j)
!$   else
    vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(5,i,j)*vt(i,j+1)+aM(4,i,j)*vt(Ic+1-i,j))/aD(i,j)
!$   end if
  else
!$   cond(3)=(tid<rpt.and.mod(Jc-j,npt+1)==1)
!$   cond(4)=(tid>=rpt.and.(npt==1.or.mod(Jc-j-rpt*(npt+1),npt)==1))
!$   if(cond(1).or.cond(2).and.(cond(3).or.cond(4))) then
!$    if(i==1) then
!$     vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(2,i,j)*vt(Ic,j)+aM(5,i,j)*vt0(i,j+1))/aD(i,j)
!$    else if(i==Ic) then
!$     vt(i,j)=a*omega*(vt(i,j)+aM(5,i,j)*vt0(i,j+1))/aD(i,j)
!$    else
!$     vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(5,i,j)*vt0(i,j+1))/aD(i,j)
!$    end if
!$   else
    if(i==1) then
     vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(2,i,j)*vt(Ic,j)+aM(5,i,j)*vt(i,j+1))/aD(i,j)
    else if(i==Ic) then
     vt(i,j)=a*omega*(vt(i,j)+aM(5,i,j)*vt(i,j+1))/aD(i,j)
    else
     vt(i,j)=a*omega*(vt(i,j)+aM(3,i,j)*vt(i+1,j)+aM(5,i,j)*vt(i,j+1))/aD(i,j)
    end if
!$   end if
  end if
 end DO
end DO
!!$OMP END SINGLE
!$OMP END DO
end if
end Subroutine SSorpcond

end Subroutine CGSTAB

Subroutine Residual(aM,b,F,F0,rms,a,Ic,Jc,Ib1,Ib2,isTe,isTw)
implicit none
integer i,j,Ic,Jc,Ib1,Ib2,Is,Ie
real(8) a,FW,FE,FS,FN
real(8) aM(5,Ic,Jc),b(Ic,Jc),F(Ic,Jc),F0(Ic,Jc),rms(Ic,Jc)
logical(1) isTe,isTw

if(Ib1>1.and.Ib2<Ic) then
 Is=2
 Ie=Ic-1
else
 Is=1
 Ie=Ic
end if

!$OMP WORKSHARE
rms=0
!$OMP END WORKSHARE
if(isTe.or.isTw) then
!$OMP DO
DO j=1,Jc-1
 DO i=Is,Ie
  if(i==1) then
   FW=F(Ic,j)
  else
   FW=F(i-1,j)
  end if
  if(i==Ic) then
   FE=F(1,j)
  else
   FE=F(i+1,j)
  end if
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   FS=F(Ic+1-i,j)
  else if(j==1) then
   FS=F(i,j)
  else
   FS=F(i,j-1)
  end if
  FN=F(i,j+1)
  if(j==1.and.(i>Ib2.or.i<Ib1).or.j>1) then
   rms(i,j)=aM(3,i,j)*FE+aM(2,i,j)*FW+aM(4,i,j)*FS+aM(5,i,j)*FN+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-aM(1,i,j)*F(i,j)/a
  end if
 end DO
end DO
!$OMP END DO
else
!$OMP DO
DO j=1,Jc-1
 DO i=Is,Ie
  if(i==1) then
   FW=F(Ic,j)
  else
   FW=F(i-1,j)
  end if
  if(i==Ic) then
   FE=F(1,j)
  else
   FE=F(i+1,j)
  end if
  if(j==1.and.(i>Ib2.or.i<Ib1)) then
   FS=F(Ic+1-i,j)
  else if(j==1) then
   FS=F(i,j)
  else
   FS=F(i,j-1)
  end if
  FN=F(i,j+1)
  rms(i,j)=aM(3,i,j)*FE+aM(2,i,j)*FW+aM(4,i,j)*FS+aM(5,i,j)*FN+b(i,j)+(1-a)*aM(1,i,j)*F0(i,j)/a-aM(1,i,j)*F(i,j)/a
 end DO
end DO
!$OMP END DO
end if
end Subroutine Residual

Subroutine DILU(aM,aD,Ic,Jc,Ib2)
implicit none
integer i,j,Ic,Jc,Ib2,Is,Ie
real(8) aM(5,Ic,Jc),aD(Ic,Jc)

if(Ib2<Ic) then
 Is=2
 Ie=Ic-1
else
 Is=1
 Ie=Ic
end if

!$OMP DO
DO j=1,Jc-1
 DO i=Is,Ie
  if(j==1) then
   if(i>Ib2) then
    aD(i,j)=aD(i,j)-(aM(2,i,j)*aM(3,i-1,j)/aD(i-1,j)+aM(4,i,j)*aM(4,Ic+1-i,j)/aD(Ic+1-i,j))
   else if(i==Ic) then
    aD(i,j)=aD(i,j)-(aM(2,i,j)*aM(3,i-1,j)/aD(i-1,j)+aM(3,i,j)*aM(2,1,j)/aD(1,j))
   else if(i/=1) then
    aD(i,j)=aD(i,j)-aM(2,i,j)*aM(3,i-1,j)/aD(i-1,j)
   end if
  else
   if(i==1) then
    aD(i,j)=aD(i,j)-aM(4,i,j)*aM(5,i,j-1)/aD(i,j-1)
   else if(i==Ic) then
    aD(i,j)=aD(i,j)-(aM(2,i,j)*aM(3,i-1,j)/aD(i-1,j)+aM(3,i,j)*aM(2,1,j)/aD(1,j)+aM(4,i,j)*aM(5,i,j-1)/aD(i,j-1))
   else
    aD(i,j)=aD(i,j)-(aM(2,i,j)*aM(3,i-1,j)/aD(i-1,j)+aM(4,i,j)*aM(5,i,j-1)/aD(i,j-1))
   end if
  end if
 end DO
end DO
!$OMP END DO
end Subroutine DILU

Subroutine Multmatrixvector(aM,u,v,a,Ic,Jc,Ib1,Ib2,isTe,isTw)
implicit none
integer i,j,Ic,Jc,Ib1,Ib2,Is,Ie
real(8) a,uw,ue,us,un
real(8) aM(5,Ic,Jc),u(Ic,Jc),v(Ic,Jc)
logical(1) isTe,isTw

if(Ib1>1.and.Ib2<Ic) then
 Is=2
 Ie=Ic-1
else
 Is=1
 Ie=Ic
end if

!$OMP WORKSHARE
v=u
!$OMP END WORKSHARE
if(isTe.or.isTw) then
!$OMP DO
DO j=1,Jc-1
  DO i=Is,Ie
   if(i==1) then
    uw=u(Ic,j)
   else
    uw=u(i-1,j)
   end if
   if(i==Ic) then
    ue=u(1,j)
   else
    ue=u(i+1,j)
   end if
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    us=u(Ic+1-i,j)
   else if(j==1) then
    us=u(i,j)
   else
    us=u(i,j-1)
   end if
   un=u(i,j+1)
   if(j==1.and.(i>Ib2.or.i<Ib1).or.j>1) then
    v(i,j)=aM(1,i,j)*u(i,j)/a-aM(3,i,j)*ue-aM(2,i,j)*uw-aM(4,i,j)*us-aM(5,i,j)*un
   end if
  end DO
end DO
!$OMP END DO
else
!$OMP DO
DO j=1,Jc-1
  DO i=Is,Ie
   if(i==1) then
    uw=u(Ic,j)
   else
    uw=u(i-1,j)
   end if
   if(i==Ic) then
    ue=u(1,j)
   else
    ue=u(i+1,j)
   end if
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    us=u(Ic+1-i,j)
   else if(j==1) then
    us=u(i,j)
   else
    us=u(i,j-1)
   end if
   un=u(i,j+1)
   v(i,j)=aM(1,i,j)*u(i,j)/a-aM(3,i,j)*ue-aM(2,i,j)*uw-aM(4,i,j)*us-aM(5,i,j)*un
  end DO
end DO
!$OMP END DO
end if
end Subroutine Multmatrixvector

Subroutine Sorprecond(aM,aD,b,b0,a,Ic,Jc,Ib1,Ib2,scalar,pretype)
!$ use omp_lib
implicit none
integer i,j,Ic,Jc,Ib1,Ib2,nt,tid,npt,rpt,Is,Ie
character(*) scalar,pretype
real(8) aM(5,Ic,Jc),aD(Ic,Jc),b(Ic,Jc),b0(Ic,Jc)
real(8) a,omega
logical(1) cond(4)

if(pretype=='SSOR'.and.scalar=='dP') then
 omega=1.5
else
 omega=1.0
end if

if(Ib1>1.and.Ib2<Ic) then
 Is=2
 Ie=Ic-1
else
 Is=1
 Ie=Ic
end if

nt=1
tid=0
!$ nt=OMP_GET_NUM_THREADS()
!$ tid=OMP_GET_THREAD_NUM()
npt=(Jc-1)/nt
rpt=mod(Jc-1,nt)

if(Ib1>1.and.Ib2<Ic) then
!$OMP WORKSHARE
 b(1,:)=omega*b(1,:)
 b(Ic,:)=omega*b(Ic,:)
!$OMP END WORKSHARE
end if
!$OMP WORKSHARE
b(:,Jc)=omega*b(:,Jc)
!$ b0=b
!$OMP END WORKSHARE

cond(1)=(npt==0)
cond(2)=(nt>1.and.tid>0)
!$OMP DO
!!$OMP SINGLE
DO j=1,Jc-1
 DO i=Is,Ie
  if(j>1) then
!$   cond(3)=(tid<rpt.and.mod(j,npt+1)==1)
!$   cond(4)=(tid>=rpt.and.(npt==1.or.mod(j-rpt*(npt+1),npt)==1))
!$   if(cond(1).or.cond(2).and.(cond(3).or.cond(4))) then
!$    if(i==1) then
!$     b(i,j)=a*omega*(b(i,j)+aM(4,i,j)*b0(i,j-1))/aD(i,j)
!$    else if(i==Ic) then
!$     b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(3,i,j)*b(1,j)+aM(4,i,j)*b0(i,j-1))/aD(i,j)
!$    else
!$     b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(4,i,j)*b0(i,j-1))/aD(i,j)
!$    end if
!$   else
    if(i==1) then
     b(i,j)=a*omega*(b(i,j)+aM(4,i,j)*b(i,j-1))/aD(i,j)
    else if(i==Ic) then
     b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(3,i,j)*b(1,j)+aM(4,i,j)*b(i,j-1))/aD(i,j)
    else
     b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(4,i,j)*b(i,j-1))/aD(i,j)
    end if
!$   end if
  else if(j==1.and.i>Ib2) then
   b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(4,i,j)*b(Ic+1-i,j))/aD(i,j)
  else if(i>=Ib1.and.(scalar=='Te'.or.scalar=='Tw')) then
   b(i,j)=omega*b(i,j)
  else
   if(i==1) then
    b(i,j)=a*omega*b(i,j)/aD(i,j)
   else if(i==Ic) then
    b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j)+aM(3,i,j)*b(1,j))/aD(i,j)
   else
    b(i,j)=a*omega*(b(i,j)+aM(2,i,j)*b(i-1,j))/aD(i,j)
   end if
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
 DO i=Is,Ie
  if(.not.(j==1.and.i>=Ib1.and.i<=Ib2.and.(scalar=='Te'.or.scalar=='Tw'))) then
   b(i,j)=aD(i,j)*b(i,j)/a
  end if
 end DO
end DO
!$OMP END DO

if(Ib1>1.and.Ib2<Ic) then
!$OMP WORKSHARE
 b(1,:)=omega*b(1,:)
 b(Ic,:)=omega*b(Ic,:)
!$OMP END WORKSHARE
end if
!$OMP WORKSHARE
b(:,Jc)=omega*b(:,Jc)
!$ b0=b
!$OMP END WORKSHARE

!$OMP DO
!!$OMP SINGLE
DO j=Jc-1,1,-1
 DO i=Ie,Is,-1
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
!$    if(i==1) then
!$     b(i,j)=a*omega*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(2,i,j)*b(Ic,j)+aM(5,i,j)*b0(i,j+1))/aD(i,j)
!$    else if(i==Ic) then
!$     b(i,j)=a*omega*(b(i,j)+aM(5,i,j)*b0(i,j+1))/aD(i,j)
!$    else
!$     b(i,j)=a*omega*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(5,i,j)*b0(i,j+1))/aD(i,j)
!$    end if
!$   else
    if(i==1) then
     b(i,j)=a*omega*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(2,i,j)*b(Ic,j)+aM(5,i,j)*b(i,j+1))/aD(i,j)
    else if(i==Ic) then
     b(i,j)=a*omega*(b(i,j)+aM(5,i,j)*b(i,j+1))/aD(i,j)
    else
     b(i,j)=a*omega*(b(i,j)+aM(3,i,j)*b(i+1,j)+aM(5,i,j)*b(i,j+1))/aD(i,j)
    end if
!$   end if
  end if
 end DO
end DO
!!$OMP END SINGLE
!$OMP END DO
end Subroutine Sorprecond