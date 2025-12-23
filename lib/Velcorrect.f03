Subroutine Velcorrect
use Aero2DCOM
implicit none
integer i,j
real(8) Sf
logical(1) isInOut
isInOut=Fstype=='vinpout'

!$OMP PARALLEL
if(solctrl=='SIMPLE') then
   !$OMP DO PRIVATE(i)
   DO j=1,Jc-1
     DO i=Is,Ie
       U(i,j)=U(i,j)-Rau*dPx(i,j)*Vol(i,j)/auP(i,j)
       V(i,j)=V(i,j)-Rau*dPy(i,j)*Vol(i,j)/auP(i,j)
     end DO
   end DO
   !$OMP END DO
else if(solctrl=='SIMPLEC') then
   !$OMP DO PRIVATE(i)
   DO j=1,Jc-1
     DO i=Is,Ie
       U(i,j)=U(i,j)-Rau*dPx(i,j)*Vol(i,j)/(auP(i,j)-Rau*auNB(i,j))
       V(i,j)=V(i,j)-Rau*dPy(i,j)*Vol(i,j)/(auP(i,j)-Rau*auNB(i,j))
     end DO
   end DO
   !$OMP END DO
end if
!$OMP DO PRIVATE(i,Sf)
DO j=1,Jc-1
  DO i=Is,Ie+1
   Sf=sqrt(Xfk(i,j)**2+Yfk(i,j)**2)
   if(i==1) then
    Unk(i,j)=Unk(i,j)+duk(i,j)*(dP(Ic,j)-dP(i,j))*Sf/dkd(i,j)
   else if(i==Ip) then
    Unk(i,j)=Unk(i,j)+duk(i,j)*(dP(i-1,j)-dP(1,j))*Sf/dkd(i,j)
   else
    Unk(i,j)=Unk(i,j)+duk(i,j)*(dP(i-1,j)-dP(i,j))*Sf/dkd(i,j)
   end if
  end DO
end DO
!$OMP END DO
!$OMP DO PRIVATE(i,Sf)
DO j=1,Jc
  DO i=Is,Ie
   Sf=sqrt(Xfa(i,j)**2+Yfa(i,j)**2)
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    Vna(i,j)=Vna(i,j)+dva(i,j)*(dP(Ic+1-i,j)-dP(i,j))*Sf/dad(i,j)
   else if(j==1) then
    Vna(i,j)=Vna(i,j)
   else
    Vna(i,j)=Vna(i,j)+dva(i,j)*(dP(i,j-1)-dP(i,j))*Sf/dad(i,j)
   end if
  end DO
end DO
!$OMP END DO

if(isInOut.and.Is>1) then
 !$OMP WORKSHARE
 U(1,1:Jc-1)=U(2,1:Jc-1)
 U(Ic,1:Jc-1)=U(Ic-1,1:Jc-1)
 V(1,1:Jc-1)=V(2,1:Jc-1)
 V(Ic,1:Jc-1)=V(Ic-1,1:Jc-1)
 !$OMP END WORKSHARE
end if
if(isInOut) then
 !$OMP DO
 DO i=1,Ic
  if(Vna(i,Jc)>0) then
   U(i,Jc)=U(i,Jc-1)
   V(i,Jc)=V(i,Jc-1)
  end if
 end DO
 !$OMP END DO
end if
!$OMP END PARALLEL
end Subroutine Velcorrect
