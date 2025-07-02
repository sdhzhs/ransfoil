Subroutine Velcorrect
use Aero2DCOM
implicit none
integer i,j
if(solctrl=='SIMPLE') then
   DO j=1,Jc-1
     DO i=Is,Ie
       U(i,j)=U(i,j)-Rau*dPx(i,j)*Jg(i,j)*dx*dy/auP(i,j)
       V(i,j)=V(i,j)-Rau*dPy(i,j)*Jg(i,j)*dx*dy/auP(i,j)
     end DO
   end DO
else if(solctrl=='SIMPLEC') then
   DO j=1,Jc-1
     DO i=Is,Ie
       U(i,j)=U(i,j)-Rau*dPx(i,j)*Jg(i,j)*dx*dy/(auP(i,j)-Rau*auNB(i,j))
       V(i,j)=V(i,j)-Rau*dPy(i,j)*Jg(i,j)*dx*dy/(auP(i,j)-Rau*auNB(i,j))
     end DO
   end DO
end if
DO j=1,Jc-1
  DO i=Is,Ie+1
   if(i==1) then
    Unk(i,j)=Unk(i,j)+duk(i,j)*(dP(Ic,j)-dP(i,j))
   else if(i==Ip) then
    Unk(i,j)=Unk(i,j)+duk(i,j)*(dP(i-1,j)-dP(1,j))
   else
    Unk(i,j)=Unk(i,j)+duk(i,j)*(dP(i-1,j)-dP(i,j))
   end if
  end DO
end DO
DO j=1,Jc
  DO i=Is,Ie
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
    Vna(i,j)=Vna(i,j)+dva(i,j)*(dP(Ic+1-i,j)-dP(i,j))
   else if(j==1) then
    Vna(i,j)=Vna(i,j)
   else
    Vna(i,j)=Vna(i,j)+dva(i,j)*(dP(i,j-1)-dP(i,j))
   end if
  end DO
end DO
Un=U*Yga-V*Xga
Vn=V*Xgk-U*Ygk
end Subroutine Velcorrect
