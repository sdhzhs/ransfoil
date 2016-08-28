Subroutine Velcorrect
use Aero2DCOM
implicit none
integer i,j
if(solctrl=='SIMPLE') then
   DO j=1,Jc-1
     DO i=2,Ic-1
        U(i,j)=U(i,j)-Rau*dPx(i,j)*Jg(i,j)*dx*dy/auP(i,j)
        V(i,j)=V(i,j)-Rau*dPy(i,j)*Jg(i,j)*dx*dy/auP(i,j)
     end DO
   end DO
else if(solctrl=='SIMPLEC') then
    DO j=1,Jc-1
     DO i=2,Ic-1
        U(i,j)=U(i,j)-Rau*dPx(i,j)*Jg(i,j)*dx*dy/(auP(i,j)-Rau*auNB(i,j))
        V(i,j)=V(i,j)-Rau*dPy(i,j)*Jg(i,j)*dx*dy/(auP(i,j)-Rau*auNB(i,j))
     end DO
    end DO
end if
DO j=1,Jc-1
  DO i=2,Ic-1
   Unw(i,j)=Unw(i,j)+wdu(i,j)*(dP(i-1,j)-dP(i,j))
   Une(i,j)=Une(i,j)+edu(i,j)*(dP(i,j)-dP(i+1,j))
   if(j==1.and.(i>Ib2.or.i<Ib1)) then
   Vns(i,j)=Vns(i,j)+sdv(i,j)*(dP(Ic+1-i,j)-dP(i,j))
   else if(j==1) then
   Vns(i,j)=Vns(i,j)
   else
   Vns(i,j)=Vns(i,j)+sdv(i,j)*(dP(i,j-1)-dP(i,j))
   end if
   Vnn(i,j)=Vnn(i,j)+ndv(i,j)*(dP(i,j)-dP(i,j+1))
  end DO
end DO
DO j=1,Jc
   DO i=1,Ic
   Un(i,j)=U(i,j)*Yga(i,j)-V(i,j)*Xga(i,j)
   end DO
end DO
DO j=1,Jc
   DO i=1,Ic
   Vn(i,j)=V(i,j)*Xgk(i,j)-U(i,j)*Ygk(i,j)
   end DO
end DO
end Subroutine Velcorrect
