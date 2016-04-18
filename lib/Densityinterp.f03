Subroutine Densityinterp
use Aero2DCOM
implicit none
integer i,j
real(8) rwp,rwm,rep,rem,rsp,rsm,rnp,rnm,Psiwp,Psiwm,Psiep,Psiem,Psisp,Psism,Psinp,Psinm
  DO j=1,Jc-1
   DO i=2,Ic-1
   we(i,j)=sign(5d-1,Une(i,j))
   ww(i,j)=sign(5d-1,Unw(i,j))
   wn(i,j)=sign(5d-1,Vnn(i,j))
   ws(i,j)=sign(5d-1,Vns(i,j))
   end DO
  end DO
  DO j=1,Jc-1
   DO i=2,Ic-1
   if(Proctrl=='incom') then
   rouw(i,j)=rou(i,j)
   roue(i,j)=rou(i,j)
   roun(i,j)=rou(i,j)
   rous(i,j)=rou(i,j)
   else if(Proctrl=='com') then
   if(denface=='center') then
   rouw(i,j)=0.5*(rou(i,j)+rou(i-1,j))
   roue(i,j)=0.5*(rou(i,j)+rou(i+1,j))
   roun(i,j)=0.5*(rou(i,j)+rou(i,j+1))
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
   rous(i,j)=0.5*(rou(i,j)+rou(Ic+1-i,j))
   else if(j==1) then
   rous(i,j)=rou(i,j)
   else
   rous(i,j)=0.5*(rou(i,j)+rou(i,j-1))
   end if
   else if(denface=='1upwind') then
   rouw(i,j)=(0.5-ww(i,j))*rou(i,j)+(0.5+ww(i,j))*rou(i-1,j)
   roue(i,j)=(0.5+we(i,j))*rou(i,j)+(0.5-we(i,j))*rou(i+1,j)
   roun(i,j)=(0.5+wn(i,j))*rou(i,j)+(0.5-wn(i,j))*rou(i,j+1)
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
   rous(i,j)=(0.5-ws(i,j))*rou(i,j)+(0.5+ws(i,j))*rou(Ic+1-i,j)
   else if(j==1) then
   rous(i,j)=rou(i,j)
   else
   rous(i,j)=(0.5-ws(i,j))*rou(i,j)+(0.5+ws(i,j))*rou(i,j-1)
   end if
   else if(denface=='2upwind') then
   if(i==2) then
   rouw(i,j)=(0.5-ww(i,j))*(1.5*rou(i,j)-0.5*rou(i+1,j))+(0.5+ww(i,j))*rou(i-1,j)
   else
   rouw(i,j)=(0.5-ww(i,j))*(1.5*rou(i,j)-0.5*rou(i+1,j))+(0.5+ww(i,j))*(1.5*rou(i-1,j)-0.5*rou(i-2,j))
   end if
   if(i==Ic-1) then
   roue(i,j)=(0.5+we(i,j))*(1.5*rou(i,j)-0.5*rou(i-1,j))+(0.5-we(i,j))*rou(i+1,j)
   else
   roue(i,j)=(0.5+we(i,j))*(1.5*rou(i,j)-0.5*rou(i-1,j))+(0.5-we(i,j))*(1.5*rou(i+1,j)-0.5*rou(i+2,j))
   end if
   if(j==Jc-1) then
   roun(i,j)=(0.5+wn(i,j))*(1.5*rou(i,j)-0.5*rou(i,j-1))+(0.5-wn(i,j))*rou(i,j+1)
   else if(j==1.and.(i<Ib1.or.i>Ib2)) then
   roun(i,j)=(0.5+wn(i,j))*(1.5*rou(i,j)-0.5*rou(Ic+1-i,j))+(0.5-wn(i,j))*(1.5*rou(i,j+1)-0.5*rou(i,j+2))
   else if(j==1) then
   roun(i,j)=(0.5+wn(i,j))*rou(i,j)+(0.5-wn(i,j))*(1.5*rou(i,j+1)-0.5*rou(i,j+2))
   else
   roun(i,j)=(0.5+wn(i,j))*(1.5*rou(i,j)-0.5*rou(i,j-1))+(0.5-wn(i,j))*(1.5*rou(i,j+1)-0.5*rou(i,j+2))
   end if
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
   rous(i,j)=(0.5-ws(i,j))*(1.5*rou(i,j)-0.5*rou(i,j+1))+(0.5+ws(i,j))*(1.5*rou(Ic+1-i,j)-0.5*rou(Ic+1-i,j+1))
   else if(j==1) then
   rous(i,j)=rou(i,j)
   else if(j==2.and.(i<Ib1.or.i>Ib2)) then
   rous(i,j)=(0.5-ws(i,j))*(1.5*rou(i,j)-0.5*rou(i,j+1))+(0.5+ws(i,j))*(1.5*rou(i,j-1)-0.5*rou(Ic+1-i,j-1))
   else if(j==2) then
   rous(i,j)=(0.5-ws(i,j))*(1.5*rou(i,j)-0.5*rou(i,j+1))+(0.5+ws(i,j))*rou(i,j-1)
   else
   rous(i,j)=(0.5-ws(i,j))*(1.5*rou(i,j)-0.5*rou(i,j+1))+(0.5+ws(i,j))*(1.5*rou(i,j-1)-0.5*rou(i,j-2))
   end if
   else if(denface=='Quick') then
   if(i==2) then
   rouw(i,j)=(0.5-ww(i,j))*(6*rou(i,j)+3*rou(i-1,j)-rou(i+1,j))/8+(0.5+ww(i,j))*rou(i-1,j)
   else
   rouw(i,j)=(0.5-ww(i,j))*(6*rou(i,j)+3*rou(i-1,j)-rou(i+1,j))/8+(0.5+ww(i,j))*(6*rou(i-1,j)+3*rou(i,j)-rou(i-2,j))/8
   end if
   if(i==Ic-1) then
   roue(i,j)=(0.5+we(i,j))*(6*rou(i,j)+3*rou(i+1,j)-rou(i-1,j))/8+(0.5-we(i,j))*rou(i+1,j)
   else
   roue(i,j)=(0.5+we(i,j))*(6*rou(i,j)+3*rou(i+1,j)-rou(i-1,j))/8+(0.5-we(i,j))*(6*rou(i+1,j)+3*rou(i,j)-rou(i+2,j))/8
   end if
   if(j==Jc-1) then
   roun(i,j)=(0.5+wn(i,j))*(6*rou(i,j)+3*rou(i,j+1)-rou(i,j-1))/8+(0.5-wn(i,j))*rou(i,j+1)
   else if(j==1.and.(i<Ib1.or.i>Ib2)) then
   roun(i,j)=(0.5+wn(i,j))*(6*rou(i,j)+3*rou(i,j+1)-rou(Ic+1-i,j))/8+(0.5-wn(i,j))*(6*rou(i,j+1)+3*rou(i,j)-rou(i,j+2))/8
   else if(j==1) then
   roun(i,j)=(0.5+wn(i,j))*rou(i,j)+(0.5-wn(i,j))*(6*rou(i,j+1)+3*rou(i,j)-rou(i,j+2))/8
   else
   roun(i,j)=(0.5+wn(i,j))*(6*rou(i,j)+3*rou(i,j+1)-rou(i,j-1))/8+(0.5-wn(i,j))*(6*rou(i,j+1)+3*rou(i,j)-rou(i,j+2))/8
   end if
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
   rous(i,j)=(0.5-ws(i,j))*(6*rou(i,j)+3*rou(Ic+1-i,j)-rou(i,j+1))/8+(0.5+ws(i,j))*(6*rou(Ic+1-i,j)+3*rou(i,j)-rou(Ic+1-i,j+1))/8
   else if(j==1) then
   rous(i,j)=rou(i,j)
   else if(j==2.and.(i<Ib1.or.i>Ib2)) then
   rous(i,j)=(0.5-ws(i,j))*(6*rou(i,j)+3*rou(i,j-1)-rou(i,j+1))/8+(0.5+ws(i,j))*(6*rou(i,j-1)+3*rou(i,j)-rou(Ic+1-i,j-1))/8
   else if(j==2) then
   rous(i,j)=(0.5-ws(i,j))*(6*rou(i,j)+3*rou(i,j-1)-rou(i,j+1))/8+(0.5+ws(i,j))*rou(i,j-1)
   else
   rous(i,j)=(0.5-ws(i,j))*(6*rou(i,j)+3*rou(i,j-1)-rou(i,j+1))/8+(0.5+ws(i,j))*(6*rou(i,j-1)+3*rou(i,j)-rou(i,j-2))/8
   end if
   else if(denface=='tvd') then
   if(i==2.or.abs(rou(i,j)-rou(i-1,j))<1d-30) then
   rwp=0
   else
   rwp=(rou(i-1,j)-rou(i-2,j))/(rou(i,j)-rou(i-1,j))
   end if
   if(abs(rou(i,j)-rou(i-1,j))<1d-30) then
   rwm=0
   else
   rwm=(rou(i+1,j)-rou(i,j))/(rou(i,j)-rou(i-1,j))
   end if
   if(abs(rou(i+1,j)-rou(i,j))<1d-30) then
   rep=0
   else
   rep=(rou(i,j)-rou(i-1,j))/(rou(i+1,j)-rou(i,j))
   end if
   if(i==Ic-1.or.abs(rou(i+1,j)-rou(i,j))<1d-30) then
   rem=0
   else
   rem=(rou(i+2,j)-rou(i+1,j))/(rou(i+1,j)-rou(i,j))
   end if
   if(j==1.and.(i>=Ib1.and.i<=Ib2).or.abs(rou(i,j)-rou(Ic+1-i,j))<1d-30) then
   rsp=0
   else if(j==1) then
   rsp=(rou(Ic+1-i,j)-rou(Ic+1-i,j+1))/(rou(i,j)-rou(Ic+1-i,j))
   else if(j==2.and.(i>=Ib1.and.i<=Ib2).or.abs(rou(i,j)-rou(i,j-1))<1d-30) then
   rsp=0
   else if(j==2) then
   rsp=(rou(i,j-1)-rou(Ic+1-i,j-1))/(rou(i,j)-rou(i,j-1))
   else
   rsp=(rou(i,j-1)-rou(i,j-2))/(rou(i,j)-rou(i,j-1))
   end if
   if(j==1.and.(i>=Ib1.and.i<=Ib2).or.abs(rou(i,j)-rou(Ic+1-i,j))<1d-30) then
   rsm=0
   else if(j==1) then
   rsm=(rou(i,j+1)-rou(i,j))/(rou(i,j)-rou(Ic+1-i,j))
   else if(abs(rou(i,j)-rou(i,j-1))<1d-30) then
   rsm=0
   else
   rsm=(rou(i,j+1)-rou(i,j))/(rou(i,j)-rou(i,j-1))
   end if
   if(j==1.and.(i>=Ib1.and.i<=Ib2).or.abs(rou(i,j+1)-rou(i,j))<1d-30) then
   rnp=0
   else if(j==1) then
   rnp=(rou(i,j)-rou(Ic+1-i,j))/(rou(i,j+1)-rou(i,j))
   else
   rnp=(rou(i,j)-rou(i,j-1))/(rou(i,j+1)-rou(i,j))
   end if
   if(j==Jc-1.or.abs(rou(i,j+1)-rou(i,j))<1d-30) then
   rnm=0
   else
   rnm=(rou(i,j+2)-rou(i,j+1))/(rou(i,j+1)-rou(i,j))
   end if
   Psiwp=(rwp+rwp**2)/(1+rwp**2)
   Psiwm=(rwm+rwm**2)/(1+rwm**2)
   Psiep=(rep+rep**2)/(1+rep**2)
   Psiem=(rem+rem**2)/(1+rem**2)
   Psisp=(rsp+rsp**2)/(1+rsp**2)
   Psism=(rsm+rsm**2)/(1+rsm**2)
   Psinp=(rnp+rnp**2)/(1+rnp**2)
   Psinm=(rnm+rnm**2)/(1+rnm**2)
   rouw(i,j)=(0.5-ww(i,j))*(rou(i,j)+Psiwm*(rou(i-1,j)-rou(i,j))/2)+(0.5+ww(i,j))*(rou(i-1,j)+Psiwp*(rou(i,j)-rou(i-1,j))/2)
   roue(i,j)=(0.5+we(i,j))*(rou(i,j)+Psiep*(rou(i+1,j)-rou(i,j))/2)+(0.5-we(i,j))*(rou(i+1,j)+Psiem*(rou(i,j)-rou(i+1,j))/2)
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
   rous(i,j)=(0.5-ws(i,j))*(rou(i,j)+Psism*(rou(Ic+1-i,j)-rou(i,j))/2)+(0.5+ws(i,j))*(rou(Ic+1-i,j)+Psisp*(rou(i,j)-rou(Ic+1-i,j))/2)
   else if(j==1) then
   rous(i,j)=rou(i,j)
   else
   rous(i,j)=(0.5-ws(i,j))*(rou(i,j)+Psism*(rou(i,j-1)-rou(i,j))/2)+(0.5+ws(i,j))*(rou(i,j-1)+Psisp*(rou(i,j)-rou(i,j-1))/2)
   end if
   roun(i,j)=(0.5+wn(i,j))*(rou(i,j)+Psinp*(rou(i,j+1)-rou(i,j))/2)+(0.5-wn(i,j))*(rou(i,j+1)+Psinm*(rou(i,j)-rou(i,j+1))/2)
   end if
   end if
   end DO
  end DO
end Subroutine Densityinterp
