Subroutine Cell2node(Fg,Fc,scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Fg(Ip,Jp),Fc(Ic,Jc)
character(*) scalar
DO j=1,Jc
  DO i=2,Ic
  if(j==1.and.i>=Ib1.and.i<=Ib2+1) then
  if((scalar=='U'.or.scalar=='V').and.Turmod/='inv'.or.scalar=='Tn'.or.scalar=='Tk') then
  Fg(i,j)=0
  else if(scalar=='T') then
  Fg(i,j)=Tf
  else if(scalar=='rou'.and.Proctrl=='incom') then
  Fg(i,j)=roui
  else
  Fg(i,j)=0.5*(Fc(i,j)+Fc(i-1,j))
  end if
  else if(j==1) then
  Fg(i,j)=(Fc(i,j)+Fc(i-1,j)+Fc(Ic+1-i,j)+Fc(Ic+2-i,j))/4
  else
  Fg(i,j)=(Fc(i,j)+Fc(i-1,j)+Fc(i,j-1)+Fc(i-1,j-1))/4
  end if
  end DO
end DO
if(scalar=='U') then
  Fg(:,Jp)=Ui
  Fg(1,:)=Ui
  Fg(Ip,:)=Ui
else if(scalar=='V') then
  Fg(:,Jp)=Vi
  Fg(1,:)=Vi
  Fg(Ip,:)=Vi
else if(scalar=='P') then
  Fg(:,Jp)=0
  Fg(1,:)=0
  Fg(Ip,:)=0
else if(scalar=='T') then
  Fg(:,Jp)=Ta
  Fg(1,:)=Ta
  Fg(Ip,:)=Ta
else if(scalar=='rou') then
  Fg(:,Jp)=roui
  Fg(1,:)=roui
  Fg(Ip,:)=roui
else if(scalar=='Tn') then
  Fg(:,Jp)=Tni
  Fg(1,:)=Tni
  Fg(Ip,:)=Tni
else if(scalar=='Tk') then
  Fg(:,Jp)=Tki
  Fg(1,:)=Tki
  Fg(Ip,:)=Tki
else if(scalar=='Te') then
  Fg(:,Jp)=Tei
  Fg(1,:)=Tei
  Fg(Ip,:)=Tei
else if(scalar=='Tw') then
  Fg(:,Jp)=Twi
  Fg(1,:)=Twi
  Fg(Ip,:)=Twi
end if
end Subroutine Cell2node
