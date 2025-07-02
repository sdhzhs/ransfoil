Subroutine Cell2node(Fg,Fc,scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Fg(Ip,Jp),Fc(Ic,Jc)
character(*) scalar
DO j=1,Jc
  DO i=Is,Ie+1
  if(j==1.and.i>=Ib1.and.i<=Ib2+1) then
   if((scalar=='U'.or.scalar=='V').and.Turmod/='inv'.or.scalar=='Tn'.or.scalar=='Tk') then
    Fg(i,j)=0
   else if(scalar=='T') then
    Fg(i,j)=Tf
   else if(scalar=='rho'.and.Proctrl=='incom') then
    Fg(i,j)=rhoi
   else
    if(i==1) then
     Fg(i,j)=0.5*(Fc(i,j)+Fc(Ic,j))
    else if(i==Ip) then
     Fg(i,j)=0.5*(Fc(1,j)+Fc(i-1,j))
    else
     Fg(i,j)=0.5*(Fc(i,j)+Fc(i-1,j))
    end if
   end if
  else if(j==1) then
   Fg(i,j)=(Fc(i,j)+Fc(i-1,j)+Fc(Ic+1-i,j)+Fc(Ic+2-i,j))/4
  else
   if(i==1) then
    Fg(i,j)=(Fc(i,j)+Fc(Ic,j)+Fc(i,j-1)+Fc(Ic,j-1))/4
   else if(i==Ip) then
    Fg(i,j)=(Fc(1,j)+Fc(i-1,j)+Fc(1,j-1)+Fc(i-1,j-1))/4
   else
    Fg(i,j)=(Fc(i,j)+Fc(i-1,j)+Fc(i,j-1)+Fc(i-1,j-1))/4
   end if
  end if
  end DO
end DO
if(scalar=='U') then
 Fg(:,Jp)=Ui
 if(Is>1) then 
  Fg(1,:)=Ui
  Fg(Ip,:)=Ui
 end if
else if(scalar=='V') then
 Fg(:,Jp)=Vi
 if(Is>1) then
  Fg(1,:)=Vi
  Fg(Ip,:)=Vi
 end if
else if(scalar=='P') then
 Fg(:,Jp)=0
 if(Is>1) then
  Fg(1,:)=0
  Fg(Ip,:)=0
 end if
else if(scalar=='T') then
 Fg(:,Jp)=Ta
 if(Is>1) then
  Fg(1,:)=Ta
  Fg(Ip,:)=Ta
 end if
else if(scalar=='rho') then
 Fg(:,Jp)=rhoi
 if(Is>1) then
  Fg(1,:)=rhoi
  Fg(Ip,:)=rhoi
 end if
else if(scalar=='Tn') then
 Fg(:,Jp)=Tni
 if(Is>1) then
  Fg(1,:)=Tni
  Fg(Ip,:)=Tni
 end if
else if(scalar=='Tk') then
 Fg(:,Jp)=Tki
 if(Is>1) then
  Fg(1,:)=Tki
  Fg(Ip,:)=Tki
 end if
else if(scalar=='Te') then
 Fg(:,Jp)=Tei
 if(Is>1) then
  Fg(1,:)=Tei
  Fg(Ip,:)=Tei
 end if
else if(scalar=='Tw') then
 Fg(:,Jp)=Twi
 if(Is>1) then
  Fg(1,:)=Twi
  Fg(Ip,:)=Twi
 end if
end if
end Subroutine Cell2node
