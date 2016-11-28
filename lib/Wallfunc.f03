Subroutine Wallfunc
use Aero2DCOM
implicit none
integer i,j
real(8) Ym,Yt
real(8) Dwplus,fai1,F1
real(8) ksplus(Ib1:Ib2),deltaB(Ib1:Ib2),Prough(Ib1:Ib2),lamda(Ib1:Ib2),Uplusl(Ib1:Ib2),Uplust(Ib1:Ib2),Tplusl(Ib1:Ib2),&
Tplust(Ib1:Ib2),Twplusl(Ib1:Ib2),Twplust(Ib1:Ib2),Twplus(Ib1:Ib2),betai(Ib1:Ib2)
Ym=11.225
Yt=11.8
DO i=Ib1,Ib2
  if(Turmod=='sst') then
  ustar(i)=((miu(i,1)*sqrt(U(i,1)**2+V(i,1)**2)/(rou(i,1)*Yp(i)))**2+(alpha1*Tk(i,1))**2)**0.25
  else if(Turmod=='sa') then
  ustar(i)=((miu(i,1)*sqrt(U(i,1)**2+V(i,1)**2)/(rou(i,1)*Yp(i)))**2+(miut(i,1)/(rou(i,1)*kapa*Yp(i)))**4)**0.25
  else if(Turmod=='ke') then
  ustar(i)=Cu**0.25*Tk(i,1)**0.5
  end if
  ksplus(i)=rou(i,1)*ks(i)*ustar(i)/miu(i,1)
  if(Turmod=='sa'.or.Turmod=='sst') then
  Yplus(i)=rou(i,1)*Yp(i)*ustar(i)/miu(i,1)
  if(Walltreat=='wf'.and.ks(i)>1e-10) then
  Yplus(i)=max(Yplus(i),ksplus(i)/2,2.5)
  end if
  else if(Turmod=='ke') then
  Ystar(i)=rou(i,1)*Yp(i)*ustar(i)/miu(i,1)
  Ystar(i)=max(Ystar(i),ksplus(i)/2,Ym)
  end if
  if(ksplus(i)<=2.25) then
  deltaB(i)=0
  else if(ksplus(i)>2.25.and.ksplus(i)<=90) then
  deltaB(i)=log((ksplus(i)-2.25)/87.25+Cks*ksplus(i))*sin(0.4258*(log(ksplus(i))-0.811))/kapa
  else
  deltaB(i)=log(1+Cks*ksplus(i))/kapa
  end if
  if(ks(i)<1e-10) then
  Prough(i)=Pc(i,1)
  else
  Prough(i)=3.15*Pr(i,1)**0.695*(exp(kapa*deltaB(i))/Ep-1/Ep)**0.359+(1/exp(kapa*deltaB(i)))**0.6*Pc(i,1)
  end if
end DO
if((Turmod=='sa'.or.Turmod=='sst').and.Walltreat=='wf') then
 lamda=-0.01*Yplus**4/(1+5.*Yplus)
 Uplusl=Yplus
 Uplust=log(Ep*Yplus)/kapa-deltaB
 Uplus=exp(lamda)*Uplusl+exp(1./lamda)*Uplust
 lamda=-0.01*(Pr(Ib1:Ib2,1)*Yplus)**4/(1+5.*Pr(Ib1:Ib2,1)**3*Yplus)
 DO i=Ib1,Ib2
  if(visheat=='Y') then
  Tplusl(i)=Pr(i,1)*Uplusl(i)
  Tplust(i)=Prt*(Uplust(i)+Prough(i))
  Tplus(i)=exp(lamda(i))*Tplusl(i)+exp(1./lamda(i))*Tplust(i)
  Q(i)=(ca*rou(i,1)*ustar(i)*(Tf-T(i,1))-exp(lamda(i))*0.5*rou(i,1)*ustar(i)*Pr(i,1)*(U(i,1)**2+V(i,1)**2)-&
  exp(1./lamda(i))*0.5*rou(i,1)*ustar(i)*(Prt*(U(i,1)**2+V(i,1)**2)+(Pr(i,1)-Prt)*(Yt*ustar(i))**2))/Tplus(i)
  Tplusl(i)=Pr(i,1)*Uplusl(i)+0.5*rou(i,1)*ustar(i)*Pr(i,1)*(U(i,1)**2+V(i,1)**2)/Q(i)
  Tplust(i)=Prt*(Uplust(i)+Prough(i))+0.5*rou(i,1)*ustar(i)*(Prt*(U(i,1)**2+V(i,1)**2)+(Pr(i,1)-Prt)*(Yt*ustar(i))**2)/Q(i)
  Tplus(i)=exp(lamda(i))*Tplusl(i)+exp(1./lamda(i))*Tplust(i)
  else
  Tplusl(i)=Pr(i,1)*Uplusl(i)
  Tplust(i)=Prt*(Uplust(i)+Prough(i))
  Tplus(i)=exp(lamda(i))*Tplusl(i)+exp(1./lamda(i))*Tplust(i)
  Q(i)=ca*rou(i,1)*ustar(i)*(Tf-T(i,1))/Tplus(i)
  end if
 end DO
else if(Turmod=='ke') then
DO i=Ib1,Ib2
 Uplus(i)=log(Ep*Ystar(i))/kapa-deltaB(i)
 if(visheat=='Y') then
 Tplus(i)=Prt*(Uplus(i)+Prough(i))
 Q(i)=(ca*rou(i,1)*ustar(i)*(Tf-T(i,1))-0.5*rou(i,1)*ustar(i)*(Prt*(U(i,1)**2+V(i,1)**2)+&
 (Pr(i,1)-Prt)*(Yt*ustar(i))**2))/Tplus(i)
 Tplus(i)=Prt*(Uplus(i)+Prough(i))+0.5*rou(i,1)*ustar(i)*(Prt*(U(i,1)**2+V(i,1)**2)+(Pr(i,1)-Prt)*(Yt*ustar(i))**2)/Q(i)
 else
 Tplus(i)=Prt*(Uplus(i)+Prough(i))
 Q(i)=ca*rou(i,1)*ustar(i)*(Tf-T(i,1))/Tplus(i)
 end if
end DO
end if
Tplus=max(Tplus,0.1)
if(Turmod=='sst') then
DO j=1,Jc-1
  DO i=2,Ic-1
   Dwplus=max(2*rou(i,j)*(Tkx(i,j)*Twx(i,j)+Tky(i,j)*Twy(i,j))/(sigmaw2*Tw(i,j)),1e-10)
   fai1=min(max(sqrt(Tk(i,j))/(0.09*Tw(i,j)*d(i,j)),500*miu(i,j)/(rou(i,j)*d(i,j)**2*Tw(i,j))),&
   4*rou(i,j)*Tk(i,j)/(sigmaw2*Dwplus*d(i,j)**2))
   F1=tanh(fai1**4)
   sigmatk(i,j)=1/(F1/sigmak1+(1-F1)/sigmak2)
   sigmatw(i,j)=1/(F1/sigmaw1+(1-F1)/sigmaw2)
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
   betai(i)=F1*betai1+(1-F1)*betai2
   end if
  end DO
end DO
sigmatk(1,:)=sigmatk(2,:)
sigmatk(Ic,:)=sigmatk(Ic-1,:)
sigmatk(:,Jc)=sigmatk(:,Jc-1)
sigmatw(1,:)=sigmatw(2,:)
sigmatw(Ic,:)=sigmatw(Ic-1,:)
sigmatw(:,Jc)=sigmatw(:,Jc-1)
if(Walltreat=='wf') then
lamda=-0.01*Yplus**4/(1+5.*Yplus)
Twplusl=6./(betai*Yplus**2)
Twplust=1/(sqrt(betastarf)*kapa*Yplus)
Twplus=exp(lamda)*Twplusl+exp(1./lamda)*Twplust
else if(Walltreat=='lr') then
DO i=Ib1,Ib2
 if(ksplus(i)>1e-10) then
 ksplus(i)=max(1.0,ksplus(i))
 if(ksplus(i)<25) then
 Twplusl(i)=(50/ksplus(i))**2
 else if(ksplus(i)>=25) then
 Twplusl(i)=100/ksplus(i)
 end if
 Twplus(i)=min(Twplusl(i),6./(betai(i)*Yplus(i)**2))
 else
 Twplus(i)=6./(betai(i)*Yplus(i)**2)
 end if
end DO
end if
Tw(Ib1:Ib2,1)=rou(Ib1:Ib2,1)*ustar**2*Twplus/miu(Ib1:Ib2,1)
else if(Turmod=='ke') then
Te(Ib1:Ib2,1)=ustar**3/(kapa*Yp)
end if
end Subroutine Wallfunc
