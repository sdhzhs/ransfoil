Subroutine Condiff(scalar)
use Aero2DCOM
implicit none
integer i,j
real(8) Dnow,Dnoe,Dnos,Dnon,Faw,Fae,Fks,Fkn,Fwallw,Fwalle,Dwplus,fai1,F1,betaistar,Fmt,alphaf,betai,fniu2,Sv,rm,gm
real(8) Fwall(Ib1:Ib2)
real(8) F(Ic,Jc),Ga(Ic,Jc)
real(8) Fw(Ic,Jc),Fe(Ic,Jc),Fs(Ic,Jc),Fn(Ic,Jc),Dw(Ic,Jc),De(Ic,Jc),Ds(Ic,Jc),Dn(Ic,Jc),DF(Ic,Jc),bno(Ic,Jc),cor(Ic,Jc)
real(8) betastar(Ic,Jc),alpha(Ic,Jc),beta(Ic,Jc),Dwt(Ic,Jc),Sm(Ic,Jc),fw1(Ic,Jc),C3e(Ic,Jc)
character(*) scalar
if(scalar=='U'.or.scalar=='V') then
Fwall=0
Ga=miu+miut
if(scalar=='U') then
F=U
else if(scalar=='V') then
F=V
end if
else if(scalar=='T') then
F=T
Fwall=Tf
Ga=ka/ca+miut/Prt
else if(scalar=='Tn') then
F=Tn
Fwall=1.5*Tn(Ib1:Ib2,1)-0.5*Tn(Ib1:Ib2,2)
Ga=(miu+rou*Tn)/sigman
else if(scalar=='Tk'.and.Turmod=='ke') then
F=Tk
Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
Ga=miu+miut/sigmak
else if(scalar=='Tk'.and.Turmod=='sst') then
F=Tk
Fwall=1.5*Tk(Ib1:Ib2,1)-0.5*Tk(Ib1:Ib2,2)
Ga=miu+miut/sigmatk
else if(scalar=='Te') then
F=Te
Fwall=1.5*Te(Ib1:Ib2,1)-0.5*Te(Ib1:Ib2,2)
Ga=miu+miut/sigmae
else if(scalar=='Tw') then
F=Tw
Fwall=1.5*Tw(Ib1:Ib2,1)-0.5*Tw(Ib1:Ib2,2)
Ga=miu+miut/sigmatw
end if
DO j=1,Jc-1
  DO i=2,Ic-1
   Dw(i,j)=interpl(a1(i,j)*Ga(i,j),a1(i-1,j)*Ga(i-1,j),dk(i,j),dk(i-1,j))*dy/dx
   De(i,j)=interpl(a1(i,j)*Ga(i,j),a1(i+1,j)*Ga(i+1,j),dk(i,j),dk(i+1,j))*dy/dx
   Dn(i,j)=interpl(y1(i,j)*Ga(i,j),y1(i,j+1)*Ga(i,j+1),da(i,j),da(i,j+1))*dx/dy
   Dnow=interpl(b1(i,j)*Ga(i,j),b1(i-1,j)*Ga(i-1,j),dk(i,j),dk(i-1,j))*dy
   Dnoe=interpl(b1(i,j)*Ga(i,j),b1(i+1,j)*Ga(i+1,j),dk(i,j),dk(i+1,j))*dy
   Dnon=interpl(b1(i,j)*Ga(i,j),b1(i,j+1)*Ga(i,j+1),da(i,j),da(i,j+1))*dx
   if(j==1.and.(i<Ib1.or.i>Ib2)) then
   Ds(i,j)=interpl(y1(i,j)*Ga(i,j),y1(Ic+1-i,j)*Ga(Ic+1-i,j),da(i,j),da(Ic+1-i,j))*dx/dy
   Dnos=interpl(b1(i,j)*Ga(i,j),b1(Ic+1-i,j)*Ga(Ic+1-i,j),da(i,j),da(Ic+1-i,j))*dx
   Faw=(F(i-1,j+1)+F(i,j+1)-F(Ic+2-i,j)-F(Ic+1-i,j))/(4*dy)
   Fae=(F(i+1,j+1)+F(i,j+1)-F(Ic-i,j)-F(Ic+1-i,j))/(4*dy)
   Fks=(F(Ic-i,j)+F(i+1,j)-F(Ic+2-i,j)-F(i-1,j))/(4*dx)
   Fkn=(F(i+1,j+1)+F(i+1,j)-F(i-1,j+1)-F(i-1,j))/(4*dx)
   else if(j==1) then
   Ds(i,j)=y1(i,j)*Ga(i,j)*dx/dy
   Dnos=b1(i,j)*Ga(i,j)*dx
   if(i==Ib1) then
   Fwallw=Fwall(i)
   else
   Fwallw=0.5*(Fwall(i)+Fwall(i-1))
   end if
   if(i==Ib2) then
   Fwalle=Fwall(i)
   else
   Fwalle=0.5*(Fwall(i)+Fwall(i+1))
   end if
   Faw=(0.25*(F(i-1,j+1)+F(i,j+1)+F(i-1,j)+F(i,j))-Fwallw)/dy
   Fae=(0.25*(F(i+1,j+1)+F(i,j+1)+F(i+1,j)+F(i,j))-Fwalle)/dy
   Fks=(Fwalle-Fwallw)/dx
   Fkn=(F(i+1,j+1)+F(i+1,j)-F(i-1,j+1)-F(i-1,j))/(4*dx)
   else
   Ds(i,j)=interpl(y1(i,j)*Ga(i,j),y1(i,j-1)*Ga(i,j-1),da(i,j),da(i,j-1))*dx/dy
   Dnos=interpl(b1(i,j)*Ga(i,j),b1(i,j-1)*Ga(i,j-1),da(i,j),da(i,j-1))*dx
   Faw=(F(i-1,j+1)+F(i,j+1)-F(i-1,j-1)-F(i,j-1))/(4*dy)
   Fae=(F(i+1,j+1)+F(i,j+1)-F(i+1,j-1)-F(i,j-1))/(4*dy)
   Fks=(F(i+1,j-1)+F(i+1,j)-F(i-1,j-1)-F(i-1,j))/(4*dx)
   Fkn=(F(i+1,j+1)+F(i+1,j)-F(i-1,j+1)-F(i-1,j))/(4*dx)
   end if
   bno(i,j)=Dnow*Faw-Dnoe*Fae+Dnos*Fks-Dnon*Fkn
  end DO
end DO
aP=1
aW=0
aE=0
aS=0
aN=0
DO j=1,Jc-1
  DO i=2,Ic-1
    Fw(i,j)=dy*rouk(i,j)*Unk(i,j)
    Fe(i,j)=dy*rouk(i+1,j)*Unk(i+1,j)
    Fs(i,j)=dx*roua(i,j)*Vna(i,j)
    Fn(i,j)=dx*roua(i,j+1)*Vna(i,j+1)
    DF(i,j)=Fe(i,j)-Fw(i,j)+Fn(i,j)-Fs(i,j)
    aW(i,j)=Dw(i,j)+max(Fw(i,j),0.0)
    aE(i,j)=De(i,j)+max(-Fe(i,j),0.0)
    aN(i,j)=Dn(i,j)+max(-Fn(i,j),0.0)
    if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    aS(i,j)=0
    if((Turmod=='sa'.and.Walltreat=='lr').or.(Turmod=='sst'.and.Walltreat=='lr').or.Turmod=='lam'.or.Turmod=='inv') then
     if(scalar=='Tk') then
      aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)
     else if(scalar=='Tn'.and.ks(i)>1e-10) then
      aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)+2*Ds(i,j)*Yp(i)/d(i,j)
     else if(scalar=='Tw') then
      aW(i,j)=0
      aE(i,j)=0
      aN(i,j)=0
     else
      aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)+2*Ds(i,j)
     end if
    else if(Turmod=='sa'.and.Walltreat=='wf'.or.(Turmod=='sst'.and.Walltreat=='wf').or.Turmod=='ke') then
     if(scalar=='U'.or.scalar=='V') then
      aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)+rou(i,j)*ustar(i)*DR(i)/Uplus(i)
     else if(scalar=='T') then
      aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)+rou(i,j)*ustar(i)*DR(i)/Tplus(i)
     else if(scalar=='Tn') then
      aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)+2*Ds(i,j)
     else if(Turmod=='sst'.and.scalar=='Tk') then
      aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)+rou(i,j)*ustar(i)**3*Jg(i,j)*dx*dy/(Tk(i,j)*kapa*Yp(i))
     else if(Turmod=='ke'.and.scalar=='Tk') then
      aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)+rou(i,j)*ustar(i)**3*Jg(i,j)*dx*dy/(Tk(i,j)*kapa*Yp(i))
     else if(scalar=='Te'.or.scalar=='Tw') then
      aW(i,j)=0
      aE(i,j)=0
      aN(i,j)=0
     end if
    end if
    else
    aS(i,j)=Ds(i,j)+max(Fs(i,j),0.0)
    aP(i,j)=aW(i,j)+aE(i,j)+aS(i,j)+aN(i,j)+DF(i,j)
    end if
  end DO
end DO
Call Defercorrect(F,cor,Fw,Fe,Fs,Fn)
if(scalar=='T'.or.scalar=='Tk'.or.scalar=='Tw'.or.scalar=='Te') then
St=sqrt(2*(Ux**2+Vy**2)+(Uy+Vx)**2)
end if
if(Turmod=='sst'.and.(scalar=='Tk'.or.scalar=='Tw')) then
 DO j=1,Jc-1
   DO i=2,Ic-1
   Dwplus=max(2*rou(i,j)*(Tkx(i,j)*Twx(i,j)+Tky(i,j)*Twy(i,j))/(sigmaw2*Tw(i,j)),1e-10)
   fai1=min(max(sqrt(Tk(i,j))/(0.09*Tw(i,j)*d(i,j)),500*miu(i,j)/(rou(i,j)*d(i,j)**2*Tw(i,j))),&
   4*rou(i,j)*Tk(i,j)/(sigmaw2*Dwplus*d(i,j)**2))
   F1=tanh(fai1**4)
   if(sqrt(2*Tk(i,j)/(gama*R*T(i,j)/Ma))<=Mt0) then
   Fmt=0
   else
   Fmt=2*Tk(i,j)/(gama*R*T(i,j)/Ma)-Mt0**2
   end if
   betaistar=betastarf*(4./15+(Ret(i,j)/Rbeta)**4)/(1+(Ret(i,j)/Rbeta)**4)
   betastar(i,j)=betaistar*(1+Zetastar*Fmt)
   alphaf=F1*(betai1/betastarf-kapa**2/(sigmaw1*sqrt(betastarf)))+&
   (1-F1)*(betai2/betastarf-kapa**2/(sigmaw2*sqrt(betastarf)))
   alpha(i,j)=alphaf/alphastar(i,j)*(alpha0+Ret(i,j)/Rw)/(1+Ret(i,j)/Rw)
   betai=F1*betai1+(1-F1)*betai2
   beta(i,j)=betai*(1-betaistar/betai*Zetastar*Fmt)
   Dwt(i,j)=2*(1-F1)*rou(i,j)*sigmaw2*(Tkx(i,j)*Twx(i,j)+Tky(i,j)*Twy(i,j))/Tw(i,j)
   end DO
 end DO
else if(scalar=='Tn') then
 DO j=1,Jc-1
   DO i=2,Ic-1
   fniu2=1-Tn(i,j)/(miu(i,j)/rou(i,j)+Tn(i,j)*fniu1(i,j))
   Sv=abs(Uy(i,j)-Vx(i,j))
   Sm(i,j)=Sv+Tn(i,j)*fniu2/(kapa*d(i,j))**2
   rm=Tn(i,j)/(Sm(i,j)*(kapa*d(i,j))**2)
   gm=rm+Cw2*(rm**6-rm)
   fw1(i,j)=gm*((1+Cw3**6)/(gm**6+Cw3**6))**(1./6)
   end DO
 end DO
else if(scalar=='Te') then
 C3e=tanh(abs(V/U))
end if
b=F
DO j=1,Jc-1
  DO i=2,Ic-1
  if(scalar=='U') then
   if(Proctrl=='com') then
    b(i,j)=Jg(i,j)*(-Px(i,j)-2*(muxx(i,j)+mvyx(i,j))/3+muxx(i,j)+mvxy(i,j))*dx*dy
   else if(Proctrl=='incom') then
    b(i,j)=Jg(i,j)*(-Px(i,j)+muxx(i,j)+mvxy(i,j))*dx*dy
   end if
  else if(scalar=='V') then
   if(Proctrl=='com') then
    b(i,j)=Jg(i,j)*(-Py(i,j)-2*(muxy(i,j)+mvyy(i,j))/3+muyx(i,j)+mvyy(i,j))*dx*dy
   else if(Proctrl=='incom') then
    b(i,j)=Jg(i,j)*(-Py(i,j)+muyx(i,j)+mvyy(i,j))*dx*dy
   end if
  else if(scalar=='T') then
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    if((Turmod=='sa'.and.Walltreat=='lr').or.(Turmod=='sst'.and.Walltreat=='lr').or.Turmod=='lam'.or.Turmod=='inv') then
     if(Proctrl=='com') then
      if(visheat=='Y') then
      b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j)-2*(miu(i,j)+miut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+&
      (miu(i,j)+miut(i,j))*St(i,j)**2)*dx*dy/ca+2*Ds(i,j)*Tf
      else if(visheat=='N') then
      b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))*dx*dy/ca+2*Ds(i,j)*Tf
      end if
     else if(Proctrl=='incom') then
      if(visheat=='Y') then
      b(i,j)=Jg(i,j)*(miu(i,j)+miut(i,j))*St(i,j)**2*dx*dy/ca+2*Ds(i,j)*Tf
      else if(visheat=='N') then
      b(i,j)=2*Ds(i,j)*Tf
      end if
     end if
    else if(Turmod=='sst'.and.Walltreat=='wf'.or.(Turmod=='sa'.and.Walltreat=='wf').or.Turmod=='ke') then
     if(Proctrl=='com') then
      b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))*dx*dy/ca+rou(i,j)*ustar(i)*Tf*DR(i)/Tplus(i)
     else if(Proctrl=='incom') then
      b(i,j)=rou(i,j)*ustar(i)*Tf*DR(i)/Tplus(i)
     end if
    end if
   else
    if(Proctrl=='com') then
     if(visheat=='Y') then
     b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j)-2*(miu(i,j)+miut(i,j))*(Ux(i,j)+Vy(i,j))**2/3+&
     (miu(i,j)+miut(i,j))*St(i,j)**2)*dx*dy/ca
     else if(visheat=='N') then
     b(i,j)=Jg(i,j)*(U(i,j)*Px(i,j)+V(i,j)*Py(i,j))*dx*dy/ca
     end if
    else if(Proctrl=='incom') then
     if(visheat=='Y') then
     b(i,j)=Jg(i,j)*(miu(i,j)+miut(i,j))*St(i,j)**2*dx*dy/ca
     else if(visheat=='N') then
     b(i,j)=0
     end if
    end if
   end if
  else if(scalar=='Tn') then
   b(i,j)=Cb2*rou(i,j)*(Tnx(i,j)**2+Tny(i,j)**2)*Jg(i,j)*dx*dy/sigman+rou(i,j)*Cb1*Sm(i,j)*Tn(i,j)*Jg(i,j)*dx*dy-&
   rou(i,j)*Cw1*fw1(i,j)*(Tn(i,j)/d(i,j))**2*Jg(i,j)*dx*dy
  else if(scalar=='Tk'.and.Turmod=='ke') then
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    if(Proctrl=='com') then
    b(i,j)=g*miut(i,j)*rouy(i,j)*Jg(i,j)*dx*dy/(rou(i,j)*Prt)-2*rou(i,j)*Te(i,j)*Tk(i,j)*Jg(i,j)*dx*dy/(gama*R*T(i,j)/Ma)+&
    rou(i,j)*ustar(i)**2*sqrt(U(i,j)**2+V(i,j)**2)*Jg(i,j)*dx*dy/(Uplus(i)*kapa*Yp(i))
    else if(Proctrl=='incom') then
    b(i,j)=rou(i,j)*ustar(i)**2*sqrt(U(i,j)**2+V(i,j)**2)*Jg(i,j)*dx*dy/(Uplus(i)*kapa*Yp(i))
    end if
   else
    if(Proctrl=='com') then
    b(i,j)=(miu(i,j)+miut(i,j))*St(i,j)**2*Jg(i,j)*dx*dy+g*miut(i,j)*rouy(i,j)*Jg(i,j)*dx*dy/(rou(i,j)*Prt)-&
    2*rou(i,j)*Te(i,j)*Tk(i,j)*Jg(i,j)*dx*dy/(gama*R*T(i,j)/Ma)-rou(i,j)*Te(i,j)*Jg(i,j)*dx*dy
    else if(Proctrl=='incom') then
    b(i,j)=(miu(i,j)+miut(i,j))*St(i,j)**2*Jg(i,j)*dx*dy-rou(i,j)*Te(i,j)*Jg(i,j)*dx*dy
    end if
   end if
  else if(scalar=='Tk'.and.Turmod=='sst') then
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    if(Walltreat=='wf') then
      b(i,j)=rou(i,j)*ustar(i)**2*sqrt(U(i,j)**2+V(i,j)**2)*Jg(i,j)*dx*dy/(Uplus(i)*kapa*Yp(i))
    else if(Walltreat=='lr') then
      b(i,j)=min(miut(i,j)*St(i,j)**2,10*rou(i,j)*betastar(i,j)*Tk(i,j)*Tw(i,j))*Jg(i,j)*dx*dy-&
      rou(i,j)*betastar(i,j)*Tk(i,j)*Tw(i,j)*Jg(i,j)*dx*dy
    end if
   else
    b(i,j)=min(miut(i,j)*St(i,j)**2,10*rou(i,j)*betastar(i,j)*Tk(i,j)*Tw(i,j))*Jg(i,j)*dx*dy-&
    rou(i,j)*betastar(i,j)*Tk(i,j)*Tw(i,j)*Jg(i,j)*dx*dy
   end if
  else if(scalar=='Te') then
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    cycle
   else
    if(Proctrl=='com') then
    b(i,j)=Te(i,j)*(C1e*((miu(i,j)+miut(i,j))*St(i,j)**2+C3e(i,j)*g*miut(i,j)*rouy(i,j)/(rou(i,j)*Prt))-&
    C2e*rou(i,j)*Te(i,j))*Jg(i,j)*dx*dy/Tk(i,j)
    else if(Proctrl=='incom') then
    b(i,j)=Te(i,j)*(C1e*(miu(i,j)+miut(i,j))*St(i,j)**2-C2e*rou(i,j)*Te(i,j))*Jg(i,j)*dx*dy/Tk(i,j)
    end if
   end if
  else if(scalar=='Tw') then
   if(j==1.and.(i>=Ib1.and.i<=Ib2)) then
    cycle
   else
    b(i,j)=rou(i,j)*alpha(i,j)*St(i,j)**2*Jg(i,j)*dx*dy-rou(i,j)*beta(i,j)*Tw(i,j)**2*Jg(i,j)*dx*dy+Dwt(i,j)*Jg(i,j)*dx*dy
   end if
  end if
  if(scalar=='Tk'.or.scalar=='Te'.or.scalar=='Tw') then
   b(i,j)=b(i,j)+bno(i,j)
  else
   b(i,j)=b(i,j)+bno(i,j)+cor(i,j)
  end if
  end DO
end DO
if(scalar=='U') then
auP=aP
auNB=aW+aE+aS+aN
end if
end Subroutine Condiff
