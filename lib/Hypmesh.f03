Subroutine hypmeshgen(Xb,Yb,Xg,Yg,ra,fd,Ip,Jp,trailconfig,gtype)
implicit none
real(8),parameter::Pi=3.1415926535897932d+0
integer i,j,l,Ip,Jp,jtran,Jmax
real(8) fd,dis,nua,nuj,theta,epsic
real(8) Xb(Ip),Yb(Ip),ra(Jp)
real(8) Xg(Ip,Jp),Yg(Ip,Jp)
real(8) axE1,axW1,axP1,ayE1,ayW1,ayP1,b1,axE2,axW2,axP2,ayE2,ayW2,ayP2,b2,df,N,Se,cosa,angle,a,stab1,stab2,&
tx,ty,nx,ny,dsp,dsm,dsp0,dsm0,Xgkm,Ygkm,Xgam,Ygam,Detm,epsi,txplus,txminus,typlus,tyminus,Xgk,Ygk,Xga,Yga,Det
real(8) V(Ip),V0(Ip),Vo(Ip),Xg0(Ip),Yg0(Ip),Xg00(Ip),Yg00(Ip),d(Ip),d0(Ip)
real(8) Ap(2,2,Ip),Aw(2,2,Ip),Ae(2,2,Ip),B(2,Ip),X(2,Ip)
character(*) trailconfig,gtype

Xg(:,1)=Xb
Yg(:,1)=Yb
d=0
d0=0
Jmax=Jp
jtran=int(3*Jmax/4)
dis=fd
nua=2./3
if(trailconfig=='blunt') then
 theta=6
else
 theta=1
end if
epsic=1
DO j=2,Jp
 dis=dis*ra(j)
 nuj=2**(2-j)
 Xg0=Xg(:,j-1)
 Yg0=Yg(:,j-1)
 if(j>=3) then
  Xg00=Xg(:,j-2)
  Yg00=Yg(:,j-2)
  Call checkJacobi(Xg0,Yg0,Xg00,Yg00,Ip)
 else
  Xg00=Xg0
  Yg00=Yg0
 end if
 DO i=1,Ip
  if(i==1) then
   V(i)=sqrt((Xg0(i+1)-Xg0(Ip-1))**2+(Yg0(i+1)-Yg0(Ip-1))**2)*dis/2
  else if(i==Ip) then
   V(i)=sqrt((Xg0(2)-Xg0(i-1))**2+(Yg0(2)-Yg0(i-1))**2)*dis/2
  else
   V(i)=sqrt((Xg0(i+1)-Xg0(i-1))**2+(Yg0(i+1)-Yg0(i-1))**2)*dis/2
  end if
 end DO
 DO l=1,j-2
  Vo=V
  if(gtype=='O') then
   DO i=1,Ip
    if(i==1) then
     V(i)=(1-nua)*Vo(i)+nua*0.5*(Vo(i+1)+Vo(Ip-1))
    else if(i==Ip) then
     V(i)=(1-nua)*Vo(i)+nua*0.5*(Vo(2)+Vo(i-1))
    else
     V(i)=(1-nua)*Vo(i)+nua*0.5*(Vo(i+1)+Vo(i-1))
    end if
   end DO
  else if(gtype=='C') then
   DO i=2,Ip-1
    if(i==2) then
     V(i)=(1-nua)*Vo(i)+nua*0.5*(Vo(i+1)+Vo(i))
    else if(i==Ip-1) then
     V(i)=(1-nua)*Vo(i)+nua*0.5*(Vo(i)+Vo(i-1))
    else
     V(i)=(1-nua)*Vo(i)+nua*0.5*(Vo(i+1)+Vo(i-1))
    end if
   end DO
  end if
 end DO
 if(j==2) then
  V0=V
 end if
 DO i=1,Ip
  if(i==Ip) then
   dsp=sqrt((Xg0(2)-Xg0(i))**2+(Yg0(2)-Yg0(i))**2)
   dsp0=sqrt((Xg00(2)-Xg00(i))**2+(Yg00(2)-Yg00(i))**2)
  else
   dsp=sqrt((Xg0(i+1)-Xg0(i))**2+(Yg0(i+1)-Yg0(i))**2)
   dsp0=sqrt((Xg00(i+1)-Xg00(i))**2+(Yg00(i+1)-Yg00(i))**2)
  end if
  if(i==1) then
   dsm=sqrt((Xg0(Ip-1)-Xg0(i))**2+(Yg0(Ip-1)-Yg0(i))**2)
   dsm0=sqrt((Xg00(Ip-1)-Xg00(i))**2+(Yg00(Ip-1)-Yg00(i))**2)
  else
   dsm=sqrt((Xg0(i-1)-Xg0(i))**2+(Yg0(i-1)-Yg0(i))**2)
   dsm0=sqrt((Xg00(i-1)-Xg00(i))**2+(Yg00(i-1)-Yg00(i))**2)
  end if
  d(i)=(dsp0+dsm0)/(dsp+dsm)
 end DO
 if(j-1>jtran.and.maxval(d)-maxval(d0)<0.and.jtran==int(3*Jmax/4)) then
  jtran=j-1
 end if
 if(j==2) then
  Se=sqrt(1./(Jmax-1))
 else if(j-1<=jtran) then
  Se=sqrt((j-2.)/(Jmax-1))
 else if(j-1>jtran) then
  Se=sqrt((jtran-1.)/(Jmax-1))
 end if
 DO i=1,Ip
  if(i==Ip) then
   dsp=sqrt((Xg0(2)-Xg0(i))**2+(Yg0(2)-Yg0(i))**2)
   txplus=(Xg0(2)-Xg0(i))/dsp
   typlus=(Yg0(2)-Yg0(i))/dsp
  else
   dsp=sqrt((Xg0(i+1)-Xg0(i))**2+(Yg0(i+1)-Yg0(i))**2)
   txplus=(Xg0(i+1)-Xg0(i))/dsp
   typlus=(Yg0(i+1)-Yg0(i))/dsp
  end if
  if(i==1) then
   dsm=sqrt((Xg0(Ip-1)-Xg0(i))**2+(Yg0(Ip-1)-Yg0(i))**2)
   txminus=(Xg0(Ip-1)-Xg0(i))/dsm
   tyminus=(Yg0(Ip-1)-Yg0(i))/dsm
  else
   dsm=sqrt((Xg0(i-1)-Xg0(i))**2+(Yg0(i-1)-Yg0(i))**2)
   txminus=(Xg0(i-1)-Xg0(i))/dsm
   tyminus=(Yg0(i-1)-Yg0(i))/dsm
  end if
  if(i==1) then
   Xgk=(Xg0(i+1)-Xg0(Ip-1))/2
   Ygk=(Yg0(i+1)-Yg0(Ip-1))/2
  else if(i==Ip) then
   Xgk=(Xg0(2)-Xg0(i-1))/2
   Ygk=(Yg0(2)-Yg0(i-1))/2
  else
   Xgk=(Xg0(i+1)-Xg0(i-1))/2
   Ygk=(Yg0(i+1)-Yg0(i-1))/2
  end if
  Det=Xgk**2+Ygk**2
  if(Det==0d+0) then
   Xga=-V0(i)
   Yga=V0(i)
  else
   Xga=-Ygk*V0(i)/Det
   Yga=Xgk*V0(i)/Det
  end if
  tx=txplus-txminus
  ty=typlus-tyminus
  Xgkm=(dsp+dsm)*tx/2
  Ygkm=(dsp+dsm)*ty/2
  Detm=Xgkm**2+Ygkm**2
  if(Detm==0d+0) then
   Xgam=-V0(i)
   Ygam=V0(i)
  else
   Xgam=-Ygkm*V0(i)/Detm
   Ygam=Xgkm*V0(i)/Detm
  end if
  Xga=nuj*Xgam+(1-nuj)*Xga
  Yga=nuj*Ygam+(1-nuj)*Yga
  df=max(d(i)**(2/Se),1d-1)
  if(Det==0d+0) then
   N=0
  else
   N=sqrt((Xga**2+Yga**2)/Det)
  end if
  if(tx==0d+0.and.ty==0d+0) then
   nx=0
   ny=0
  else
   nx=-ty/sqrt(tx**2+ty**2)
   ny=tx/sqrt(tx**2+ty**2)
  end if
  cosa=nx*txplus+ny*typlus
  angle=acos(cosa)
  if(angle>=0.and.angle<=Pi/2) then
   a=1/(1-cosa**2)
  else
   a=1
  end if
  epsi=epsic*N*Se*df*a
  if(i==1) then
   stab1=epsi*((Xg0(i+1)+Xg0(Ip-1)-2*Xg0(i))-(Xg00(i+1)+Xg00(Ip-1)-2*Xg00(i)))
   stab2=epsi*((Yg0(i+1)+Yg0(Ip-1)-2*Yg0(i))-(Yg00(i+1)+Yg00(Ip-1)-2*Yg00(i)))
  else if(i==Ip) then
   stab1=epsi*((Xg0(2)+Xg0(i-1)-2*Xg0(i))-(Xg00(2)+Xg00(i-1)-2*Xg00(i)))
   stab2=epsi*((Yg0(2)+Yg0(i-1)-2*Yg0(i))-(Yg00(2)+Yg00(i-1)-2*Yg00(i)))
  else
   stab1=epsi*((Xg0(i+1)+Xg0(i-1)-2*Xg0(i))-(Xg00(i+1)+Xg00(i-1)-2*Xg00(i)))
   stab2=epsi*((Yg0(i+1)+Yg0(i-1)-2*Yg0(i))-(Yg00(i+1)+Yg00(i-1)-2*Yg00(i)))
  end if
  if(Det==0d+0) then
   ayW1=0
   ayE1=0
   axW1=2*epsi
   axE1=2*epsi
   b1=-stab1
   ayW2=2*epsi
   ayE2=2*epsi
   axW2=0
   axE2=0
   b2=-stab2
  else
   ayW1=(1+theta)*(Xgk*Yga+Ygk*Xga)/Det/2
   ayE1=-(1+theta)*(Xgk*Yga+Ygk*Xga)/Det/2
   axW1=(1+theta)*(Xgk*Xga-Ygk*Yga)/Det/2+2*epsi
   axE1=-(1+theta)*(Xgk*Xga-Ygk*Yga)/Det/2+2*epsi
   b1=-Ygk*V(i)/Det-stab1
   ayW2=(1+theta)*(Ygk*Yga-Xgk*Xga)/Det/2+2*epsi
   ayE2=-(1+theta)*(Ygk*Yga-Xgk*Xga)/Det/2+2*epsi
   axW2=(1+theta)*(Xgk*Yga+Ygk*Xga)/Det/2
   axE2=-(1+theta)*(Xgk*Yga+Ygk*Xga)/Det/2
   b2=Xgk*V(i)/Det-stab2
  end if
  axP1=1+4*epsi
  ayP1=0
  axP2=0
  ayP2=1+4*epsi
  Ap(1,1,i)=axP1
  Aw(1,1,i)=axW1
  Ae(1,1,i)=axE1
  Ap(1,2,i)=ayP1
  Aw(1,2,i)=ayW1
  Ae(1,2,i)=ayE1
  Ap(2,1,i)=axP2
  Aw(2,1,i)=axW2
  Ae(2,1,i)=axE2
  Ap(2,2,i)=ayP2
  Aw(2,2,i)=ayW2
  Ae(2,2,i)=ayE2
  B(1,i)=b1
  B(2,i)=b2
 end DO
 if(gtype=='C') then
  Ap(1,1,1)=1
  Aw(1,1,1)=0
  Ae(1,1,1)=0
  Ap(1,2,1)=0
  Aw(1,2,1)=0
  Ae(1,2,1)=0
  Ap(2,1,1)=0
  Aw(2,1,1)=0
  Ae(2,1,1)=0
  Ap(2,2,1)=1
  Aw(2,2,1)=0
  Ae(2,2,1)=1
  Ap(1,1,Ip)=1
  Aw(1,1,Ip)=0
  Ae(1,1,Ip)=0
  Ap(1,2,Ip)=0
  Aw(1,2,Ip)=0
  Ae(1,2,Ip)=0
  Ap(2,1,Ip)=0
  Aw(2,1,Ip)=0
  Ae(2,1,Ip)=0
  Ap(2,2,Ip)=1
  Aw(2,2,Ip)=1
  Ae(2,2,Ip)=0
  B(1,1)=0
  B(2,1)=0
  B(1,Ip)=0
  B(2,Ip)=0
  Call BlockTDMA(Ap,Aw,Ae,B,X,2,Ip)
 else if(gtype=='O') then
  Call BlockCTDMA(Ap(:,:,1:Ip-1),Aw(:,:,1:Ip-1),Ae(:,:,1:Ip-1),B(:,1:Ip-1),X(:,1:Ip-1),2,Ip-1)
  X(:,Ip)=X(:,1)
 end if
 Xg(:,j)=Xg0+X(1,:)
 Yg(:,j)=Yg0+X(2,:)
 V0=V
 d0=d
end DO

end Subroutine hypmeshgen

Subroutine checkJacobi(Xg,Yg,Xg0,Yg0,Ip)
integer i,Ip,Ic
real(8) Xgk,Ygk,Xga,Yga,Jg
real(8) Xg(Ip),Yg(Ip),Xg0(Ip),Yg0(Ip)

Ic=Ip-1

DO i=1,Ic
  Xgk=((Xg(i+1)+Xg0(i+1))-(Xg(i)+Xg0(i)))/2
  Ygk=((Yg(i+1)+Yg0(i+1))-(Yg(i)+Yg0(i)))/2
  Xga=((Xg(i)+Xg(i+1))-(Xg0(i)+Xg0(i+1)))/2
  Yga=((Yg(i)+Yg(i+1))-(Yg0(i)+Yg0(i+1)))/2
  Jg=Xgk*Yga-Xga*Ygk
  if(Jg<=0) stop 'Error: The Jacobi of generated mesh is less than zero!' 
end DO
  
end Subroutine checkJacobi
