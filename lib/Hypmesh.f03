Subroutine hypmeshgen(Xb,Yb,Xg,Yg,ra,fd,Ip,Jp)
implicit none
real(8),parameter::Pi=3.1415926535897932d+0
integer i,j,l,Ip,Jp,jtran,Jmax
real(8) fd,dis,niua,niuj,theta,epsic
real(8) Xb(Ip),Yb(Ip),ra(Jp)
real(8) Xg(Ip,Jp),Yg(Ip,Jp)
real(8) axE1,axW1,axP1,ayE1,ayW1,ayP1,b1,axE2,axW2,axP2,ayE2,ayW2,ayP2,b2,Se,a,cosa,df,txplus,txminus,typlus,tyminus,tx,ty,nx,ny,&
dsp,dsm,dsp0,dsm0,Xgk,Ygk,Xga,Yga,Det,Xgkm,Ygkm,Xgam,Ygam,Detm,epsi,stab1,stab2,N
real(8) V(Ip),V0(Ip),Xg0(Ip),Yg0(Ip),Xg00(Ip),Yg00(Ip),d(Ip),d0(Ip)
real(8) Ap(2,2,Ip),Aw(2,2,Ip),Ae(2,2,Ip),B(2,Ip),X(2,Ip)

Xg(:,1)=Xb
Yg(:,1)=Yb
Jmax=Jp
jtran=int(3*Jmax/4)
dis=fd
niua=2./3
theta=1
epsic=10
DO j=2,Jp
 dis=dis*ra(j)
 niuj=2**(2-j)
 Xg0=Xg(:,j-1)
 Yg0=Yg(:,j-1)
 if(j>=3) then
 Xg00=Xg(:,j-2)
 Yg00=Yg(:,j-2)
 else
 Xg00=Xg0
 Yg00=Yg0
 end if
 DO i=2,Ip-1
 V(i)=sqrt((Xg0(i+1)-Xg0(i-1))**2+(Yg0(i+1)-Yg0(i-1))**2)*dis/2
 end DO
 DO l=1,j-2
 DO i=2,Ip-1
 if(i==2) then
 V(i)=(1-niua)*V(i)+niua*0.5*(V(i+1)+V(i))
 else if(i==Ip-1) then
 V(i)=(1-niua)*V(i)+niua*0.5*(V(i)+V(i-1))
 else
 V(i)=(1-niua)*V(i)+niua*0.5*(V(i+1)+V(i-1))
 end if
 end DO
 end DO
 if(j==2) then
 V0=V
 end if
 DO i=2,Ip-1
 dsp=sqrt((Xg0(i+1)-Xg0(i))**2+(Yg0(i+1)-Yg0(i))**2)
 dsm=sqrt((Xg0(i-1)-Xg0(i))**2+(Yg0(i-1)-Yg0(i))**2)
 dsp0=sqrt((Xg00(i+1)-Xg00(i))**2+(Yg00(i+1)-Yg00(i))**2)
 dsm0=sqrt((Xg00(i-1)-Xg00(i))**2+(Yg00(i-1)-Yg00(i))**2)
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
 DO i=2,Ip-1
 dsp=sqrt((Xg0(i+1)-Xg0(i))**2+(Yg0(i+1)-Yg0(i))**2)
 dsm=sqrt((Xg0(i-1)-Xg0(i))**2+(Yg0(i-1)-Yg0(i))**2)
 Xgk=(Xg0(i+1)-Xg0(i-1))/2
 Ygk=(Yg0(i+1)-Yg0(i-1))/2
 Det=Xgk**2+Ygk**2
 Xga=-Ygk*V0(i)/Det
 Yga=Xgk*V0(i)/Det
 txplus=(Xg0(i+1)-Xg0(i))/dsp
 typlus=(Yg0(i+1)-Yg0(i))/dsp
 txminus=(Xg0(i-1)-Xg0(i))/dsm
 tyminus=(Yg0(i-1)-Yg0(i))/dsm
 tx=txplus-txminus
 ty=typlus-tyminus
 Xgkm=(dsp+dsm)*tx/2
 Ygkm=(dsp+dsm)*ty/2
 Detm=Xgkm**2+Ygkm**2
 Xgam=-Ygkm*V0(i)/Detm
 Ygam=Xgkm*V0(i)/Detm
 Xga=niuj*Xgam+(1-niuj)*Xga
 Yga=niuj*Ygam+(1-niuj)*Yga
 df=max(d(i)**(2/Se),0.1)
 N=sqrt((Xga**2+Yga**2)/Det)
 nx=-ty/sqrt(tx**2+ty**2)
 ny=tx/sqrt(tx**2+ty**2)
 cosa=nx*txplus+ny*typlus
 if(acos(cosa)>=0.and.acos(cosa)<=Pi/2) then
 a=1/(1-cosa**2)
 else
 a=1
 end if
 epsi=epsic*N*Se*df*a
 stab1=epsi*((Xg0(i+1)+Xg0(i-1)-2*Xg0(i))-(Xg00(i+1)+Xg00(i-1)-2*Xg00(i)))
 stab2=epsi*((Yg0(i+1)+Yg0(i-1)-2*Yg0(i))-(Yg00(i+1)+Yg00(i-1)-2*Yg00(i)))
 ayW1=(1+theta)*(Xgk*Yga+Ygk*Xga)/Det/2
 ayE1=-(1+theta)*(Xgk*Yga+Ygk*Xga)/Det/2
 axW1=(1+theta)*(Xgk*Xga-Ygk*Yga)/Det/2+2*epsi
 axE1=-(1+theta)*(Xgk*Xga-Ygk*Yga)/Det/2+2*epsi
 axP1=1+4*epsi
 ayP1=0
 b1=-Ygk*V(i)/Det-stab1
 ayW2=(1+theta)*(Ygk*Yga-Xgk*Xga)/Det/2+2*epsi
 ayE2=-(1+theta)*(Ygk*Yga-Xgk*Xga)/Det/2+2*epsi
 axW2=(1+theta)*(Xgk*Yga+Ygk*Xga)/Det/2
 axE2=-(1+theta)*(Xgk*Yga+Ygk*Xga)/Det/2
 axP2=0
 ayP2=1+4*epsi
 b2=Xgk*V(i)/Det-stab2
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
 Xg(:,j)=Xg0+X(1,:)
 Yg(:,j)=Yg0+X(2,:)
 V0=V
 d0=d
end DO

end Subroutine hypmeshgen
