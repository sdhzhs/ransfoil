Subroutine Genmesh(libmod)
use Aero2DCOM
implicit none
integer i,j,l,Iw1,Iw2,Iw3,maxl
real(8) err,fb,ltrail,lfar,ratio,ratio0
real(8),allocatable,dimension(:)::Xt,fac
character(*) libmod
if(libmod=='S'.or.libmod=='I') then
open(unit=1,file=filename(1),status='old')
read(1,*) Iwd
allocate(Xwd(Iwd),Ywd(Iwd))
DO i=1,Iwd
read(1,*) Xwd(i),Ywd(i)
end DO
read(1,*) Iwu
allocate(Xwu(Iwu),Ywu(Iwu))
DO i=1,Iwu
read(1,*) Xwu(i),Ywu(i)
end DO
close(1)
print *,'Read airfoil coordinates completed!'
else if(libmod=='C') then
allocate(Xwd(Iwd),Ywd(Iwd))
allocate(Xwu(Iwu),Ywu(Iwu))
Call C_F_POINTER(cXwd,fXwd,(/Iwd/))
Call C_F_POINTER(cYwd,fYwd,(/Iwd/))
Call C_F_POINTER(cXwu,fXwu,(/Iwu/))
Call C_F_POINTER(cYwu,fYwu,(/Iwu/))
Xwd=fXwd
Ywd=fYwd
Xwu=fXwu
Ywu=fYwu
end if
Iw1=max(Iwd,Iwu)
Iw2=Iw1+Iwd-1
Iw3=Iw2+Iwu-1
if(libmod=='S'.or.libmod=='I') then
Ip=Iw3+Iw1-1
Ic=Ip-1
Jc=Jp-1
Ib1=Iw1
Ib2=Iw3-1
Call Allocarray(libmod)
end if
allocate(Xt(Iw1),fac(Jp))
Xg(Iw2:Iw1:-1,1)=Xwd
Yg(Iw2:Iw1:-1,1)=Ywd
Xg(Iw2:Iw3,1)=Xwu
Yg(Iw2:Iw3,1)=Ywu
fb=sqrt((Xwd(Iwd)-Xwd(Iwd-1))**2+(Ywd(Iwd)-Ywd(Iwd-1))**2)
ltrail=11
Call tanhgridline(ltrail,fb,Xt,Iw1)
DO i=Iw3+1,Ip
j=i-Iw3+1
Xg(i,1)=Xg(Iw3,1)+Xt(j)
Yg(i,1)=Yg(Iw3,1)
end DO
Xg(1:Iw1-1,1)=Xg(Ip:Iw3+1:-1,1)
Yg(1:Iw1-1,1)=Yg(Ip:Iw3+1:-1,1)
lfar=10
ratio=2.0
maxl=100
err=1e-10
DO l=1,maxl
ratio0=ratio
ratio=ratio-(1-ratio**(Jp-Ifd)-(lfar/fd-Ifd+1)*(1-ratio))/(-(Jp-Ifd)*ratio**(Jp-Ifd-1)+lfar/fd-Ifd+1)
if(abs(ratio-ratio0)<err) exit
end DO
DO j=2,Jp
if(j<=Ifd) then
 fac(j)=1.0
else
 fac(j)=ratio
end if
end DO
Call hypmeshgen(Xg(:,1),Yg(:,1),Xg,Yg,fac,fd,Ip,Jp)
print *,'Generate 2D C-type structured mesh completed!'
deallocate(Xt,fac)
end Subroutine Genmesh

Subroutine Allocarray(libmod)
use Aero2DCOM
implicit none
integer i
character(*) libmod
if(libmod=='C') then
DO i=1,8
if(cTurmod(i)(1:1)/=C_NULL_CHAR) then
Turmod(i:i)=cTurmod(i)(1:1)
else
Turmod(i:i)=' '
end if
end DO
end if
allocate(U0(Ic,Jc),V0(Ic,Jc),T0(Ic,Jc),U(Ic,Jc),V(Ic,Jc),T(Ic,Jc),rou(Ic,Jc),miu(Ic,Jc),P(Ic,Jc),dP(Ic,Jc),miut(Ic,Jc),&
Pr(Ic,Jc),Pc(Ic,Jc),auP(Ic,Jc),auNB(Ic,Jc),aP(Ic,Jc),aW(Ic,Jc),aE(Ic,Jc),aS(Ic,Jc),aN(Ic,Jc),b(Ic,Jc),Xg(Ip,Jp),Yg(Ip,Jp),&
Xc(Ic,Jc),Yc(Ic,Jc),Xga(Ic,Jc),Xgk(Ic,Jc),Yga(Ic,Jc),Ygk(Ic,Jc),dk(Ic,Jc),da(Ic,Jc),Jg(Ic,Jc),a1(Ic,Jc),y1(Ic,Jc),b1(Ic,Jc),&
Un(Ic,Jc),Vn(Ic,Jc),Unk(Ip,Jc),Vna(Ic,Jp),duk(Ip,Jc),dva(Ic,Jp),Ux(Ic,Jc),Uy(Ic,Jc),Vx(Ic,Jc),Vy(Ic,Jc),Px(Ic,Jc),Py(Ic,Jc),&
dPx(Ic,Jc),dPy(Ic,Jc),muxx(Ic,Jc),muxy(Ic,Jc),muyx(Ic,Jc),mvxy(Ic,Jc),mvyx(Ic,Jc),mvyy(Ic,Jc),rouk(Ip,Jc),roua(Ic,Jp),&
d(Ic,Jc),St(Ic,Jc))
if(Turmod=='sa') then
allocate(Tn0(Ic,Jc),Tn(Ic,Jc),Tnx(Ic,Jc),Tny(Ic,Jc),Xi(Ic,Jc),fniu1(Ic,Jc))
else if(Turmod=='ke') then
allocate(Tk0(Ic,Jc),Te0(Ic,Jc),Tk(Ic,Jc),Te(Ic,Jc),roux(Ic,Jc),rouy(Ic,Jc))
else if(Turmod=='sst') then
allocate(Tk0(Ic,Jc),Tw0(Ic,Jc),Tk(Ic,Jc),Tw(Ic,Jc),Tkx(Ic,Jc),Tky(Ic,Jc),Twx(Ic,Jc),Twy(Ic,Jc),fai2(Ic,Jc),F2(Ic,Jc),&
Ret(Ic,Jc),alphastar(Ic,Jc),sigmatk(Ic,Jc),sigmatw(Ic,Jc))
end if
allocate(Xw(Ib1:Ib2),Yw(Ib1:Ib2),Yp(Ib1:Ib2),DR(Ib1:Ib2),Sw(Ib1:Ib2),ks(Ib1:Ib2),Q(Ib1:Ib2),Yplus(Ib1:Ib2),Ystar(Ib1:Ib2),&
ustar(Ib1:Ib2),Uplus(Ib1:Ib2),Tplus(Ib1:Ib2),hcv(Ib1:Ib2),Ax(Ib1:Ib2),Ay(Ib1:Ib2))
end Subroutine Allocarray
