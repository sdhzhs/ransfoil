Subroutine Genmesh(libmod)
use Aero2DCOM
implicit none
integer i,j,l,Iw1,Iw2,Iw3,maxl
real(8) err,fb,ltrail,ratio,ratio0
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
Jp=75
Ic=Ip-1
Jc=Jp-1
Ib1=Iw1
Ib2=Iw3-1
Call Allocarray
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
ratio=2.0
maxl=100
err=1d-10
DO l=1,maxl
ratio0=ratio
ratio=ratio-(1-ratio**(Jp-6)-(10/fd-5)*(1-ratio))/(-(Jp-6)*ratio**(Jp-7)+10/fd-5)
if(abs(ratio-ratio0)<err) exit
end DO
DO j=2,Jp
if(j<=6) then
 fac(j)=1.0
else
 fac(j)=ratio
end if
end DO
Call hypmeshgen(Xg(:,1),Yg(:,1),Xg,Yg,fac,fd,Ip,Jp)
print *,'Generate 2D C-type structured mesh completed!'
deallocate(Xt,fac)
end Subroutine Genmesh

Subroutine Allocarray
use Aero2DCOM
implicit none
allocate(U0(Ic,Jc),V0(Ic,Jc),T0(Ic,Jc),Tn0(Ic,Jc),Tk0(Ic,Jc),Te0(Ic,Jc),Tw0(Ic,Jc),rou(Ic,Jc),miu(Ic,Jc),P(Ic,Jc),dP(Ic,Jc),&
U(Ic,Jc),V(Ic,Jc),T(Ic,Jc),Tn(Ic,Jc),Tk(Ic,Jc),Te(Ic,Jc),Tw(Ic,Jc),miut(Ic,Jc),Pr(Ic,Jc),Pc(Ic,Jc),auP(Ic,Jc),auW(Ic,Jc),&
auE(Ic,Jc),auS(Ic,Jc),auN(Ic,Jc),bu(Ic,Jc),bv(Ic,Jc),apP(Ic,Jc),apW(Ic,Jc),apE(Ic,Jc),apS(Ic,Jc),apN(Ic,Jc),bp(Ic,Jc),atP(Ic,Jc),&
atW(Ic,Jc),atE(Ic,Jc),atS(Ic,Jc),atN(Ic,Jc),bt(Ic,Jc),anP(Ic,Jc),anW(Ic,Jc),anE(Ic,Jc),anS(Ic,Jc),anN(Ic,Jc),bn(Ic,Jc),akP(Ic,Jc),&
akW(Ic,Jc),akE(Ic,Jc),akS(Ic,Jc),akN(Ic,Jc),bk(Ic,Jc),aeP(Ic,Jc),aeW(Ic,Jc),aeE(Ic,Jc),aeS(Ic,Jc),aeN(Ic,Jc),be(Ic,Jc),awP(Ic,Jc),&
awW(Ic,Jc),awE(Ic,Jc),awS(Ic,Jc),awN(Ic,Jc),bw(Ic,Jc),Xg(Ip,Jp),Yg(Ip,Jp),Xc(Ic,Jc),Yc(Ic,Jc),Xga(Ic,Jc),Xgk(Ic,Jc),Yga(Ic,Jc),&
Ygk(Ic,Jc),Jg(Ic,Jc),a1(Ic,Jc),y1(Ic,Jc),b1(Ic,Jc),Un(Ic,Jc),Vn(Ic,Jc),Unw(Ic,Jc),Une(Ic,Jc),Vns(Ic,Jc),Vnn(Ic,Jc),wdu(Ic,Jc),&
edu(Ic,Jc),sdv(Ic,Jc),ndv(Ic,Jc),Ux(Ic,Jc),Uy(Ic,Jc),Vx(Ic,Jc),Vy(Ic,Jc),Tnx(Ic,Jc),Tny(Ic,Jc),Tkx(Ic,Jc),Tky(Ic,Jc),Twx(Ic,Jc),&
Twy(Ic,Jc),Px(Ic,Jc),Py(Ic,Jc),dPx(Ic,Jc),dPy(Ic,Jc),roux(Ic,Jc),rouy(Ic,Jc),Rkx(Ic,Jc),Rky(Ic,Jc),muxx(Ic,Jc),muxy(Ic,Jc),&
muyx(Ic,Jc),mvxy(Ic,Jc),mvyx(Ic,Jc),mvyy(Ic,Jc),ww(Ic,Jc),we(Ic,Jc),ws(Ic,Jc),wn(Ic,Jc),rouw(Ic,Jc),roue(Ic,Jc),rous(Ic,Jc),&
roun(Ic,Jc),Xi(Ic,Jc),fniu1(Ic,Jc),d(Ic,Jc),fai2(Ic,Jc),F2(Ic,Jc),St(Ic,Jc),Ret(Ic,Jc),alphastar(Ic,Jc),sigmatk(Ic,Jc),&
sigmatw(Ic,Jc))
allocate(Xw(Ib1:Ib2),Yw(Ib1:Ib2),Dyp(Ib1:Ib2),DR(Ib1:Ib2),Sw(Ib1:Ib2),ks(Ib1:Ib2),Q(Ib1:Ib2),Yplus(Ib1:Ib2),Ystar(Ib1:Ib2),&
ustar(Ib1:Ib2),Uplus(Ib1:Ib2),Tplus(Ib1:Ib2),hcv(Ib1:Ib2),Ax(Ib1:Ib2),Ay(Ib1:Ib2))
end Subroutine Allocarray
