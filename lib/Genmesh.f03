Subroutine Genmesh(libmod)
use Aero2DCOM
implicit none
integer i,j,l,Iw1,Iw2,Iw3,maxl,It
logical(1) opentrail
real(8) err,ft,ltrail,lfar,ratio,ratio0
real(8)::tol=1e-8
real(8),allocatable,dimension(:)::Xt,fac
character(*) libmod

opentrail=.false.
if(libmod=='S'.or.libmod=='I') then
 Call Readfoil
else if(libmod=='M') then
 Call Readmesh(libmod)
 return
else if(libmod=='C') then
 if(Pntctrl=='N') then
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
 else
  allocate(Xwp0(Iw0),Ywp0(Iw0))
  Call C_F_POINTER(cXwp0,fXwp0,(/Iw0/))
  Call C_F_POINTER(cYwp0,fYwp0,(/Iw0/))
  Xwp0=fXwp0
  Ywp0=fYwp0
 end if
end if
if(Pntctrl=='Y') then
 allocate(Xwp(Iw),Ywp(Iw))
 Call Connector2(Xwp,Ywp,Xwp0,Ywp0,fb,eb,Iw,Iw0)
 Iwd=(Iw+1)/2
 Iwu=Iwd
 if(abs(Ywp(Iw)-Ywp(1))>tol) then
  opentrail=.true.
  It=6
  Iwu=Iwu+It
 end if
 allocate(Xwd(Iwd),Ywd(Iwd))
 allocate(Xwu(Iwu),Ywu(Iwu))
 Xwd=Xwp((Iw+1)/2:Iw)
 Ywd=Ywp((Iw+1)/2:Iw)
 Xwu(1:(Iw+1)/2)=Xwp((Iw+1)/2:1:-1)
 Ywu(1:(Iw+1)/2)=Ywp((Iw+1)/2:1:-1)
 if(abs(Ywp(Iw)-Ywp(1))>tol) then
  DO i=1,It
   Xwu((Iw+1)/2+i)=Xwp(Iw)
   Ywu((Iw+1)/2+i)=Ywp(1)+i*(Ywp(Iw)-Ywp(1))/It
  end DO
 end if
else
 if(abs(Xwd(Iwd)-Xwd(Iwd-1))<tol.or.abs(Xwu(Iwu)-Xwu(Iwu-1))<tol) then
  opentrail=.true.
 end if
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
ft=sqrt((Xwd(Iwd)-Xwd(Iwd-1))**2+(Ywd(Iwd)-Ywd(Iwd-1))**2)
ltrail=11
Call tanhgridline(ltrail,ft,Xt,Iw1)
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
Call hypmeshgen(Xg(:,1),Yg(:,1),Xg,Yg,fac,fd,Ip,Jp,opentrail)
print *,'Generate 2D C-type structured mesh completed!'
deallocate(Xt,fac)

end Subroutine Genmesh

Subroutine Readfoil
use Aero2DCOM
implicit none
integer i,stat
character(128) ioerrmsg

open(unit=1,file=filename(1),status='old',IOSTAT=stat,IOMSG=ioerrmsg)
if(stat>0) stop ioerrmsg
 if(Pntctrl=='N') then
  read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) Iwd
  allocate(Xwd(Iwd),Ywd(Iwd))
  DO i=1,Iwd
   read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) Xwd(i),Ywd(i)
  end DO
  read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) Iwu
  allocate(Xwu(Iwu),Ywu(Iwu))
  DO i=1,Iwu
   read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) Xwu(i),Ywu(i)
  end DO
 else
  read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) Iw0
  allocate(Xwp0(Iw0),Ywp0(Iw0))
  DO i=1,Iw0
   read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) Xwp0(i),Ywp0(i)
  end DO
 end if
if(stat>0) stop filename(1)//ioerrmsg
close(1)
print *,'Read airfoil coordinates completed!'
if(Pntctrl=='N') then
 if(Ywd(2)>=Ywu(2)) then
  stop 'Error: The coordinates of lower airfoil should be given first in the xyz file!'
 end if
else
 if(Ywp0(2)<Ywp0(Iw0-1)) then
  stop 'Error: The coordinates of upper airfoil should be given first in the cpt file!'
 end if
end if

end Subroutine Readfoil

Subroutine Readmesh(libmod)
use Aero2DCOM
implicit none
integer blocks,i,j,stat
character(128) ioerrmsg
character(*) libmod

open(unit=1,file=filename(1),status='old',IOSTAT=stat,IOMSG=ioerrmsg)
if(stat>0) stop ioerrmsg
 read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) blocks
 read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) Ip,Jp,Ib1,Ib2
 if(stat>0) stop filename(1)//ioerrmsg
 Ic=Ip-1
 Jc=Jp-1
 Call Allocarray(libmod)
 read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) ((Xg(i,j),i=1,Ip),j=1,Jp)
 read(1,*,IOSTAT=stat,IOMSG=ioerrmsg) ((Yg(i,j),i=1,Ip),j=1,Jp)
 if(stat>0) stop filename(1)//ioerrmsg
close(1)
print *,'Read airfoil 2D C-type mesh completed!'

end Subroutine Readmesh
