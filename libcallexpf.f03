Program libaero2dcaller
  use aero2dcom
  implicit none
  integer i,Iw1,Iw2,Iw3
  character(1) mode
  character(64) filenamei
  mode='A'
  filenamei='case1/NACA0012.xyz'
  Proctrl='incom'
  Energy='Y'
  visheat='N'
  Turmod='sa'
  Walltreat='wf'
  solctrl='SIMPLE'
  Discret='2upwind'
  Init='N'
  Stag='N'
  maxs=2000
  delta=1d-5
  Rau=5d-1
  Rap=3d-1
  Rae=7d-1
  Rat=3d-1
  fd=1d-3
  Ifd=10
  Jp=75
  c=0.5334
  Vfar=75
  AoA=4
  Ta=263.15
  Tf=273.15
  Po=100000
  tvr=50
  ksi=0
  open(unit=1,file=filenamei,status='old')
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
  Iw1=max(Iwd,Iwu)
  Iw2=Iw1+Iwd-1
  Iw3=Iw2+Iwu-1
  Ip=Iw3+Iw1-1
  Ic=Ip-1
  Jc=Jp-1
  Ib1=Iw1
  Ib2=Iw3-1
  Call Allocarray
  Call Aero2D(mode,0,'')
  print *,Cl,Cd,Cf,Cm
  print *,Xpc,Ypc
  print *,maxval(hcv),minval(hcv)
  AoA=5
  Init='A'
  Call Aero2D(mode,0,'')
  print *,Cl,Cd,Cf,Cm
  print *,Xpc,Ypc
  print *,maxval(hcv),minval(hcv)
  Call Deallocarray
  deallocate(Xwd,Ywd,Xwu,Ywu)
end Program libaero2dcaller
