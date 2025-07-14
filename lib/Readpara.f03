Subroutine Readpara(libmod,scptname)
use Aero2DCOM
implicit none
character(1) ch
character(*) libmod,scptname
character(128) ioerrmsg
integer stat
if(libmod=='S'.or.libmod=='M') then
 open(unit=10,file=scptname,status='old',IOSTAT=stat,IOMSG=ioerrmsg)
 if(stat>0) stop ioerrmsg
  read(10,*) ch
  read(10,*) filename(1)
  read(10,*) ch
  read(10,*) Pntctrl
  read(10,*) ch
  read(10,*) Proctrl
  read(10,*) ch
  read(10,*) Energy
  if(Energy=='Y') then
   read(10,*) ch
   read(10,*) visheat
  end if
  read(10,*) ch
  read(10,*) Turmod
  if(Turmod=='sa'.or.Turmod=='sst') then
   read(10,*) ch
   read(10,*) Walltreat
  end if
  read(10,*) ch
  read(10,*) solctrl
  read(10,*) ch
  read(10,*) Discret
  if(Proctrl=='com') then
   read(10,*) ch
   read(10,*) denface
  end if
  read(10,*) ch
  read(10,*) Linsol
  if(Energy=='Y') then
   read(10,*) ch
   read(10,*) Tmptype
  end if
  read(10,*) ch
  read(10,*) Init
  if(Init=='Y') then
   read(10,*) ch
   read(10,*) filename(9)
  else if(Init=='N') then
   read(10,*) ch
   read(10,*) Stag
  end if
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) maxs
  if(stat>0) STOP 'max iter num: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) delta
  if(stat>0) STOP 'min residue: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Rau
  if(stat>0) STOP 'relax u: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Rap
  if(stat>0) STOP 'relax p: '//ioerrmsg
  if(Energy=='Y') then
   read(10,*) ch
   read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Rae
   if(stat>0) STOP 'relax e: '//ioerrmsg
  end if
  if(Turmod/='inv'.and.Turmod/='lam') then
   read(10,*) ch
   read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Rat
   if(stat>0) STOP 'relax tur: '//ioerrmsg
  end if
  if(Pntctrl=='Y') then
   read(10,*) ch
   read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Iw
   if(stat>0) STOP 'nump on profile: '//ioerrmsg
   read(10,*) ch
   read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) fb
   if(stat>0) STOP 'fd on LE: '//ioerrmsg
   read(10,*) ch
   read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) eb
   if(stat>0) STOP 'fd on TE: '//ioerrmsg
  end if
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) fd
  if(stat>0) STOP 'fd normal: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Ifd
  if(stat>0) STOP 'numlay normal: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Jp
  if(stat>0) STOP 'nump normal: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) c
  if(stat>0) STOP 'chord length: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Vfar
  if(stat>0) STOP 'velocity: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) AoA
  if(stat>0) STOP 'AoA: '//ioerrmsg
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Ta
  if(stat>0) STOP 'Free Temp: '//ioerrmsg
  if(Tmptype=='fixed') then
    read(10,*) ch
    read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Tf
    if(stat>0) STOP 'Wall Temp: '//ioerrmsg
  else if(Tmptype=='flux') then
    read(10,*) ch
    read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Qf
    if(stat>0) STOP 'Wall Flux: '//ioerrmsg
  end if
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Po
  if(stat>0) STOP 'ambient press: '//ioerrmsg
  if(Turmod=='ke'.or.Turmod=='sst') then
   read(10,*) ch
   read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) Itur
   if(stat>0) STOP 'tur intensity: '//ioerrmsg
   read(10,*) ch
   read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) tvr
   if(stat>0) STOP 'tur ratio: '//ioerrmsg
  else if(Turmod=='sa') then
   read(10,*) ch
   read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) tvr
   if(stat>0) STOP 'tur ratio: '//ioerrmsg
  end if
  read(10,*) ch
  read(10,*,IOSTAT=stat,IOMSG=ioerrmsg) ksi
  if(stat>0) STOP 'rough: '//ioerrmsg
  read(10,*) ch
  read(10,*) dir
  read(10,*,IOSTAT=stat) ch
  if(is_iostat_end(stat)) then
   gtype=''
  else
   read(10,*) gtype
  end if
 close(10)
else if(libmod=='I') then
 print *,'Input a name of airfoil coordinates file(1D XYZ/CPT):'
 read *,filename(1)
 print *,'Using control points based parametric spline(Y/N)?'
 read *,Pntctrl
 print *,'Select air property(com/incom):'
 read *,Proctrl
 print *,'Solve Energy equation(Y/N)?'
 read *,Energy
 if(Energy=='Y') then
  print *,'Include Viscous Heating(Y/N)?'
  read *,visheat
 end if
 print *,'Select turbulence model(inv/lam/sa/ke/sst):'
 read *,Turmod
 if(Turmod=='sa'.or.Turmod=='sst') then
  print *,'Select wall treatment method(wf/lr):'
  read *,Walltreat
 end if
 print *,'Select coupling algorithm(SIMPLE/SIMPLEC):'
 read *,solctrl
 print *,'Select discretization scheme(1upwind/2upwind/Quick/tvd):'
 read *,Discret
 if(Proctrl=='com') then
  print *,'Select density interpolation scheme(center/1upwind/2upwind/Quick/tvd):'
  read *,denface
 end if
 print *,'Select type of linear solver for convective-diffusion equations(sor/pbicg):'
 read *,Linsol
 if(Energy=='Y') then
   print *,'Set boundary condition type of tempreture(fixed/flux):'
   read *, Tmptype
 end if
 print *,'Initialize from a file(Y/N)?'
 read *,Init
 if(Init=='Y') then
  print *,'Input a file name for initialization:'
  read *,filename(9)
 else if(Init=='N') then
  print *,'Initialize using stagnation values(Y/N)?'
  read *,Stag
 end if
 print *,'Input maximum iteration steps:'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) maxs
 if(stat>0) STOP ioerrmsg
 print *,'Input minimum residual of iteration:'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) delta
 if(stat>0) STOP ioerrmsg
 print *,'Input relaxation factor of velocity:'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Rau
 if(stat>0) STOP ioerrmsg
 print *,'Input relaxation factor of pressure:'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Rap
 if(stat>0) STOP ioerrmsg
 if(Energy=='Y') then
  print *,'Input relaxation factor of temperature:'
  read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Rae
  if(stat>0) STOP ioerrmsg
 end if
 if(Turmod/='inv'.and.Turmod/='lam') then
  print *,'Input relaxation factor of turbulence:'
  read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Rat
  if(stat>0) STOP ioerrmsg
 end if
 if(Pntctrl=='Y') then
   print *,'Input number of grid points on whole airfoil(interpolation based on parametric spline):'
   read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Iw
   if(stat>0) STOP ioerrmsg
   print *,'Input dimensionless mesh spacing near leading edge:'
   read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) fb
   if(stat>0) STOP ioerrmsg
   print *,'Input dimensionless mesh spacing near trailing edge:'
   read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) eb
   if(stat>0) STOP ioerrmsg
 end if
 print *,'Input dimensionless near wall mesh spacing:'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) fd
 if(stat>0) STOP ioerrmsg
 print *,'Input layers of uniform near wall mesh:'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Ifd
 if(stat>0) STOP ioerrmsg
 print *,'Input total number of grid points normal to the airfoil surface:'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Jp
 if(stat>0) STOP ioerrmsg
 print *,'Input chord length(m):'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) c
 if(stat>0) STOP ioerrmsg
 print *,'Input velocity of free stream(m/s):'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Vfar
 if(stat>0) STOP ioerrmsg
 print *,'Input angle of attack(deg):'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) AoA
 if(stat>0) STOP ioerrmsg
 print *,'Input temperature of free stream(K):'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Ta
 if(stat>0) STOP ioerrmsg
 if(Tmptype=='fixed') then
   print *,'Input temperature of wall(K):'
   read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Tf
   if(stat>0) STOP ioerrmsg
 else if(Tmptype=='flux') then
   print *,'Input heat flux of wall(J/m^2s):'
   read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Qf
   if(stat>0) STOP ioerrmsg
 end if
 print *,'Input ambient pressure(Pa):'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Po
 if(stat>0) STOP ioerrmsg
 if(Turmod=='ke'.or.Turmod=='sst') then
  print *,'Input turbulence intensity of free stream:'
  read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) Itur
  if(stat>0) STOP ioerrmsg
  print *,'Input turbulent viscosity ratio of free stream:'
  read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) tvr
  if(stat>0) STOP ioerrmsg
 else if(Turmod=='sa') then
  print *,'Input turbulent viscosity ratio of free stream:'
  read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) tvr
  if(stat>0) STOP ioerrmsg
 end if
 print *,'Input wall roughness(m):'
 read (*,*,IOSTAT=stat,IOMSG=ioerrmsg) ksi
 if(stat>0) STOP ioerrmsg
 print *,'Input a directory of output files:'
 read *,dir
 print *,'Input topological type of generating mesh(C/O):'
 read *,gtype
end if
print *,'Read control parameters completed!'
end Subroutine Readpara
