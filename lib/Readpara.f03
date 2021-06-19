Subroutine Readpara(libmod,scptname)
use Aero2DCOM
implicit none
character(1) ch
character(*) libmod,scptname
if(libmod=='S') then
 open(unit=10,file=scptname,status='old')
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
  read(10,*) maxs
  read(10,*) ch
  read(10,*) delta
  read(10,*) ch
  read(10,*) Rau
  read(10,*) ch
  read(10,*) Rap
  if(Energy=='Y') then
   read(10,*) ch
   read(10,*) Rae
  end if
  if(Turmod/='inv'.and.Turmod/='lam') then
   read(10,*) ch
   read(10,*) Rat
  end if
  if(Pntctrl=='Y') then
   read(10,*) ch
   read(10,*) Iw
   read(10,*) ch
   read(10,*) fb
   read(10,*) ch
   read(10,*) eb
  end if
  read(10,*) ch
  read(10,*) fd
  read(10,*) ch
  read(10,*) Ifd
  read(10,*) ch
  read(10,*) Jp
  read(10,*) ch
  read(10,*) c
  read(10,*) ch
  read(10,*) Vfar
  read(10,*) ch
  read(10,*) AoA
  read(10,*) ch
  read(10,*) Ta
  read(10,*) ch
  read(10,*) Tf
  read(10,*) ch
  read(10,*) Po
  if(Turmod=='ke'.or.Turmod=='sst') then
   read(10,*) ch
   read(10,*) Itur
   read(10,*) ch
   read(10,*) tvr
  else if(Turmod=='sa') then
   read(10,*) ch
   read(10,*) tvr
  end if
  read(10,*) ch
  read(10,*) ksi
  read(10,*) ch
  read(10,*) dir
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
 read *,maxs
 print *,'Input minimum residual of iteration:'
 read *,delta
 print *,'Input relaxation factor of velocity:'
 read *,Rau
 print *,'Input relaxation factor of pressure:'
 read *,Rap
 if(Energy=='Y') then
  print *,'Input relaxation factor of temperature:'
  read *,Rae
 end if
 if(Turmod/='inv'.and.Turmod/='lam') then
  print *,'Input relaxation factor of turbulence:'
  read *,Rat
 end if
 if(Pntctrl=='Y') then
   print *,'Input number of grid points on whole airfoil(interpolation based on parametric spline):'
   read *,Iw
   print *,'Input dimensionless mesh spacing near leading edge:'
   read *,fb
   print *,'Input dimensionless mesh spacing near trailing edge:'
   read *,eb
 end if
 print *,'Input dimensionless near wall mesh spacing:'
 read *,fd
 print *,'Input layers of uniform near wall mesh:'
 read *,Ifd
 print *,'Input total number of grid points normal to the airfoil surface:'
 read *,Jp
 print *,'Input chord length(m):'
 read *,c
 print *,'Input velocity of free stream(m/s):'
 read *,Vfar
 print *,'Input angle of attack(deg):'
 read *,AoA
 print *,'Input temperature of free stream(K):'
 read *,Ta
 print *,'Input temperature of wall(K):'
 read *,Tf
 print *,'Input ambient pressure(Pa):'
 read *,Po
 if(Turmod=='ke'.or.Turmod=='sst') then
  print *,'Input turbulence intensity of free stream:'
  read *,Itur
  print *,'Input turbulent viscosity ratio of free stream:'
  read *,tvr
 else if(Turmod=='sa') then
  print *,'Input turbulent viscosity ratio of free stream:'
  read *,tvr
 end if
 print *,'Input wall roughness(m):'
 read *,ksi
 print *,'Input a directory of output files:'
 read *,dir
end if
print *,'Read control parameters completed!'
end Subroutine Readpara