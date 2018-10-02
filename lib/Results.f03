Subroutine Saveresults
use Aero2DCOM
implicit none
integer i,j,ioerr
real(8) Ug(Ip,Jp),Vg(Ip,Jp),Pg(Ip,Jp),rhog(Ip,Jp),Tg(Ip,Jp),Tng(Ip,Jp),Tkg(Ip,Jp),Teg(Ip,Jp),Twg(Ip,Jp)
Call Cell2node(Ug,U,'U')
Call Cell2node(Vg,V,'V')
Call Cell2node(Pg,P,'P')
Call Cell2node(rhog,rho,'rho')
Call Cell2node(Tg,T,'T')
if(Turmod=='sa') then
 Call Cell2node(Tng,Tn,'Tn')
else if(Turmod=='ke') then
 Call Cell2node(Tkg,Tk,'Tk')
 Call Cell2node(Teg,Te,'Te')
else if(Turmod=='sst') then
 Call Cell2node(Tkg,Tk,'Tk')
 Call Cell2node(Twg,Tw,'Tw')
end if
filename(3)=trim(dir)//'/Autosave.dat'
filename(4)=trim(dir)//'/Grid.xyz'
filename(5)=trim(dir)//'/Solutions.dat'
filename(6)=trim(dir)//'/Solnames.nam'
filename(7)=trim(dir)//'/Wallsol.dat'
filename(8)=trim(dir)//'/Aeroreport.dat'
open(unit=3,file=filename(3),form='unformatted',status='replace')
  write(3,iostat=ioerr) Xc
  write(3,iostat=ioerr) Yc
  write(3,iostat=ioerr) U
  write(3,iostat=ioerr) V
  write(3,iostat=ioerr) P
  write(3,iostat=ioerr) T
  write(3,iostat=ioerr) rho
  if(Turmod=='sa') then
   write(3,iostat=ioerr) Tn
  else if(Turmod=='ke') then
   write(3,iostat=ioerr) Tk
   write(3,iostat=ioerr) Te
  else if(Turmod=='sst') then
   write(3,iostat=ioerr) Tk
   write(3,iostat=ioerr) Tw
  end if
close(3)
open(unit=3,file=filename(4),status='replace')
 write(3,*) 1
 write(3,*) Ip,Jp
 write(3,*) ((Xg(I,J),I=1,Ip),J=1,Jp)
 write(3,*) ((Yg(I,J),I=1,Ip),J=1,Jp)
close(3)
open(unit=4,file=filename(5),status='replace')
 write(4,*) 1
 if(Turmod=='sa') then
  write(4,*) Ip,Jp,6
 else if(Turmod=='ke'.or.Turmod=='sst') then
  write(4,*) Ip,Jp,7
 else if(Turmod=='lam'.or.Turmod=='inv') then
  write(4,*) Ip,Jp,5
 end if
 write(4,*) ((Ug(i,j),i=1,Ip),j=1,Jp)
 write(4,*) ((Vg(i,j),i=1,Ip),j=1,Jp)
 write(4,*) ((Pg(i,j),i=1,Ip),j=1,Jp)
 write(4,*) ((Tg(i,j),i=1,Ip),j=1,Jp)
 write(4,*) ((rhog(i,j),i=1,Ip),j=1,Jp)
 if(Turmod=='sa') then
  write(4,*) ((Tng(i,j),i=1,Ip),j=1,Jp)
 else if(Turmod=='ke') then
  write(4,*) ((Tkg(i,j),i=1,Ip),j=1,Jp)
  write(4,*) ((Teg(i,j),i=1,Ip),j=1,Jp)
 else if(Turmod=='sst') then
  write(4,*) ((Tkg(i,j),i=1,Ip),j=1,Jp)
  write(4,*) ((Twg(i,j),i=1,Ip),j=1,Jp)
 end if
close(4)
open(unit=4,file=filename(6),status='replace')
 write(4,*) 'X-velocity'
 write(4,*) 'Y-velocity'
 write(4,*) 'static pressure'
 write(4,*) 'static temperature'
 write(4,*) 'density'
 if(Turmod=='sa') then
  write(4,*) 'modified turbulent kinematic viscosity'
 else if(Turmod=='ke') then
  write(4,*) 'turbulence kinetic energy'
  write(4,*) 'dissipation rate'
 else if(Turmod=='sst') then
  write(4,*) 'turbulence kinetic energy'
  write(4,*) 'specific dissipation rate'
 end if
close(4)
filename(6)=trim(dir)//'/GridSol.vtk'
open(unit=4,file=filename(6),status='replace')
 write(4,'(A26)') '# vtk DataFile Version 4.0'
 write(4,'(A36)') '2D Structured Grid and Solution Data'
 write(4,'(A5)') 'ASCII'
 write(4,'(A23)') 'DATASET STRUCTURED_GRID'
 write(4,'(A10,1X,I4,1X,I4,1X,I4)') 'DIMENSIONS',Ip,Jp,1
 write(4,'(A6,1X,I6,1X,A6)') 'POINTS',Ip*Jp,'double'
 DO j=1,Jp
  DO i=1,Ip
   write(4,*) Xg(i,j),Yg(i,j),0.0
  end DO
 end DO
 write(4,'(A10,1X,I6)') 'POINT_DATA',Ip*Jp
 write(4,'(A25)') 'SCALARS Pressure double 1'
 write(4,'(A20)') 'LOOKUP_TABLE default'
 DO j=1,Jp
  DO i=1,Ip
   write(4,*) Pg(i,j)
  end DO
 end DO
 write(4,'(A28)') 'SCALARS Temperature double 1'
 write(4,'(A20)') 'LOOKUP_TABLE default'
 DO j=1,Jp
  DO i=1,Ip
   write(4,*) Tg(i,j)
  end DO
 end DO
 write(4,'(A24)') 'SCALARS Density double 1'
 write(4,'(A20)') 'LOOKUP_TABLE default'
 DO j=1,Jp
  DO i=1,Ip
   write(4,*) rhog(i,j)
  end DO
 end DO
 if(Turmod=='sa') then
  write(4,'(A23)') 'SCALARS Turvis double 1'
  write(4,'(A20)') 'LOOKUP_TABLE default'
  DO j=1,Jp
   DO i=1,Ip
    write(4,*) Tng(i,j)
   end DO
  end DO
 else if(Turmod=='ke') then
  write(4,'(A22)') 'SCALARS Turke double 1'
  write(4,'(A20)') 'LOOKUP_TABLE default'
  DO j=1,Jp
   DO i=1,Ip
    write(4,*) Tkg(i,j)
   end DO
  end DO
  write(4,'(A22)') 'SCALARS Turdr double 1'
  write(4,'(A20)') 'LOOKUP_TABLE default'
  DO j=1,Jp
   DO i=1,Ip
    write(4,*) Teg(i,j)
   end DO
  end DO
 else if(Turmod=='kw') then
  write(4,'(A22)') 'SCALARS Turke double 1'
  write(4,'(A20)') 'LOOKUP_TABLE default'
  DO j=1,Jp
   DO i=1,Ip
    write(4,*) Tkg(i,j)
   end DO
  end DO
  write(4,'(A23)') 'SCALARS Tursdr double 1'
  write(4,'(A20)') 'LOOKUP_TABLE default'
  DO j=1,Jp
   DO i=1,Ip
    write(4,*) Twg(i,j)
   end DO
  end DO
 end if
 write(4,'(A23)') 'VECTORS Velocity double'
 DO j=1,Jp
  DO i=1,Ip
   write(4,*) Ug(i,j),Vg(i,j),0.0
  end DO
 end DO
close(4)
open(unit=4,file=filename(7),status='replace')
  write(4,'(11(A11,2X))',advance='no') '#        X=','Y=','S=','P=','U=','V=','T=','mut=','hcv=','Ax=','Ay='
  if(Turmod=='sa'.or.Turmod=='sst') then
   write(4,'(A11)') 'Yplus='
  else if(Turmod=='ke') then
   write(4,'(A11)') 'Ystar='
  else
   write(4,*) ''
  end if
  DO i=Ib1,Ib2
   write(4,'(11(ES11.4,2X))',advance='no') Xw(i),Yw(i),Sw(i),P(i,1),U(i,1),V(i,1),T(i,1),mut(i,1),hcv(i),Ax(i),Ay(i)
   if(Turmod=='sa'.or.Turmod=='sst') then
    write(4,'(ES11.4)') Yplus(i)
   else if(Turmod=='ke') then
    write(4,'(ES11.4)') Ystar(i)
   else
    write(4,*) ''
   end if
  end DO
close(4)
open(unit=4,file=filename(8),status='replace')
 write(4,*) 'File name of airfoil coordinates: ',filename(1)
 write(4,*) 'Property: ',Proctrl
 write(4,*) 'Solving energy equation: ',Energy
 if(Energy=='Y') then
  write(4,*) 'Including viscous heating: ',visheat
 end if
 write(4,*) 'Turbulence model: ',Turmod
 if(Turmod=='sa'.or.Turmod=='sst') then
  write(4,*) 'Wall treatment method: ',Walltreat
 end if
 write(4,*) 'Coupling algorithm: ',solctrl
 write(4,*) 'Spatial discretization scheme: ',Discret
 if(Proctrl=='com') then
  write(4,*) 'Density interpolation scheme: ',denface
 end if
 write(4,*) 'Initializing from a file: ',Init
 if(Init=='Y') then
  write(4,*) 'File name for initialization: ',filename(9)
 else if(Init=='N') then
  write(4,*) 'Initializing using stagnation values: ',Stag
 end if
 write(4,*) 'Momentum relaxation factor:',Rau
 write(4,*) 'Pressure relaxation factor:',Rap
 if(Energy=='Y') then
  write(4,*) 'Energy relaxation factor:',Rae
 end if
 if(Turmod/='inv'.and.Turmod/='lam') then
  write(4,*) 'Turbulence relaxation factor:',Rat
 end if
 write(4,*) '----------------------------------------------'
 write(4,*) 'The aerodynamic parameters of this airfoil are:'
 write(4,*) 'Reynolds number:',Re
 write(4,*) 'Mach number:',Mach
 write(4,*) 'Angle of attack:',AoA
 if(Turmod=='sa'.or.Turmod=='sst') then
  write(4,*) 'Average Y+:',sum(Yplus)/(Ib2-Ib1+1)
 else if(Turmod=='ke') then
  write(4,*) 'Average Y*:',sum(Ystar)/(Ib2-Ib1+1)
 end if
 write(4,*) 'Lift coefficient:',Cl
 write(4,*) 'Drag coefficient:',Cd
 write(4,*) 'Friction coefficient:',Cf
 write(4,*) 'Pitching moment coefficient (1/4 chord):',Cm
 write(4,*) 'Pressure center (unit chord):',Xpc,Ypc
close(4)
end Subroutine Saveresults
