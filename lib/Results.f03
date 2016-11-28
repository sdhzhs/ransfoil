Subroutine Results
use Aero2DCOM
implicit none
integer i,j,ioerr
real(8) Ug(Ip,Jp),Vg(Ip,Jp),Pg(Ip,Jp),roug(Ip,Jp),Tg(Ip,Jp),Tng(Ip,Jp),Tkg(Ip,Jp),Teg(Ip,Jp),Twg(Ip,Jp)
DO j=1,Jc
  DO i=2,Ic
  if(j==1.and.i>=Ib1.and.i<=Ib2+1) then
  if(Turmod/='inv') then
  Ug(i,j)=0
  Vg(i,j)=0
  else
  Ug(i,j)=0.5*(U(i,j)+U(i-1,j))
  Vg(i,j)=0.5*(V(i,j)+V(i-1,j))
  end if
  Pg(i,j)=0.5*(P(i,j)+P(i-1,j))
  Tg(i,j)=Tf
  if(Proctrl=='com') then
  roug(i,j)=(Pg(i,j)+Po)*Ma/(R*Tg(i,j))
  else if(Proctrl=='incom') then
  roug(i,j)=roui
  end if
  if(Turmod=='sa') then
  Tng(i,j)=0
  else if(Turmod=='ke') then
  Tkg(i,j)=0
  Teg(i,j)=0.5*(Te(i,j)+Te(i-1,j))
  else if(Turmod=='sst') then
  Tkg(i,j)=0
  Twg(i,j)=0.5*(Tw(i,j)+Tw(i-1,j))
  end if
  else if(j==1) then
  Ug(i,j)=(U(i,j)+U(i-1,j)+U(Ic+1-i,j)+U(Ic+2-i,j))/4
  Vg(i,j)=(V(i,j)+V(i-1,j)+V(Ic+1-i,j)+V(Ic+2-i,j))/4
  Pg(i,j)=(P(i,j)+P(i-1,j)+P(Ic+1-i,j)+P(Ic+2-i,j))/4
  Tg(i,j)=(T(i,j)+T(i-1,j)+T(Ic+1-i,j)+T(Ic+2-i,j))/4
  roug(i,j)=(rou(i,j)+rou(i-1,j)+rou(Ic+1-i,j)+rou(Ic+2-i,j))/4
  if(Turmod=='sa') then
  Tng(i,j)=(Tn(i,j)+Tn(i-1,j)+Tn(Ic+1-i,j)+Tn(Ic+2-i,j))/4
  else if(Turmod=='ke') then
  Tkg(i,j)=(Tk(i,j)+Tk(i-1,j)+Tk(Ic+1-i,j)+Tk(Ic+2-i,j))/4
  Teg(i,j)=(Te(i,j)+Te(i-1,j)+Te(Ic+1-i,j)+Te(Ic+2-i,j))/4
  else if(Turmod=='sst') then
  Tkg(i,j)=(Tk(i,j)+Tk(i-1,j)+Tk(Ic+1-i,j)+Tk(Ic+2-i,j))/4
  Twg(i,j)=(Tw(i,j)+Tw(i-1,j)+Tw(Ic+1-i,j)+Tw(Ic+2-i,j))/4
  end if
  else
  Ug(i,j)=(U(i,j)+U(i-1,j)+U(i,j-1)+U(i-1,j-1))/4
  Vg(i,j)=(V(i,j)+V(i-1,j)+V(i,j-1)+V(i-1,j-1))/4
  Pg(i,j)=(P(i,j)+P(i-1,j)+P(i,j-1)+P(i-1,j-1))/4
  Tg(i,j)=(T(i,j)+T(i-1,j)+T(i,j-1)+T(i-1,j-1))/4
  roug(i,j)=(rou(i,j)+rou(i-1,j)+rou(i,j-1)+rou(i-1,j-1))/4
  if(Turmod=='sa') then
  Tng(i,j)=(Tn(i,j)+Tn(i-1,j)+Tn(i,j-1)+Tn(i-1,j-1))/4
  else if(Turmod=='ke') then
  Tkg(i,j)=(Tk(i,j)+Tk(i-1,j)+Tk(i,j-1)+Tk(i-1,j-1))/4
  Teg(i,j)=(Te(i,j)+Te(i-1,j)+Te(i,j-1)+Te(i-1,j-1))/4
  else if(Turmod=='sst') then
  Tkg(i,j)=(Tk(i,j)+Tk(i-1,j)+Tk(i,j-1)+Tk(i-1,j-1))/4
  Twg(i,j)=(Tw(i,j)+Tw(i-1,j)+Tw(i,j-1)+Tw(i-1,j-1))/4
  end if
  end if
  end DO
end DO
  Ug(:,Jp)=Ui
  Vg(:,Jp)=Vi
  Pg(:,Jp)=0
  Tg(:,Jp)=Ta
  roug(:,Jp)=roui
  Ug(1,:)=Ui
  Ug(Ip,:)=Ui
  Vg(1,:)=Vi
  Vg(Ip,:)=Vi
  Pg(1,:)=0
  Pg(Ip,:)=0
  Tg(1,:)=Ta
  Tg(Ip,:)=Ta
  roug(1,:)=roui
  roug(Ip,:)=roui
if(Turmod=='sa') then
  Tng(:,Jp)=Tni
  Tng(1,:)=Tni
  Tng(Ip,:)=Tni
else if(Turmod=='ke') then
  Tkg(:,Jp)=Tki
  Teg(:,Jp)=Tei
  Tkg(1,:)=Tki
  Tkg(Ip,:)=Tki
  Teg(1,:)=Tei
  Teg(Ip,:)=Tei
else if(Turmod=='sst') then
  Tkg(:,Jp)=Tki
  Twg(:,Jp)=Twi
  Tkg(1,:)=Tki
  Tkg(Ip,:)=Tki
  Twg(1,:)=Twi
  Twg(Ip,:)=Twi
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
    write(3,iostat=ioerr) rou
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
write(4,*) ((roug(i,j),i=1,Ip),j=1,Jp)
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
  write(4,*) roug(i,j)
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
   write(4,'(11(A11,2X))',advance='no') '#        X=','Y=','S=','P=','U=','V=','T=','miut=','hcv=','Ax=','Ay='
   if(Turmod=='sa'.or.Turmod=='sst') then
   write(4,'(A11)') 'Yplus='
   else if(Turmod=='ke') then
   write(4,'(A11)') 'Ystar='
   else
   write(4,*) ''
   end if
   DO i=Ib1,Ib2
    write(4,'(11(ES11.4,2X))',advance='no') Xw(i),Yw(i),Sw(i),P(i,1),U(i,1),V(i,1),T(i,1),miut(i,1),hcv(i),Ax(i),Ay(i)
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
end Subroutine Results
