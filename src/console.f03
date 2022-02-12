Program RANSFOIL
  implicit none
  integer argc
  character(1) mode,termi
  character(64) argv(3)
  character(64) scptname
  argc=command_argument_count()
  Call get_command_argument(0,argv(1))
  if(argc==0) then
   print *,"usage: ",trim(argv(1))," [options] [configuration file]"
   print *,"Console program to mainly calculate aerodynamic parameters of an airfoil by numerically solving the RANS equations. &
   &Optionally the heat transfer coefficient on airfoil surface can also be obtained. This program reads airfoil coordinates of grid points (or control points) &
   &from a 1D XYZ/CPT file, then generates grid and solution files in 2D PLOT3D format to record grid and airflow data and a report file to show aerodynamic &
   &parameters."
   print *,"More information can be found in the README file."
   print *,"Options:"
   print *,"--stdin        interactive mode, read parameters from standard input"
   print *,"--script       batch mode, read parameters from a configuration file"
   print *,"--version      show the current version"
   print *,"-h, --help     show this message"
  else if(argc==1) then
   Call get_command_argument(1,argv(2))
   if(trim(argv(2))=="--stdin") then
    mode='I'
    scptname="Null"
    termi='N'
    DO while(termi=='N')
     Call Aero2D(mode,1,scptname)
     print *,'Terminate the interactive mode?(Y/N)'
     read *,termi
    end DO
   else if(trim(argv(2))=="--script") then
    print *,"should input a name of the configuration file."
    print *,"usage: ",trim(argv(1))," [options] [configuration file]"
   else if(trim(argv(2))=="--version") then
    print *,"ransfoil version 2.2.10"
    print *,"Copyright (c) 2012-2022, Hou Shuo"
   else if(trim(argv(2))=="-h".or.trim(argv(2))=="--help") then
    print *,"usage: ",trim(argv(1))," [options] [configuration file]"
    print *,"Console program to mainly calculate aerodynamic parameters of an airfoil by numerically solving the RANS equations. &
    &Optionally the heat transfer coefficient on airfoil surface can also be obtained. This program reads airfoil coordinates of grid points (or control points) &
    &from a 1D XYZ/CPT file, then generates grid and solution files in 2D PLOT3D format to record grid and airflow data and a report file to show aerodynamic &
    &parameters."
    print *,"More information can be found in the README file."
    print *,"Options:"
    print *,"--stdin        interactive mode, read parameters from standard input"
    print *,"--script       batch mode, read parameters from a configuration file"
    print *,"--version      show the current version"
    print *,"-h, --help     show this message"
   else
    print *,"This option does not exist."
   end if
  else if(argc==2) then
   Call get_command_argument(1,argv(2))
   Call get_command_argument(2,argv(3))
   if(trim(argv(2))=="--script") then
    mode='S'
    scptname=trim(argv(3))
    Call Aero2D(mode,1,scptname)
   else
    print *,"This option does not exist."
   end if
  end if
end Program RANSFOIL
