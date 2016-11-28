Module Aero2DCOM
use ISO_C_BINDING
implicit none
real(8),parameter::Pi=3.1415926535897932d+0,R=8.31,Ma=0.029,miu0=1.716e-5,Ti=273.11,Si=110.56,g=9.8,gama=1.4,Prt=0.85,kapa=0.4187,&
Ep=9.793,Cks=0.5,Cb1=0.1355,Cb2=0.622,sigman=2.0/3,Cw2=0.3,Cw3=2.0,Cniu1=7.1,Cw1=Cb1/kapa**2+(1+Cb2)/sigman,C1e=1.44,C2e=1.92,&
Cu=0.09,sigmak=1.0,sigmae=1.3,sigmak1=1.176,sigmak2=1.0,sigmaw1=2.0,sigmaw2=1.168,alpha1=0.31,betai1=0.075,betai2=0.0828,&
alphastarf=1.0,alpha0=1./9,betastarf=0.09,Rbeta=8.0,Rk=6.0,Rw=2.95,Zetastar=1.5,Mt0=0.25
save
integer(C_INT),bind(C)::Ip,Jp,Ic,Jc,Ib1,Ib2,Iwd,Iwu,Ifd,maxs
character(8) Proctrl,Energy,visheat,Turmod,Walltreat,solctrl,Discret,denface,Init,Stag
character(64) filename(9),dir
character(C_CHAR),bind(C)::cProctrl(8),cEnergy(8),cvisheat(8),cTurmod(8),cWalltreat(8),csolctrl(8),cDiscret(8),cdenface(8),&
cInit(8),cStag(8)
character(C_CHAR),bind(C)::cfilename(64),cdir(64)
real(C_DOUBLE),bind(C)::dx,dy,fd,delta,Rau,Rap,Rae,Rat,Vfar,AoA,Ta,Tf,Po,ksi,Itur,tvr,c,Ui,Vi,roui,miui,Tki,Tei,Twi,Tni,&
ca,ka,Re,Mach,rmsu,rmsv,rmst,rmsn,rmsk,rmse,rmsw,rmsm,Cl,Cd,Cf,Cm,Xpc,Ypc
real(8),allocatable,target,dimension(:,:)::rou,miu,P,dP,U,V,T,Tn,Tk,Te,Tw,miut,U0,V0,T0,Tn0,Tk0,Te0,Tw0,Pr,Pc,auP,auNB,aP,aW,&
aE,aS,aN,b,Xg,Yg,Xc,Yc,Xga,Xgk,Yga,Ygk,Jg,a1,y1,b1,d,Un,Vn,wdu,edu,sdv,ndv,Unw,Une,Vns,Vnn,Ux,Uy,Vx,Vy,Px,Py,dPx,dPy,roux,rouy,&
Tnx,Tny,Tkx,Tky,Twx,Twy,muxx,muxy,muyx,mvxy,mvyx,mvyy,ww,we,ws,wn,rouw,roue,rous,roun,Xi,fniu1,fai2,F2,St,Ret,alphastar,sigmatk,&
sigmatw
real(8),allocatable,target,dimension(:)::Xwd,Ywd,Xwu,Ywu,Xw,Yw,Yp,DR,Sw,ks,Q,Yplus,Ystar,ustar,Uplus,Tplus,hcv,Ax,Ay
real(C_DOUBLE),pointer::fXwd(:),fYwd(:),fXwu(:),fYwu(:)
type(C_PTR),bind(C)::cXwd,cYwd,cXwu,cYwu,cXw,cYw,cSw,cYplus,cYstar,chcv,cAx,cAy
type(C_PTR),bind(C)::cXg,cYg,cXc,cYc,crou,cmiu,cP,cVx,cVy,cT,cTn,cTk,cTe,cTw,cmiut

! Below is a list of common identifiers in aero2dcom module
! =========================================================
! characters variables
! identifier name        meaning
! Proctrl                control parameter, configure compressible property of air
! Energy                 control parameter, configure whether including energy equation
! visheat                control parameter, configure whether including the viscous heating terms
! Turmod                 control parameter, configure the turbulence model
! Walltreat              control parameter, configure the wall treatment method
! solctrl                control parameter, configure the algorithm of velocity-pressure coupling
! Discret                control parameter, configure the spatial schemes of convective terms
! denface                control parameter, configure the interpolation scheme of air density
! Init                   control parameter, configure whether using a data file for initialization
! Stag                   control parameter, configure whether using the stagnation values for initialization
! dir                    directory name of output files
! -----------------------------------------------------
! characters arrays
! identifier name        meaning
! filename               file name array
! filename(1)            aerofoil coordinates filename
! filename(9)            data filename for initialization
! --------------------------------------------------------
! integer variables
! identifier name        meaning
! Ip,Jp                  number of gird points in corresponding directions
! Ic,Jc                  number of gird cells in corresponding directions
! Ib1,Ib2                lower and upper index of gird cells on airfoil
! Iwd,Iwu                number of coordinate points on lower and upper airfoil
! Ifd                    layers of near wall mesh
! maxs                   maximum iteration steps
! ----------------------------------------------
! floating-point variables
! identifier name        meaning
! dx,dy                  spatial step size in corresponding directions of computational space
! fd                     dimensionless near wall mesh spacing
! delta                  minimum residual
! Rau                    relaxation factor of velocity
! Rap                    relaxation factor of pressure
! Rae                    relaxation factor of temperature
! Rat                    relaxation factor of turbulence
! Vfar                   free stream velocity
! AoA                    angle of attack
! Ta                     static temperature of free stream
! Tf                     wall temperature
! Po                     static pressure of free stream
! ksi                    equivalent sand grain roughness
! Itur                   turbulence intensity of free stream
! tvr                    turbulent viscosity ratio of free stream
! c                      chord length
! Cl                     lift coefficient
! Cd                     drag coefficient
! Cf                     friction coefficient
! Cm                     pitching moment coefficient
! Xpc,Ypc                X-Y coordinates of pressure center
! ---------------------------------------------------------
! one-dimensional floating-point arrays
! identifier name        meaning
! Xwd,Ywd,Xwu,Ywu        X-Y coordinates of grid points on lower and upper airfoil
! Xw,Yw                  X-Y coordinates of cells center on airfoil
! Sw                     curve length on airfoil
! Yplus,Ystar            dimensionless wall distance
! hcv                    surface heat transfer coefficient
! Ax,Ay                  X-Y components of wall shear stress
! ----------------------------------------------------------
! two-dimensional floating-point arrays
! identifier name        meaning
! Xg,Yg                  X-Y coordinates of grid points in structured mesh
! Xc,Yc                  X-Y coordinates of cells center in structured mesh
! rou                    density field
! miu                    molecular viscosity field
! P                      pressure field
! dP                     correctional pressure field
! U,V                    X-Y components of velocity field
! T                      temperature field
! Tn                     modified turbulent viscosity field
! Tk                     turbulence kinetic energy field
! Te                     turbulence dissipation rate field
! Tw                     turbulence specific dissipation rate field
! miut                   turbulent viscosity field

! The identifiers with prefix `c' is the C partners of the corresponding Fortran variables

end Module Aero2DCOM
