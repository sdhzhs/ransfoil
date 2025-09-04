Subroutine hypreinit_gpu(A,b,x,solver,precond,p_values,Ic,Jc,Ib1,Ib2,solid)
use, intrinsic :: iso_c_binding
use, intrinsic :: iso_fortran_env, only: int64
use cudaf

implicit none
include 'HYPREf.h'

integer      i,Ic,Jc,Ib1,Ib2
integer(1)   solid
integer      ierr,ndims,nentries,nparts,nvars,part,var,object_type,nb
integer      ilower(2),iupper(2),stencil_indices(5),offsets(2,5),vartypes(1),bclower(2),bcupper(2),nblower(2),nbupper(2),map(2),dir(2)

integer(8)  grid
integer(8)  stencil
integer(8)  graph
integer(8)  A
integer(8)  b
integer(8)  x
integer(8)  solver,precond

integer(8)  MPI_COMM_WORLD
integer,parameter::HYPRE_SSTRUCT_VARIABLE_CELL = 0

type(c_ptr) :: p_values

integer :: stat

!stat = device_malloc_managed(int(40*Ic*Jc, int64), p_values)
stat = device_malloc(int(40*Ic*Jc, int64), p_values)

ndims = 2
nparts = 1
nvars = 1
nentries = 5
part = 0
var = 0
object_type = HYPRE_PARCSR
vartypes(1) = HYPRE_SSTRUCT_VARIABLE_CELL

Call HYPRE_Initialize(ierr)
Call HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE,ierr)
Call HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE,ierr)
Call HYPRE_SetSpGemmUseVendor(0,ierr)

Call HYPRE_SStructGridCreate(MPI_COMM_WORLD,ndims,nparts,grid,ierr)
 ilower(1) = 1
 ilower(2) = 1
 iupper(1) = Ic
 iupper(2) = Jc
 Call HYPRE_SStructGridSetExtents(grid,part,ilower,iupper,ierr)
 Call HYPRE_SStructGridSetVariables(grid,part,nvars,vartypes,ierr)
 nb = 0
 if(Ib1>1.and.Ib2<Ic) then
  bclower(1)=1
  bclower(2)=0
  bcupper(1)=Ib1-1
  bcupper(2)=0
  nblower(1)=Ic
  nblower(2)=1
  nbupper(1)=Ib2+1
  nbupper(2)=1
  map(1)=0
  map(2)=1
  dir(1)=-1
  dir(2)=-1
  Call HYPRE_SStructGridSetNeighborPart(grid,part,bclower,bcupper,nb,nblower,nbupper,map,dir,ierr)
  bclower(1)=Ib2+1
  bcupper(1)=Ic
  nblower(1)=Ib1-1
  nbupper(1)=1
  Call HYPRE_SStructGridSetNeighborPart(grid,part,bclower,bcupper,nb,nblower,nbupper,map,dir,ierr)
 else
  bclower(1)=0
  bclower(2)=1
  bcupper(1)=0
  bcupper(2)=Jc
  nblower(1)=Ic
  nblower(2)=1
  nbupper(1)=Ic
  nbupper(2)=Jc
  map(1)=0
  map(2)=1
  dir(1)=1
  dir(2)=1
  Call HYPRE_SStructGridSetNeighborPart(grid,part,bclower,bcupper,nb,nblower,nbupper,map,dir,ierr)
  bclower(1)=Ic+1
  bcupper(1)=Ic+1
  nblower(1)=1
  nbupper(1)=1
  Call HYPRE_SStructGridSetNeighborPart(grid,part,bclower,bcupper,nb,nblower,nbupper,map,dir,ierr)
 end if
Call HYPRE_SStructGridAssemble(grid,ierr)

Call HYPRE_SStructStencilCreate(ndims,nentries,stencil,ierr)
offsets(1,1) = 0
offsets(2,1) = 0
offsets(1,2) = -1
offsets(2,2) = 0
offsets(1,3) = 1
offsets(2,3) = 0
offsets(1,4) = 0
offsets(2,4) = -1
offsets(1,5) = 0
offsets(2,5) = 1
DO i=1,nentries
 stencil_indices(i) = i-1
end DO
DO i=1,nentries
 Call HYPRE_SStructStencilSetEntry(stencil,stencil_indices(i),offsets(1,i),var,ierr)
end DO

Call HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, graph, ierr)
Call HYPRE_SStructGraphSetObjectType(graph,object_type,ierr)
Call HYPRE_SStructGraphSetStencil(graph,part,var,stencil,ierr)
Call HYPRE_SStructGraphAssemble(graph,ierr)

Call HYPRE_SStructMatrixCreate(MPI_COMM_WORLD,graph,A,ierr)
Call HYPRE_SStructMatrixSetObjectTyp(A,object_type,ierr)
Call HYPRE_SStructMatrixInitialize(A,ierr)

Call HYPRE_SStructVectorCreate(MPI_COMM_WORLD,grid,b,ierr)
Call HYPRE_SStructVectorCreate(MPI_COMM_WORLD,grid,x,ierr)
Call HYPRE_SStructVectorSetObjectTyp(b,object_type,ierr)
Call HYPRE_SStructVectorSetObjectTyp(x,object_type,ierr)
Call HYPRE_SStructVectorInitialize(b,ierr)
Call HYPRE_SStructVectorInitialize(x,ierr)

if(solid==1) then
 Call HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, solver, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGCreate(MPI_COMM_WORLD, solver, ierr)
else if(solid==3) then
 !Call HYPRE_BoomerAMGCreate(solver, ierr)
else if(solid==4) then
 Call HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, solver, ierr)
 !Call HYPRE_BoomerAMGCreate(precond, ierr)
end if

Call HYPRE_SStructGridDestroy(grid, ierr)
Call HYPRE_SStructStencilDestroy(stencil, ierr)
Call HYPRE_SStructGraphDestroy(graph, ierr)

end Subroutine hypreinit_gpu

Subroutine hyprerelease_gpu(A,b,x,solver,precond,p_values,solid)
use, intrinsic :: iso_c_binding
use cudaf

implicit none
include 'HYPREf.h'

integer(1)   solid
integer      ierr

integer(8)  A
integer(8)  b
integer(8)  x
integer(8)  solver,precond

type(c_ptr) :: p_values

integer :: stat

if(solid==1) then
 Call HYPRE_SStructBiCGSTABDestroy(solver, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGDestroy(solver, ierr)
else if(solid==3) then
 !Call HYPRE_BoomerAMGDestroy(solver, ierr)
else if(solid==4) then
 !Call HYPRE_BoomerAMGDestroy(precond, ierr)
 Call HYPRE_ParCSRBiCGSTABDestroy(solver, ierr)
end if
Call HYPRE_SStructMatrixDestroy(A, ierr)
Call HYPRE_SStructVectorDestroy(b, ierr)
Call HYPRE_SStructVectorDestroy(x, ierr)

Call HYPRE_Finalize(ierr)

stat = device_free(p_values)

end Subroutine hyprerelease_gpu

Subroutine hyprecompute_gpu(A,b,x,solver,precond,p_values,aM,ba,F,F0,Ra,Ic,Jc,Ib1,Ib2,solid,scalar)
use, intrinsic :: iso_c_binding
use, intrinsic :: iso_fortran_env, only: int64
use cudaf

implicit none
include 'HYPREf.h'

integer      i,j,Ic,Jc,Ib1,Ib2,Is,Ie
integer(1)   solid
integer      ierr,nentries,part,var,itmax,prlv,iter,precond_id
integer      ilower(2),iupper(2),stencil_indices(5)
real(8)      Ra,tol,res
real(8), target:: aM(5,Ic,Jc),ba(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
character(*) scalar
logical(1) isP,isT,isTe,isTw

integer(8)  A
integer(8)  b
integer(8)  x
integer(8)  parA
integer(8)  parb
integer(8)  parx
integer(8)  solver,precond

real(8), pointer :: values(:)
type(c_ptr) :: p_values, p_F, p_aM, p_ba

integer :: stat

isP = scalar=='dP'
isT = scalar=='T'
isTe = scalar=='Te'
isTw = scalar=='Tw'

Call c_f_pointer(p_values, values, [5*Ic*Jc])
p_F = C_LOC(F)
p_aM = C_LOC(aM)
p_ba = C_LOC(ba)

if(Ib1>1.and.Ib2<Ic) then
 Is=2
 Ie=Ic-1
else
 Is=1
 Ie=Ic
end if

DO j=1,Jc
 DO i=1,Ic
  if(i>=Is.and.i<=Ie.and.j<Jc) then
   if(.not.(j==1.and.i>=Ib1.and.i<=Ib2.and.(isTe.or.isTw))) then
    ba(i,j)=ba(i,j)+(1-Ra)*aM(1,i,j)*F0(i,j)/Ra
    aM(1,i,j)=aM(1,i,j)/Ra
    aM(2,i,j)=-aM(2,i,j)
    aM(3,i,j)=-aM(3,i,j)
    aM(4,i,j)=-aM(4,i,j)
    aM(5,i,j)=-aM(5,i,j)
   end if
  end if
 end DO
end DO

stat = copy_host_to_device(int(Ic*Jc * 40, int64), p_values, p_aM)

nentries = 5
part = 0
var = 0
ilower(1) = 1
ilower(2) = 1
iupper(1) = Ic
iupper(2) = Jc

DO i=1,nentries
 stencil_indices(i) = i-1
end DO

Call HYPRE_SStructMatrixSetBoxValues(A,part,ilower,iupper,var,nentries,stencil_indices,values,ierr)

Call HYPRE_SStructMatrixAssemble(A,ierr)
Call HYPRE_SStructMatrixGetObject(A, parA, ierr)

stat = copy_host_to_device(int(Ic*Jc * 8, int64), p_values, p_ba)

Call HYPRE_SStructVectorSetBoxValues(b,part,ilower,iupper,var,values,ierr)

stat = copy_host_to_device(int(Ic*Jc * 8, int64), p_values, p_F)

Call HYPRE_SStructVectorSetBoxValues(x,part,ilower,iupper,var,values,ierr)

Call HYPRE_SStructVectorAssemble(b,ierr)
Call HYPRE_SStructVectorAssemble(x,ierr)
Call HYPRE_SStructVectorGetObject(b, parb, ierr)
Call HYPRE_SStructVectorGetObject(x, parx, ierr)

itmax = 1000
prlv = 0
if(isP) then
 tol = 1.0e-4
else if(isT) then
 tol = 1.0e-8
else
 tol = 1.0e-6
end if

if(solid==1) then
 Call HYPRE_SStructBiCGSTABSetTol(solver, tol, ierr)
 Call HYPRE_SStructBiCGSTABSetPrintLe(solver, prlv, ierr)
 Call HYPRE_SStructBiCGSTABSetMaxIter(solver, itmax, ierr)
 Call HYPRE_SStructBiCGSTABSetup(solver, A, b, x, ierr)
 Call HYPRE_SStructBiCGSTABSolve(solver, A, b, x, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGSetTol(solver, tol, ierr)
 Call HYPRE_SStructPCGSetPrintLevel(solver, prlv, ierr)
 Call HYPRE_SStructPCGSetMaxIter(solver, itmax, ierr)
 Call HYPRE_SStructPCGSetup(solver, A, b, x, ierr)
 Call HYPRE_SStructPCGSolve(solver, A, b, x, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGCreate(solver, ierr)
 Call HYPRE_BoomerAMGSetTol(solver, tol, ierr)
 Call HYPRE_BoomerAMGSetPrintLevel(solver, prlv, ierr)
 Call HYPRE_BoomerAMGSetMaxIter(solver, itmax, ierr)
 if(isP) then
  Call HYPRE_BoomerAMGSetMaxLevels(solver, 20, ierr)
 else
  Call HYPRE_BoomerAMGSetMaxLevels(solver, 1, ierr)
 end if
 !Call HYPRE_BoomerAMGSetCoarsenType(solver, 8, ierr)
 !Call HYPRE_BoomerAMGSetInterpType(solver, 6, ierr)
 !Call HYPRE_BoomerAMGSetRelaxType(solver, 6, ierr)
 !Call HYPRE_BoomerAMGSetRelaxOrder(solver, 0, ierr)
 !Call HYPRE_BoomerAMGSetAggNumLevels(solver, 5, ierr)
 !Call HYPRE_BoomerAMGSetAggInterpType(solver, 5, ierr)
 !Call HYPRE_BoomerAMGSetKeepTransp(solver, 1, ierr)
 !Call HYPRE_BoomerAMGSetRAP2(solver, 0, ierr)
 !Call HYPRE_BoomerAMGSetNumSweeps(solver, 1, ierr)
 !Call HYPRE_BoomerAMGSetSmoothType(solver, 9, ierr)
 !if(isP) then
  !Call HYPRE_BoomerAMGSetSmoothNumLvls(solver, 20, ierr)
 !else
  !Call HYPRE_BoomerAMGSetSmoothNumLvls(solver, 1, ierr)
 !end if
 !Call HYPRE_BoomerAMGSetSmoothNumSwps(solver, 1, ierr)
 Call HYPRE_BoomerAMGSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_BoomerAMGSolve(solver, parA, parb, parx, ierr)
else if(solid==4) then
 Call HYPRE_ParCSRBiCGSTABSetTol(solver, tol, ierr)
 Call HYPRE_ParCSRBiCGSTABSetPrintLev(solver, prlv, ierr)
 Call HYPRE_ParCSRBiCGSTABSetMaxIter(solver, itmax, ierr)
 
 Call HYPRE_BoomerAMGCreate(precond, ierr)
 Call HYPRE_BoomerAMGSetPrintLevel(precond, 0, ierr)
 Call HYPRE_BoomerAMGSetTol(precond, 0e+0, ierr)
 Call HYPRE_BoomerAMGSetMaxIter(precond, 1, ierr)
 if(isP) then
  Call HYPRE_BoomerAMGSetMaxLevels(precond, 20, ierr)
 else
  Call HYPRE_BoomerAMGSetMaxLevels(precond, 1, ierr)
 end if
 !Call HYPRE_BoomerAMGSetNumSweeps(precond, 1, ierr)
 !Call HYPRE_BoomerAMGSetCoarsenType(precond, 6, ierr)
 !Call HYPRE_BoomerAMGSetInterpType(precond, 0, ierr)
 !Call HYPRE_BoomerAMGSetRelaxType(precond, 6, ierr)
 
 precond_id = 2
 Call HYPRE_ParCSRBiCGSTABSetPrecond(solver, precond_id, precond, ierr)
 Call HYPRE_ParCSRBiCGSTABSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_ParCSRBiCGSTABSolve(solver, parA, parb, parx, ierr)
end if

Call HYPRE_SStructVectorGather(x,ierr)
Call HYPRE_SStructVectorGetBoxValues(x,part,ilower,iupper,var,values,ierr)

stat = copy_device_to_host(int(Ic*Jc * 8, int64), p_F, p_values)

if(solid==1) then
 Call HYPRE_SStructBiCGSTABGetNumIter(solver, iter, ierr)
 Call HYPRE_SStructBiCGSTABGetFinalRe(solver, res, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGGetNumIteration(solver, iter, ierr)
 Call HYPRE_SStructPCGGetFinalRelativ(solver, res, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGGetNumIterations(solver, iter, ierr)
 Call HYPRE_BoomerAMGGetFinalReltvRes(solver, res, ierr)
 Call HYPRE_BoomerAMGDestroy(solver, ierr)
else if(solid==4) then
 Call HYPRE_ParCSRBiCGSTABGetNumIter(solver, iter, ierr)
 Call HYPRE_ParCSRBiCGSTABGetFinalRel(solver, res, ierr)
 Call HYPRE_BoomerAMGDestroy(precond, ierr)
end if

end Subroutine hyprecompute_gpu

Subroutine hypresolve_gpu(aM,ba,F,F0,Ra,Ic,Jc,Ib1,Ib2,solid,scalar)
use, intrinsic :: iso_c_binding
use, intrinsic :: iso_fortran_env, only: int64
use cudaf

implicit none
include 'HYPREf.h'

integer      i,j,Ic,Jc,Ib1,Ib2,Is,Ie
integer(1)   solid
integer      ierr,ndims,nentries,nparts,nvars,part,var,object_type,nb,itmax,prlv,iter,precond_id
integer      ilower(2),iupper(2),stencil_indices(5),offsets(2,5),vartypes(1),bclower(2),bcupper(2),nblower(2),nbupper(2),map(2),dir(2)
real(8)      Ra,tol,res
real(8), target:: aM(5,Ic,Jc),ba(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
character(*) scalar
logical(1) isP,isT,isTe,isTw

integer(8)  grid
integer(8)  stencil
integer(8)  graph
integer(8)  A
integer(8)  b
integer(8)  x
integer(8)  parA
integer(8)  parb
integer(8)  parx
integer(8)  solver,precond

integer(8)  MPI_COMM_WORLD
integer,parameter::HYPRE_SSTRUCT_VARIABLE_CELL = 0

real(8), pointer :: values(:)
type(c_ptr) :: p_values, p_F, p_aM, p_ba

integer :: stat

isP = scalar=='dP'
isT = scalar=='T'
isTe = scalar=='Te'
isTw = scalar=='Tw'

!stat = device_malloc_managed(int(5*Ic*Jc * 8, int64), p_values)
stat = device_malloc(int(5*Ic*Jc * 8, int64), p_values)
Call c_f_pointer(p_values, values, [5*Ic*Jc])
p_F = C_LOC(F)
p_aM = C_LOC(aM)
p_ba = C_LOC(ba)

ndims = 2
nparts = 1
nvars = 1
nentries = 5
part = 0
var = 0
object_type = HYPRE_PARCSR
vartypes(1) = HYPRE_SSTRUCT_VARIABLE_CELL

Call HYPRE_Initialize(ierr)
Call HYPRE_SetMemoryLocation(HYPRE_MEMORY_DEVICE,ierr)
Call HYPRE_SetExecutionPolicy(HYPRE_EXEC_DEVICE,ierr)
Call HYPRE_SetSpGemmUseVendor(0,ierr)

Call HYPRE_SStructGridCreate(MPI_COMM_WORLD,ndims,nparts,grid,ierr)
 ilower(1) = 1
 ilower(2) = 1
 iupper(1) = Ic
 iupper(2) = Jc
 Call HYPRE_SStructGridSetExtents(grid,part,ilower,iupper,ierr)
 Call HYPRE_SStructGridSetVariables(grid,part,nvars,vartypes,ierr)
 nb = 0
 if(Ib1>1.and.Ib2<Ic) then
  bclower(1)=1
  bclower(2)=0
  bcupper(1)=Ib1-1
  bcupper(2)=0
  nblower(1)=Ic
  nblower(2)=1
  nbupper(1)=Ib2+1
  nbupper(2)=1
  map(1)=0
  map(2)=1
  dir(1)=-1
  dir(2)=-1
  Call HYPRE_SStructGridSetNeighborPart(grid,part,bclower,bcupper,nb,nblower,nbupper,map,dir,ierr)
  bclower(1)=Ib2+1
  bcupper(1)=Ic
  nblower(1)=Ib1-1
  nbupper(1)=1
  Call HYPRE_SStructGridSetNeighborPart(grid,part,bclower,bcupper,nb,nblower,nbupper,map,dir,ierr)
 else
  bclower(1)=0
  bclower(2)=1
  bcupper(1)=0
  bcupper(2)=Jc
  nblower(1)=Ic
  nblower(2)=1
  nbupper(1)=Ic
  nbupper(2)=Jc
  map(1)=0
  map(2)=1
  dir(1)=1
  dir(2)=1
  Call HYPRE_SStructGridSetNeighborPart(grid,part,bclower,bcupper,nb,nblower,nbupper,map,dir,ierr)
  bclower(1)=Ic+1
  bcupper(1)=Ic+1
  nblower(1)=1
  nbupper(1)=1
  Call HYPRE_SStructGridSetNeighborPart(grid,part,bclower,bcupper,nb,nblower,nbupper,map,dir,ierr)
 end if
Call HYPRE_SStructGridAssemble(grid,ierr)

Call HYPRE_SStructStencilCreate(ndims,nentries,stencil,ierr)
offsets(1,1) = 0
offsets(2,1) = 0
offsets(1,2) = -1
offsets(2,2) = 0
offsets(1,3) = 1
offsets(2,3) = 0
offsets(1,4) = 0
offsets(2,4) = -1
offsets(1,5) = 0
offsets(2,5) = 1
DO i=1,nentries
 stencil_indices(i) = i-1
end DO
DO i=1,nentries
 Call HYPRE_SStructStencilSetEntry(stencil,stencil_indices(i),offsets(1,i),var,ierr)
end DO

Call HYPRE_SStructGraphCreate(MPI_COMM_WORLD, grid, graph, ierr)
Call HYPRE_SStructGraphSetObjectType(graph,object_type,ierr)
Call HYPRE_SStructGraphSetStencil(graph,part,var,stencil,ierr)
Call HYPRE_SStructGraphAssemble(graph,ierr)

Call HYPRE_SStructMatrixCreate(MPI_COMM_WORLD,graph,A,ierr)
Call HYPRE_SStructMatrixSetObjectTyp(A,object_type,ierr)
Call HYPRE_SStructMatrixInitialize(A,ierr)

Call HYPRE_SStructVectorCreate(MPI_COMM_WORLD,grid,b,ierr)
Call HYPRE_SStructVectorCreate(MPI_COMM_WORLD,grid,x,ierr)
Call HYPRE_SStructVectorSetObjectTyp(b,object_type,ierr)
Call HYPRE_SStructVectorSetObjectTyp(x,object_type,ierr)
Call HYPRE_SStructVectorInitialize(b,ierr)
Call HYPRE_SStructVectorInitialize(x,ierr)

if(solid==1) then
 Call HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, solver, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGCreate(MPI_COMM_WORLD, solver, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGCreate(solver, ierr)
else if(solid==4) then
 Call HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, solver, ierr)
 Call HYPRE_BoomerAMGCreate(precond, ierr)
end if

if(Ib1>1.and.Ib2<Ic) then
 Is=2
 Ie=Ic-1
else
 Is=1
 Ie=Ic
end if

DO j=1,Jc
 DO i=1,Ic
  if(i>=Is.and.i<=Ie.and.j<Jc) then
   if(.not.(j==1.and.i>=Ib1.and.i<=Ib2.and.(isTe.or.isTw))) then
    ba(i,j)=ba(i,j)+(1-Ra)*aM(1,i,j)*F0(i,j)/Ra
    aM(1,i,j)=aM(1,i,j)/Ra
    aM(2,i,j)=-aM(2,i,j)
    aM(3,i,j)=-aM(3,i,j)
    aM(4,i,j)=-aM(4,i,j)
    aM(5,i,j)=-aM(5,i,j)
   end if
  end if
 end DO
end DO

stat = copy_host_to_device(int(Ic*Jc * 40, int64), p_values, p_aM)

Call HYPRE_SStructMatrixSetBoxValues(A,part,ilower,iupper,var,nentries,stencil_indices,values,ierr)

Call HYPRE_SStructMatrixAssemble(A,ierr)
Call HYPRE_SStructMatrixGetObject(A, parA, ierr)

stat = copy_host_to_device(int(Ic*Jc * 8, int64), p_values, p_ba)

Call HYPRE_SStructVectorSetBoxValues(b,part,ilower,iupper,var,values,ierr)

stat = copy_host_to_device(int(Ic*Jc * 8, int64), p_values, p_F)

Call HYPRE_SStructVectorSetBoxValues(x,part,ilower,iupper,var,values,ierr)

Call HYPRE_SStructVectorAssemble(b,ierr)
Call HYPRE_SStructVectorAssemble(x,ierr)
Call HYPRE_SStructVectorGetObject(b, parb, ierr)
Call HYPRE_SStructVectorGetObject(x, parx, ierr)

itmax = 1000
prlv = 0
if(isP) then
 tol = 1.0e-4
else if(isT) then
 tol = 1.0e-8
else
 tol = 1.0e-6
end if

if(solid==1) then
 Call HYPRE_SStructBiCGSTABSetTol(solver, tol, ierr)
 Call HYPRE_SStructBiCGSTABSetPrintLe(solver, prlv, ierr)
 Call HYPRE_SStructBiCGSTABSetMaxIter(solver, itmax, ierr)
 Call HYPRE_SStructBiCGSTABSetup(solver, A, b, x, ierr)
 Call HYPRE_SStructBiCGSTABSolve(solver, A, b, x, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGSetTol(solver, tol, ierr)
 Call HYPRE_SStructPCGSetPrintLevel(solver, prlv, ierr)
 Call HYPRE_SStructPCGSetMaxIter(solver, itmax, ierr)
 Call HYPRE_SStructPCGSetup(solver, A, b, x, ierr)
 Call HYPRE_SStructPCGSolve(solver, A, b, x, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGSetTol(solver, tol, ierr)
 Call HYPRE_BoomerAMGSetPrintLevel(solver, prlv, ierr)
 Call HYPRE_BoomerAMGSetMaxIter(solver, itmax, ierr)
 if(isP) then
  Call HYPRE_BoomerAMGSetMaxLevels(solver, 20, ierr)
 else
  Call HYPRE_BoomerAMGSetMaxLevels(solver, 1, ierr)
 end if
 !Call HYPRE_BoomerAMGSetCoarsenType(solver, 8, ierr)
 !Call HYPRE_BoomerAMGSetInterpType(solver, 6, ierr)
 !Call HYPRE_BoomerAMGSetRelaxType(solver, 6, ierr)
 !Call HYPRE_BoomerAMGSetRelaxOrder(solver, 0, ierr)
 !Call HYPRE_BoomerAMGSetAggNumLevels(solver, 5, ierr)
 !Call HYPRE_BoomerAMGSetAggInterpType(solver, 5, ierr)
 !Call HYPRE_BoomerAMGSetKeepTransp(solver, 1, ierr)
 !Call HYPRE_BoomerAMGSetRAP2(solver, 0, ierr)
 !Call HYPRE_BoomerAMGSetNumSweeps(solver, 1, ierr)
 !Call HYPRE_BoomerAMGSetSmoothType(solver, 9, ierr)
 !if(isP) then
  !Call HYPRE_BoomerAMGSetSmoothNumLvls(solver, 20, ierr)
 !else
  !Call HYPRE_BoomerAMGSetSmoothNumLvls(solver, 1, ierr)
 !end if
 !Call HYPRE_BoomerAMGSetSmoothNumSwps(solver, 1, ierr)
 Call HYPRE_BoomerAMGSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_BoomerAMGSolve(solver, parA, parb, parx, ierr)
else if(solid==4) then
 Call HYPRE_ParCSRBiCGSTABSetTol(solver, tol, ierr)
 Call HYPRE_ParCSRBiCGSTABSetPrintLev(solver, prlv, ierr)
 Call HYPRE_ParCSRBiCGSTABSetMaxIter(solver, itmax, ierr)
 
 Call HYPRE_BoomerAMGSetPrintLevel(precond, 0, ierr)
 Call HYPRE_BoomerAMGSetTol(precond, 0e+0, ierr)
 Call HYPRE_BoomerAMGSetMaxIter(precond, 1, ierr)
 if(isP) then
  Call HYPRE_BoomerAMGSetMaxLevels(precond, 20, ierr)
 else
  Call HYPRE_BoomerAMGSetMaxLevels(precond, 1, ierr)
 end if
 !Call HYPRE_BoomerAMGSetNumSweeps(precond, 1, ierr)
 !Call HYPRE_BoomerAMGSetCoarsenType(precond, 6, ierr)
 !Call HYPRE_BoomerAMGSetInterpType(precond, 0, ierr)
 !Call HYPRE_BoomerAMGSetRelaxType(precond, 6, ierr)
 
 precond_id = 2
 Call HYPRE_ParCSRBiCGSTABSetPrecond(solver, precond_id, precond, ierr)
 Call HYPRE_ParCSRBiCGSTABSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_ParCSRBiCGSTABSolve(solver, parA, parb, parx, ierr)
end if

Call HYPRE_SStructVectorGather(x,ierr)
Call HYPRE_SStructVectorGetBoxValues(x,part,ilower,iupper,var,values,ierr)

stat = copy_device_to_host(int(Ic*Jc * 8, int64), p_F, p_values)

if(solid==1) then
 Call HYPRE_SStructBiCGSTABGetNumIter(solver, iter, ierr)
 Call HYPRE_SStructBiCGSTABGetFinalRe(solver, res, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGGetNumIteration(solver, iter, ierr)
 Call HYPRE_SStructPCGGetFinalRelativ(solver, res, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGGetNumIterations(solver, iter, ierr)
 Call HYPRE_BoomerAMGGetFinalReltvRes(solver, res, ierr)
else if(solid==4) then
 Call HYPRE_ParCSRBiCGSTABGetNumIter(solver, iter, ierr)
 Call HYPRE_ParCSRBiCGSTABGetFinalRel(solver, res, ierr)
end if

if(solid==1) then
 Call HYPRE_SStructBiCGSTABDestroy(solver, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGDestroy(solver, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGDestroy(solver, ierr)
else if(solid==4) then
 Call HYPRE_BoomerAMGDestroy(precond, ierr)
 Call HYPRE_ParCSRBiCGSTABDestroy(solver, ierr)
end if
Call HYPRE_SStructGridDestroy(grid, ierr)
Call HYPRE_SStructStencilDestroy(stencil, ierr)
Call HYPRE_SStructGraphDestroy(graph, ierr)
Call HYPRE_SStructMatrixDestroy(A, ierr)
Call HYPRE_SStructVectorDestroy(b, ierr)
Call HYPRE_SStructVectorDestroy(x, ierr)

Call HYPRE_Finalize(ierr)

stat = device_free(p_values)

end Subroutine hypresolve_gpu
