Subroutine hypresolve(aP,aW,aE,aS,aN,ba,F,F0,Ra,Ic,Jc,Ib1,Ib2,solid,scalar)
implicit none
include 'HYPREf.h'

integer  i,j,k,l,Ic,Jc,Ib1,Ib2,solid
integer  ierr,ndims,nentries,nparts,nvars,part,var,object_type,itmax,prlv,iter,nb
integer  ilower(2),iupper(2),stencil_indices(5),offsets(2,5),vartypes(1),bclower(2),bcupper(2),nblower(2),nbupper(2),map(2),dir(2)
real*8  Ra,tol,res
real*8  aP(Ic,Jc),aW(Ic,Jc),aE(Ic,Jc),aS(Ic,Jc),aN(Ic,Jc),ba(Ic,Jc),aR(Ic,Jc),br(Ic,Jc),F(Ic,Jc),F0(Ic,Jc)
real*8  values(5*Ic*Jc)
character(*) scalar

integer*8  grid
integer*8  stencil
integer*8  graph
integer*8  A
integer*8  b
integer*8  x
integer*8  parA
integer*8  parb
integer*8  parx
integer*8  solver

integer*8  MPI_COMM_WORLD

integer,parameter::HYPRE_SSTRUCT_VARIABLE_CELL = 0

DO j=1,Jc
 DO i=1,Ic
  if(i>1.and.i<Ic.and.j<Jc) then
   if(j==1.and.i>=Ib1.and.i<=Ib2) then
    if(scalar=='Te'.or.scalar=='Tw') then
    aR(i,j)=aP(i,j)
    br(i,j)=ba(i,j)
    else
    aR(i,j)=aP(i,j)/Ra
    br(i,j)=ba(i,j)+(1-Ra)*aP(i,j)*F0(i,j)/Ra
    end if
   else
    aR(i,j)=aP(i,j)/Ra
    br(i,j)=ba(i,j)+(1-Ra)*aP(i,j)*F0(i,j)/Ra
   end if
  else
   aR(i,j)=aP(i,j)
   br(i,j)=ba(i,j)
  end if
 end DO
end DO

ndims = 2
nparts = 1
nvars = 1
nentries = 5
part = 0
var = 0
object_type = HYPRE_PARCSR
vartypes(1) = HYPRE_SSTRUCT_VARIABLE_CELL

Call HYPRE_SStructGridCreate(MPI_COMM_WORLD,ndims,nparts,grid,ierr)
 ilower(1) = 1
 ilower(2) = 1
 iupper(1) = Ic
 iupper(2) = Jc
 Call HYPRE_SStructGridSetExtents(grid,part,ilower,iupper,ierr)
 Call HYPRE_SStructGridSetVariables(grid,part,nvars,vartypes,ierr)
 nb = 0
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
 l=1
 DO j=1,Jc
  DO i=1,Ic
   DO k=1,nentries
    if(k==1) then
    values(l)=aR(i,j)
    else if(k==2) then
    values(l)=-aW(i,j)
    else if(k==3) then
    values(l)=-aE(i,j)
    else if(k==4) then
    values(l)=-aS(i,j)
    else if(k==5) then
    values(l)=-aN(i,j)
    end if
    l=l+1
   end DO
  end DO
 end DO
 Call HYPRE_SStructMatrixSetBoxValues(A,part,ilower,iupper,var,nentries,stencil_indices,values,ierr)
Call HYPRE_SStructMatrixAssemble(A,ierr)
Call HYPRE_SStructMatrixGetObject(A, parA, ierr)

Call HYPRE_SStructVectorCreate(MPI_COMM_WORLD,grid,b,ierr)
Call HYPRE_SStructVectorCreate(MPI_COMM_WORLD,grid,x,ierr)
Call HYPRE_SStructVectorSetObjectTyp(b,object_type,ierr)
Call HYPRE_SStructVectorSetObjectTyp(x,object_type,ierr)
Call HYPRE_SStructVectorInitialize(b,ierr)
Call HYPRE_SStructVectorInitialize(x,ierr)
Call HYPRE_SStructVectorSetBoxValues(b,part,ilower,iupper,var,br,ierr)
Call HYPRE_SStructVectorSetBoxValues(x,part,ilower,iupper,var,F,ierr)
Call HYPRE_SStructVectorAssemble(b,ierr)
Call HYPRE_SStructVectorAssemble(x,ierr)
Call HYPRE_SStructVectorGetObject(b, parb, ierr)
Call HYPRE_SStructVectorGetObject(x, parx, ierr)

itmax = 1000
prlv = 0
if(scalar=='dP') then
tol = 1.0D-3
else
tol = 1.0D-8
end if

if(solid==1) then
Call HYPRE_SStructBiCGSTABCreate(MPI_COMM_WORLD, solver, ierr)
Call HYPRE_SStructBiCGSTABSetTol(solver, tol, ierr)
Call HYPRE_SStructBiCGSTABSetPrintLe(solver, prlv, ierr)
Call HYPRE_SStructBiCGSTABSetMaxIter(solver, itmax, ierr)
Call HYPRE_SStructBiCGSTABSetup(solver, A, b, x, ierr)
Call HYPRE_SStructBiCGSTABSolve(solver, A, b, x, ierr)
else if(solid==2) then
Call HYPRE_SStructPCGCreate(MPI_COMM_WORLD, solver, ierr)
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
!Call HYPRE_BoomerAMGSetCoarsenType(solver, 10, ierr)
!Call HYPRE_BoomerAMGSetInterpType(solver, 6, ierr)
!Call HYPRE_BoomerAMGSetRelaxType(solver, 6, ierr)
Call HYPRE_BoomerAMGSetup(solver, parA, parb, parx, ierr)
Call HYPRE_BoomerAMGSolve(solver, parA, parb, parx, ierr)
end if

Call HYPRE_SStructVectorGather(x,ierr)
Call HYPRE_SStructVectorGetBoxValues(x,part,ilower,iupper,var,F,ierr)

if(solid==1) then
 Call HYPRE_SStructBiCGSTABGetNumIter(solver, iter, ierr)
 Call HYPRE_SStructBiCGSTABGetFinalRe(solver, res, ierr)
else if(solid==2) then
 Call HYPRE_SStructPCGGetNumIteration(solver, iter, ierr)
 Call HYPRE_SStructPCGGetFinalRelativ(solver, res, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGGetNumIterations(solver, iter, ierr)
 Call HYPRE_BoomerAMGGetFinalReltvRes(solver, res, ierr)
end if

Call HYPRE_SStructGridDestroy(grid, ierr)
Call HYPRE_SStructStencilDestroy(stencil, ierr)
Call HYPRE_SStructGraphDestroy(graph, ierr)
Call HYPRE_SStructMatrixDestroy(A, ierr)
Call HYPRE_SStructVectorDestroy(b, ierr)
Call HYPRE_SStructVectorDestroy(x, ierr)
if(solid==1) then
Call HYPRE_SStructBiCGSTABDestroy(solver, ierr)
else if(solid==2) then
Call HYPRE_SStructPCGDestroy(solver, ierr)
else if(solid==3) then
Call HYPRE_BoomerAMGDestroy(solver, ierr)
end if

end Subroutine hypresolve
