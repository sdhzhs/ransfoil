Subroutine hyprecpinit(A_ij, b_ij, x_ij, solver, precond, solid)
use Aero2DCOM
implicit none
include 'HYPREf.h'

integer :: solid
integer(8) :: MPI_COMM_WORLD
integer(8) :: A_ij, b_ij, x_ij
integer(8) :: solver, precond

integer :: ilower_ij, iupper_ij, ierr

ilower_ij = 0
iupper_ij = 3 * Ic * Jc - 1

Call HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower_ij, iupper_ij, ilower_ij, iupper_ij, A_ij, ierr)
Call HYPRE_IJMatrixSetObjectType(A_ij, HYPRE_PARCSR, ierr)
Call HYPRE_IJMatrixInitialize(A_ij, ierr)

Call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_ij, iupper_ij, b_ij, ierr)
Call HYPRE_IJVectorSetObjectType(b_ij, HYPRE_PARCSR, ierr)
Call HYPRE_IJVectorInitialize(b_ij, ierr)

Call HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower_ij, iupper_ij, x_ij, ierr)
Call HYPRE_IJVectorSetObjectType(x_ij, HYPRE_PARCSR, ierr)
Call HYPRE_IJVectorInitialize(x_ij, ierr)

if(solid==1) then
 Call HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, solver, ierr)
else if(solid==2) then
 Call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGCreate(solver, ierr)
else if(solid==4) then
 Call HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, solver, ierr)
 !Call HYPRE_BoomerAMGCreate(precond, ierr)
else if(solid==5) then
 Call HYPRE_ParCSRBiCGSTABCreate(MPI_COMM_WORLD, solver, ierr)
 Call HYPRE_EuclidCreate(MPI_COMM_WORLD, precond, ierr)
else if(solid==6) then
 Call HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, solver, ierr)
 Call HYPRE_EuclidCreate(MPI_COMM_WORLD, precond, ierr)
end if

end Subroutine hyprecpinit

Subroutine hyprecprelease(A_ij, b_ij, x_ij, solver, precond, solid)
use Aero2DCOM
implicit none
include 'HYPREf.h'

integer :: solid
integer(8) :: A_ij, b_ij, x_ij
integer(8) :: solver, precond
integer :: ierr

Call HYPRE_IJMatrixDestroy(A_ij, ierr)
Call HYPRE_IJVectorDestroy(b_ij, ierr)
Call HYPRE_IJVectorDestroy(x_ij, ierr)

if(solid==1) then
 Call HYPRE_ParCSRBiCGSTABDestroy(solver, ierr)
else if(solid==2) then
 Call HYPRE_ParCSRGMRESDestroy(solver, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGDestroy(solver, ierr)
else if(solid==4) then
 !Call HYPRE_BoomerAMGDestroy(precond, ierr)
 Call HYPRE_ParCSRBiCGSTABDestroy(solver, ierr)
else if(solid==5) then
 Call HYPRE_EuclidDestroy(precond, ierr)
 Call HYPRE_ParCSRBiCGSTABDestroy(solver, ierr)
 else if(solid==6) then
 Call HYPRE_EuclidDestroy(precond, ierr)
 Call HYPRE_ParCSRGMRESDestroy(solver, ierr)
end if

end Subroutine hyprecprelease

Subroutine hyprecpsolve(A_ij, b_ij, x_ij, solver, precond, solid)
use Aero2DCOM
implicit none
include 'HYPREf.h'

integer :: solid
integer :: i,j
integer :: ierr,itmax,prlv,iter,maxl,precond_id
real(8) :: tol,res
real(8) :: aR(Ic,Jc)
real(8),allocatable,dimension(:) :: values
integer,allocatable,dimension(:) :: cols
character(32) fn

real(8) :: vec_values(3)
integer :: vec_indices(3)

integer :: cell_idx, global_row, neighbor_idx, ncols
integer :: nrows, ncols_arr(1), rows(1)
integer :: nvalues

integer(8) :: A_ij, b_ij, x_ij
integer(8) :: parA, parb, parx
integer(8) :: solver, precond

! Setup RHS and aR arrays
DO j=1,Jc
 DO i=1,Ic
 if((Ib1==1.or.(i>1.and.i<Ic)).and.j<Jc) then
  aR(i,j)=auM(1,i,j)/Rau
  bu(i,j)=bu(i,j)+(1-Rau)*auM(1,i,j)*U0(i,j)/Rau
  bv(i,j)=bv(i,j)+(1-Rau)*auM(6,i,j)*V0(i,j)/Rau
  auM(6,i,j)=auM(6,i,j)/Rau
  b(i,j)=b(i,j)+(1-Rap)*aM(1,i,j)*P(i,j)/Rap
  aM(1,i,j)=aM(1,i,j)/Rap
 else
  aR(i,j)=auM(1,i,j)
 end if
 end DO
end DO

allocate(values(15))
allocate(cols(15))

nrows = 1

DO j=1,Jc
 DO i=1,Ic
  cell_idx = (j-1)*Ic + (i-1)
  
  ! ==========================================
  ! U-Momentum Equation (var = 0)
  ! ==========================================
  global_row = 3*cell_idx + 0
  rows(1) = global_row
  ncols = 0
  
  ! Center U
  ncols = ncols + 1
  cols(ncols) = 3*cell_idx + 0
  values(ncols) = aR(i,j)
  
  ! West U
  if(i>1) then
   neighbor_idx = (j-1)*Ic + (i-2)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -auM(2,i,j)
  end if
  
  ! East U
  if(i<Ic) then
   neighbor_idx = (j-1)*Ic + i
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -auM(3,i,j)
  end if

  if(Ib1==1.and.i==1) then
   neighbor_idx = (j-1)*Ic + (Ic-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -auM(2,i,j)
  end if

  if(Ib1==1.and.i==Ic) then
   neighbor_idx = (j-1)*Ic
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -auM(3,i,j)
  end if
  
  ! South U
  if(j>1) then
   neighbor_idx = (j-2)*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -auM(4,i,j)
  end if

  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   neighbor_idx = (j-1)*Ic + (Ic-i)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -auM(4,i,j)
  end if
  
  ! North U
  if(j<Jc) then
   neighbor_idx = j*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -auM(5,i,j)
  end if
  
  ! Center P
  ncols = ncols + 1
  cols(ncols) = 3*cell_idx + 2
  values(ncols) = aupM(1,i,j)
  
  ! West P
  if(i>1) then
   neighbor_idx = (j-1)*Ic + (i-2)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aupM(2,i,j)
  end if
  
  ! East P
  if(i<Ic) then
   neighbor_idx = (j-1)*Ic + i
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aupM(3,i,j)
  end if

  if(Ib1==1.and.i==1) then
   neighbor_idx = (j-1)*Ic + (Ic-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aupM(2,i,j)
  end if

  if(Ib1==1.and.i==Ic) then
   neighbor_idx = (j-1)*Ic
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aupM(3,i,j)
  end if
  
  ! South P
  if(j>1) then
   neighbor_idx = (j-2)*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aupM(4,i,j)
  end if

  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   neighbor_idx = (j-1)*Ic + (Ic-i)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aupM(4,i,j)
  end if
  
  ! North P
  if(j<Jc) then
   neighbor_idx = j*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aupM(5,i,j)
  end if
  
  ncols_arr(1) = ncols
  Call HYPRE_IJMatrixSetValues(A_ij, nrows, ncols_arr, rows, cols, values, ierr)

  ! ==========================================
  ! V-Momentum Equation (var = 1)
  ! ==========================================
  global_row = 3*cell_idx + 1
  rows(1) = global_row
  ncols = 0
  
  ! Center V
  ncols = ncols + 1
  cols(ncols) = 3*cell_idx + 1
  values(ncols) = auM(6,i,j)
  
  ! West V
  if(i>1) then
   neighbor_idx = (j-1)*Ic + (i-2)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -auM(2,i,j)
  end if
  
  ! East V
  if(i<Ic) then
   neighbor_idx = (j-1)*Ic + i
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -auM(3,i,j)
  end if

  if(Ib1==1.and.i==1) then
   neighbor_idx = (j-1)*Ic + (Ic-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -auM(2,i,j)
  end if

  if(Ib1==1.and.i==Ic) then
   neighbor_idx = (j-1)*Ic
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -auM(3,i,j)
  end if
  
  ! South V
  if(j>1) then
   neighbor_idx = (j-2)*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -auM(4,i,j)
  end if

  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   neighbor_idx = (j-1)*Ic + (Ic-i)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -auM(4,i,j)
  end if
  
  ! North V
  if(j<Jc) then
   neighbor_idx = j*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -auM(5,i,j)
  end if
  
  ! Center P
  ncols = ncols + 1
  cols(ncols) = 3*cell_idx + 2
  values(ncols) = avpM(1,i,j)
  
  ! West P
  if(i>1) then
   neighbor_idx = (j-1)*Ic + (i-2)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -avpM(2,i,j)
  end if
  
  ! East P
  if(i<Ic) then
   neighbor_idx = (j-1)*Ic + i
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -avpM(3,i,j)
  end if

  if(Ib1==1.and.i==1) then
   neighbor_idx = (j-1)*Ic + (Ic-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -avpM(2,i,j)
  end if

  if(Ib1==1.and.i==Ic) then
   neighbor_idx = (j-1)*Ic
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -avpM(3,i,j)
  end if
  
  ! South P
  if(j>1) then
   neighbor_idx = (j-2)*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -avpM(4,i,j)
  end if

  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   neighbor_idx = (j-1)*Ic + (Ic-i)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -avpM(4,i,j)
  end if
  
  ! North P
  if(j<Jc) then
   neighbor_idx = j*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -avpM(5,i,j)
  end if
  
  ncols_arr(1) = ncols
  Call HYPRE_IJMatrixSetValues(A_ij, nrows, ncols_arr, rows, cols, values, ierr)

  ! ==========================================
  ! Pressure Equation (var = 2)
  ! ==========================================
  global_row = 3*cell_idx + 2
  rows(1) = global_row
  ncols = 0
  
  ! Center U
  ncols = ncols + 1
  cols(ncols) = 3*cell_idx + 0
  values(ncols) = apuM(1,i,j)
  
  ! West U
  if(i>1) then
   neighbor_idx = (j-1)*Ic + (i-2)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -apuM(2,i,j)
  end if
  
  ! East U
  if(i<Ic) then
   neighbor_idx = (j-1)*Ic + i
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -apuM(3,i,j)
  end if

  if(Ib1==1.and.i==1) then
   neighbor_idx = (j-1)*Ic + (Ic-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -apuM(2,i,j)
  end if

  if(Ib1==1.and.i==Ic) then
   neighbor_idx = (j-1)*Ic
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -apuM(3,i,j)
  end if
  
  ! South U
  if(j>1) then
   neighbor_idx = (j-2)*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -apuM(4,i,j)
  end if

  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   neighbor_idx = (j-1)*Ic + (Ic-i)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -apuM(4,i,j)
  end if
  
  ! North U
  if(j<Jc) then
   neighbor_idx = j*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 0
   values(ncols) = -apuM(5,i,j)
  end if

  ! Center V
  ncols = ncols + 1
  cols(ncols) = 3*cell_idx + 1
  values(ncols) = apvM(1,i,j)
  
  ! West V
  if(i>1) then
   neighbor_idx = (j-1)*Ic + (i-2)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -apvM(2,i,j)
  end if
  
  ! East V
  if(i<Ic) then
   neighbor_idx = (j-1)*Ic + i
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -apvM(3,i,j)
  end if

  if(Ib1==1.and.i==1) then
   neighbor_idx = (j-1)*Ic + (Ic-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -apvM(2,i,j)
  end if

  if(Ib1==1.and.i==Ic) then
   neighbor_idx = (j-1)*Ic
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -apvM(3,i,j)
  end if
  
  ! South V
  if(j>1) then
   neighbor_idx = (j-2)*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -apvM(4,i,j)
  end if

  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   neighbor_idx = (j-1)*Ic + (Ic-i)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -apvM(4,i,j)
  end if
  
  ! North V
  if(j<Jc) then
   neighbor_idx = j*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 1
   values(ncols) = -apvM(5,i,j)
  end if

  ! Center P
  ncols = ncols + 1
  cols(ncols) = 3*cell_idx + 2
  values(ncols) = aM(1,i,j)
  
  ! West P
  if(i>1) then
   neighbor_idx = (j-1)*Ic + (i-2)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aM(2,i,j)
  end if
  
  ! East P
  if(i<Ic) then
   neighbor_idx = (j-1)*Ic + i
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aM(3,i,j)
  end if

  if(Ib1==1.and.i==1) then
   neighbor_idx = (j-1)*Ic + (Ic-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aM(2,i,j)
  end if

  if(Ib1==1.and.i==Ic) then
   neighbor_idx = (j-1)*Ic
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aM(3,i,j)
  end if
  
  ! South P
  if(j>1) then
   neighbor_idx = (j-2)*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aM(4,i,j)
  end if

  if(j==1.and.(i<Ib1.or.i>Ib2)) then
   neighbor_idx = (j-1)*Ic + (Ic-i)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aM(4,i,j)
  end if
  
  ! North P
  if(j<Jc) then
   neighbor_idx = j*Ic + (i-1)
   ncols = ncols + 1
   cols(ncols) = 3*neighbor_idx + 2
   values(ncols) = -aM(5,i,j)
  end if
  
  ncols_arr(1) = ncols
  Call HYPRE_IJMatrixSetValues(A_ij, nrows, ncols_arr, rows, cols, values, ierr)

 END DO
END DO

Call HYPRE_IJMatrixAssemble(A_ij, ierr)
Call HYPRE_IJMatrixGetObject(A_ij, parA, ierr)
deallocate(values)
deallocate(cols)

fn="debug_matrix"//char(0)
!Call HYPRE_ParCSRMatrixPrint(parA, fn, len_trim(fn), ierr)

nvalues = 3
DO j=1,Jc
 DO i=1,Ic
  cell_idx = (j-1)*Ic + (i-1)
  vec_indices(1) = 3*cell_idx + 0
  vec_indices(2) = 3*cell_idx + 1
  vec_indices(3) = 3*cell_idx + 2
  
  vec_values(1) = bu(i,j)
  vec_values(2) = bv(i,j)
  vec_values(3) = b(i,j)
  Call HYPRE_IJVectorSetValues(b_ij, nvalues, vec_indices, vec_values, ierr)
  
  vec_values(1) = 0
  vec_values(2) = 0
  vec_values(3) = 0
  Call HYPRE_IJVectorSetValues(x_ij, nvalues, vec_indices, vec_values, ierr)
 END DO
END DO

Call HYPRE_IJVectorAssemble(b_ij, ierr)
Call HYPRE_IJVectorGetObject(b_ij, parb, ierr)

Call HYPRE_IJVectorAssemble(x_ij, ierr)
Call HYPRE_IJVectorGetObject(x_ij, parx, ierr)

itmax = 1000
prlv = 0
tol = 1.0D-6

if(solid==1) then
 Call HYPRE_ParCSRBiCGSTABSetTol(solver, tol, ierr)
 Call HYPRE_ParCSRBiCGSTABSetPrintLev(solver, prlv, ierr)
 Call HYPRE_ParCSRBiCGSTABSetMaxIter(solver, itmax, ierr)
 Call HYPRE_ParCSRBiCGSTABSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_ParCSRBiCGSTABSolve(solver, parA, parb, parx, ierr)
else if(solid==2) then
 Call HYPRE_ParCSRGMRESSetTol(solver, tol, ierr)
 Call HYPRE_ParCSRGMRESSetPrintLevel(solver, prlv, ierr)
 Call HYPRE_ParCSRGMRESSetMaxIter(solver, itmax, ierr)
 Call HYPRE_ParCSRGMRESSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_ParCSRGMRESSolve(solver, parA, parb, parx, ierr)
else if(solid==3) then
 !Call HYPRE_BoomerAMGCreate(solver, ierr)
 Call HYPRE_BoomerAMGSetTol(solver, tol, ierr)
 Call HYPRE_BoomerAMGSetPrintLevel(solver, prlv, ierr)
 Call HYPRE_BoomerAMGSetMaxIter(solver, itmax, ierr)

 maxl = 2
 Call HYPRE_BoomerAMGSetMaxLevels(solver, maxl, ierr)
 Call HYPRE_BoomerAMGSetNumFunctions(solver, 3, ierr)
 Call HYPRE_BoomerAMGSetNodal(solver, 1, ierr)
 Call HYPRE_BoomerAMGSetNumSweeps(solver, 3, ierr)
 Call HYPRE_BoomerAMGSetSmoothType(solver, 6, ierr) !Schwarz method as smoother
 !Call HYPRE_BoomerAMGSetSmoothNumLvls(solver, maxl, ierr)
 Call HYPRE_BoomerAMGSetDomainType(solver, 1, ierr)
 Call HYPRE_BoomerAMGSetOverlap(solver, 0, ierr)
 !Call HYPRE_BoomerAMGSetCoarsenType(solver, 6, ierr)
 Call HYPRE_BoomerAMGSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_BoomerAMGSolve(solver, parA, parb, parx, ierr)
else if(solid==4) then
 Call HYPRE_ParCSRBiCGSTABSetTol(solver, tol, ierr)
 Call HYPRE_ParCSRBiCGSTABSetPrintLev(solver, prlv, ierr)
 Call HYPRE_ParCSRBiCGSTABSetMaxIter(solver, itmax, ierr)
 
 precond_id = 2
 maxl = 5
 Call HYPRE_BoomerAMGCreate(precond, ierr)
 Call HYPRE_BoomerAMGSetTol(precond, 0d+0, ierr)
 Call HYPRE_BoomerAMGSetPrintLevel(precond, 0, ierr)
 Call HYPRE_BoomerAMGSetMaxIter(precond, 1, ierr)
 Call HYPRE_BoomerAMGSetMaxLevels(precond, maxl, ierr)
 Call HYPRE_BoomerAMGSetNumFunctions(precond, 3, ierr)
 Call HYPRE_BoomerAMGSetNodal(precond, 1, ierr)
 Call HYPRE_BoomerAMGSetNumSweeps(precond, 3, ierr)
 Call HYPRE_BoomerAMGSetSmoothType(precond, 5, ierr) !ILU method as smoother
 Call HYPRE_BoomerAMGSetSmoothNumLvls(precond, maxl, ierr)
 !Call HYPRE_BoomerAMGSetDomainType(precond, 1, ierr)
 !Call HYPRE_BoomerAMGSetOverlap(precond, 0, ierr)
 !Call HYPRE_BoomerAMGSetCoarsenType(precond, 6, ierr)
 
 Call HYPRE_ParCSRBiCGSTABSetPrecond(solver, precond_id, precond, ierr)
 Call HYPRE_ParCSRBiCGSTABSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_ParCSRBiCGSTABSolve(solver, parA, parb, parx, ierr)
else if(solid==5) then
 Call HYPRE_ParCSRBiCGSTABSetTol(solver, tol, ierr)
 Call HYPRE_ParCSRBiCGSTABSetPrintLev(solver, prlv, ierr)
 Call HYPRE_ParCSRBiCGSTABSetMaxIter(solver, itmax, ierr)
 
 precond_id = 5
 Call HYPRE_EuclidSetLevel(precond, 0, ierr)
 Call HYPRE_EuclidSetBJ(precond, 1, ierr)
 
 Call HYPRE_ParCSRBiCGSTABSetPrecond(solver, precond_id, precond, ierr)
 Call HYPRE_ParCSRBiCGSTABSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_ParCSRBiCGSTABSolve(solver, parA, parb, parx, ierr)
else if(solid==6) then
 Call HYPRE_ParCSRGMRESSetTol(solver, tol, ierr)
 Call HYPRE_ParCSRGMRESSetPrintLevel(solver, prlv, ierr)
 Call HYPRE_ParCSRGMRESSetMaxIter(solver, itmax, ierr)
 Call HYPRE_ParCSRGMRESSetKDim(solver, 50, ierr)
 
 precond_id = 5
 Call HYPRE_EuclidSetLevel(precond, 0, ierr)
 Call HYPRE_EuclidSetBJ(precond, 1, ierr)
 
 Call HYPRE_ParCSRGMRESSetPrecond(solver, precond_id, precond, ierr)
 Call HYPRE_ParCSRGMRESSetup(solver, parA, parb, parx, ierr)
 Call HYPRE_ParCSRGMRESSolve(solver, parA, parb, parx, ierr)
end if

DO j=1,Jc
 DO i=1,Ic
  cell_idx = (j-1)*Ic + (i-1)
  vec_indices(1) = 3*cell_idx + 0
  vec_indices(2) = 3*cell_idx + 1
  vec_indices(3) = 3*cell_idx + 2
  Call HYPRE_IJVectorGetValues(x_ij, nvalues, vec_indices, vec_values, ierr)
  U(i,j) = vec_values(1)
  V(i,j) = vec_values(2)
  P(i,j) = vec_values(3)
 END DO
END DO

if(solid==1 .or. solid==4 .or. solid==5) then
 Call HYPRE_ParCSRBiCGSTABGetNumIter(solver, iter, ierr)
 Call HYPRE_ParCSRBiCGSTABGetFinalRel(solver, res, ierr)
else if(solid==2 .or. solid==6) then
 Call HYPRE_ParCSRGMRESGetNumIteratio(solver, iter, ierr)
 Call HYPRE_ParCSRGMRESGetFinalRelati(solver, res, ierr)
else if(solid==3) then
 Call HYPRE_BoomerAMGGetNumIterations(solver, iter, ierr)
 Call HYPRE_BoomerAMGGetFinalReltvRes(solver, res, ierr)
end if

!print *, 'Hypre coupled solver iterations: ', iter, ' Final residual: ', res

if(solid==3) then
 !Call HYPRE_BoomerAMGDestroy(solver, ierr)
else if(solid==4) then
 Call HYPRE_BoomerAMGDestroy(precond, ierr)
end if

end Subroutine hyprecpsolve
