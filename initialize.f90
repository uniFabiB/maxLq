SUBROUTINE initialize
   use, intrinsic :: iso_c_binding     ! Newly Added Jan 17
   USE global_variables
   use mpi
   IMPLICIT NONE
   INCLUDE "fftw3-mpi.f03"             ! Needed, as there is a fftw_mpi_local_size_3d command below

   INTEGER :: i,j,k
 
   C_local_alloc = fftw_mpi_local_size_3d(C_n(3), C_n(2), C_n(1), MPI_COMM_WORLD, C_local_N, C_local_k_offset)   ! Newly Added March 18, 2017, use "C_..."
   local_N = int( C_local_N )
   local_k_offset = int( C_local_k_offset )

   total_local_size = n(1)*n(2)*local_N
   
   ALLOCATE( Uvec(1:n(1),1:n(2),1:local_N,1:3) )
   ALLOCATE( Wvec(1:n(1),1:n(2),1:local_N,1:3) )

!   cxdata = fftw_alloc_complex(C_local_alloc)
!   cydata = fftw_alloc_complex(C_local_alloc)
!   czdata = fftw_alloc_complex(C_local_alloc)
!   call c_f_pointer(cxdata, fxdata, [C_n(1),C_n(2),C_local_N])
!   call c_f_pointer(cydata, fydata, [C_n(1),C_n(2),C_local_N])
!   call c_f_pointer(czdata, fzdata, [C_n(1),C_n(2),C_local_N])

   ALLOCATE( K1(1:n(1)) )
   ALLOCATE( K2(1:n(2)) )
   ALLOCATE( K3(1:n(3)) )
! ALLOCATE( fwork(0:total_local_size) )

  n_dim = 3*n(1)*n(2)*local_N
  dV = 1.0_pr/PRODUCT(REAL(n,pr))
!  Kcut = 2.0_pr*PI*REAL(n(1),pr)/3.0_pr

  Kcut = 2.0_pr*PI*1.0_pr*REAL(n(1),pr)/3.0_pr   ! I changed on Oct 8, 2017

  Kmax = 0.0_pr ! will be determined later when dealising
  testNonlinOrder = 0.0_pr
  if (rank == 0) then
      print*, "kmax_initially", PI*real(n(1),pr)
  end if

  !--Set up wavenumbers
  DO i = 0, n(1)-1
    IF (i<=n(1)/2) THEN
       K1(i+1) = 2.0_pr*PI*REAL(i,pr)
    ELSE
       K1(i+1) = 2.0_pr*PI*REAL(i-n(1),pr)
    END IF
  END DO
  ! K1(n(1)/2+1) = 0

  DO i = 0,n(2)-1
    IF (i <= n(2)/2) THEN
      K2(i+1) = 2.0_pr*PI*REAL(i,pr)
    ELSE
      K2(i+1) = 2.0_pr*PI*REAL(i-n(2),pr)
    END IF
  END DO
  ! K2(n(2)/2+1) = 0

  DO i = 0, n(3)-1
    IF (i<=n(3)/2) THEN
       K3(i+1) = 2.0_pr*PI*REAL(i,pr)
    ELSE
       K3(i+1) = 2.0_pr*PI*REAL(i-n(3),pr)
    END IF
  END DO
  ! K3(n(3)/2+1) = 0

  kappaTest = .true.
  toDealias = .true.
  add_pert = .FALSE.
  save_diag_NS = .FALSE.
  save_data_NS = .FALSE.
  save_diag_lineMin = .FALSE.
  save_data_lineMin = .FALSE.
  save_diag_Constr = .false.
  save_data_Constr = .FALSE.
  save_diag_Optim = .false.
  save_data_Optim = .false.
 
  IF (n(1)<256) THEN
     parallel_data = .FALSE.
  ELSE
     parallel_data = .TRUE.
  END IF
 
END SUBROUTINE


