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
   

   ALLOCATE( K1(1:n(1)) )
   ALLOCATE( K2(1:n(2)) )
   ALLOCATE( K3(1:n(3)) )

   dV = 1.0_pr/PRODUCT(REAL(n,pr))

   PI = 4.0_pr*ATAN2(1.0_pr,1.0_pr)            ! Set value of the constant pi

   kmax = PI*real(n(1),pr)       ! might be changed later when dealiasing depending on the powers 
   !if (rank == 0) print*, "kmax_initially", kmax

   !--Set up wavenumbers
   DO i = 0, n(1)-1
      IF (i<=n(1)/2) THEN
         K1(i+1) = 2.0_pr*PI*REAL(i,pr)
      ELSE
         K1(i+1) = 2.0_pr*PI*REAL(i-n(1),pr)
      END IF
   END DO

   DO i = 0,n(2)-1
      IF (i <= n(2)/2) THEN
      K2(i+1) = 2.0_pr*PI*REAL(i,pr)
      ELSE
      K2(i+1) = 2.0_pr*PI*REAL(i-n(2),pr)
      END IF
   END DO

   DO i = 0, n(3)-1
      IF (i<=n(3)/2) THEN
         K3(i+1) = 2.0_pr*PI*REAL(i,pr)
      ELSE
         K3(i+1) = 2.0_pr*PI*REAL(i-n(3),pr)
      END IF
   END DO


   !IF (n(1)<256) THEN     !todo temporarily changed on 16 may 2025 to avoid an error of u1 being not allocatable in save_field_R3toR3_ncdf : ALLOCATE( u1(1:n(1),1:n(2),1:n(3)) ) 
   IF (n(1)<64) THEN
      parallel_data = .FALSE.
   ELSE
      parallel_data = .TRUE.
   END IF
 
END SUBROUTINE


