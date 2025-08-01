SUBROUTINE initialize
   use, intrinsic :: iso_c_binding     ! Newly Added Jan 17
   USE global_variables
   use mpi
   IMPLICIT NONE
   INCLUDE "fftw3-mpi.f03"             ! Needed, as there is a fftw_mpi_local_size_3d command below

   CHARACTER(len=100) :: tempHostName
   logical :: hostFound = .false.
   INTEGER :: i, myIndex




   call hostnm(tempHostName)
   hostName = trim(tempHostName)
   myIndex = index(string=hostName,substring="gra")
   if(myIndex==1) then
      if(rank==0) print*, "hostName ", hostName, " -> ", "graham"
      inputDir = "/home/fabianbl/projects/rrg-bprotas/fabianbl/prog/input/"
      if(rank==0) print*, achar(9), "using input dir ", inputDir
      hostFound = .true.
   end if

   myIndex = index(string=hostName,substring=".nibi.sharcnet")
   if(myIndex>0) then
      if(hostFound .and. rank==0) print*, "ERROR: HOSTNAME FOUND NIBI BUT PREVIOUSLY FOUND ANOTHER ONE"
      if(rank==0) print*, "hostName ", hostName, " -> ", "nibi"
      inputDir = "/home/fabianbl/projects/rrg-bprotas/fabianbl/prog/input/"
      if(rank==0) print*, achar(9), "using input dir ", inputDir
      hostFound = .true.
   end if

   myIndex = index(string=hostName,substring=".int.cedar.computecanada.ca")
   if(myIndex>0) then
      if(hostFound .and. rank==0) print*, "ERROR: HOSTNAME FOUND CEDAR BUT PREVIOUSLY FOUND ANOTHER ONE"
      if(rank==0) print*, "hostName ", hostName, " -> ", "cedar"
      inputDir = "/home/fabianbl/projects/def-bprotas/fabianbl/prog/input/"
      if(rank==0) print*, achar(9), "using input dir ", inputDir
      hostFound = .true.      
   end if

   if(.not. hostFound) then
      if(rank==0) print*, "WARNING: NO HOSTNAME COULD BE ATTRIBUTED"
      inputDir = "/home/fabianbl/projects/def-bprotas/fabianbl/prog/input/"
      if(rank==0) print*, achar(9), "using standard input directory ", inputDir
   end if
   
   C_n(1) = RESOL
   C_n(2) = RESOL
   C_n(3) = RESOL
   n(1) = RESOL
   n(2) = RESOL
   n(3) = RESOL
   
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

   if(use_e_u_auto_for_q_less_4) then
      if(lebesgueQ<4.0_pr-MACH_EPSILON) then
         use_e_u_instead_of_uqMinus4 = .true.
      else
         use_e_u_instead_of_uqMinus4 = .false.
      end if
   end if
   


   !IF (n(1)<256) THEN     !todo temporarily changed on 16 may 2025 to avoid an error of u1 being not allocatable in save_field_R3toR3_ncdf : ALLOCATE( u1(1:n(1),1:n(2),1:n(3)) ) 
   IF (n(1)<16) THEN
      parallel_data = .FALSE.
   ELSE
      parallel_data = .TRUE.
   END IF
 
END SUBROUTINE


