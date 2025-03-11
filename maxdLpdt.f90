!-----------------------------------------------------!
! Program used to maximize dEdt                       !
! in the 3D Navier-Stokes system,                     !
! using 2 constraints.                                !
!                                                     !
! Parallel version, Complex-Complex FFT               !   
!                                                     !
! November, 2014.                                     !
!                                                     !
! Author: Diego Ayala                                 !
! Department of Mathematics and Statistics            !
! McMaster University                                 !
!-----------------------------------------------------!
 
   PROGRAM main

      USE global_variables
      USE data_ops
      use fftwfunction
      USE function_ops
      USE optimization
      USE mpi
      IMPLICIT NONE
      !=============================================
      ! Declare parameters
      !=============================================
      INTEGER :: RESOL

      

      !=============================================
      ! MPI Initialization; and fftw mpi initialization
      !=============================================
      CALL MPI_INIT(Statinfo)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,Statinfo)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,Statinfo)
      call fftw_mpi_init()                     ! Newly added in Jan 17, but not sure whether it should be here

      !=============================================
      ! Read in runtime parameters; input manually parameters' values.
      ! Read in code parameters; pass input values to code's parameters. 
      ! Newly added in Jan 23
      !=============================================

      !=============================================
      ! Initial Data
      ! iguess 0 = load FRT_N256E500T017_Uvec_fwdTE0220
      ! iguess 1 = load sines
      ! iguess 2 = load random a
      ! iguess 3 = load random smooth a
      ! iguess 4 = load random expSpec a
      ! iguess 5 = load random polySpec a
      ! iguess 6 = load random k a
      ! iguess 10 = load previous (not working)
      ! iguess 50 = Arnold-Beltrami-Childress, ...
      !=============================================
      iguess = 50
      RESOL = 32

      !lebesgueQlist = (/2.0, 4.0, 5.0, 7.0, 10.0/)
      lebesgueQ = 5.0_pr

      lambda1  = 1.0_pr  ! Newly added on Otc 05, 2017
      lambda2  = 0.0_pr  ! We are using H^s with 0<s<2 and do not have H^2, so we do not need this !OLD I changed this velue on March 5, 2017; lambda2 is the value in Sobolev norm


      IF (rank==0) THEN
         print*, "start"
      END If

      !=============================================
      !--If put discretization number before MPI_BARRIER, we have to MPI_BCAST n to all processors
      !=============================================
      C_n(1) = RESOL
      C_n(2) = RESOL
      C_n(3) = RESOL
      n(1) = RESOL
      n(2) = RESOL
      n(3) = RESOL

      !=============================================
      !--Initialize parameters values
      !=============================================
      CALL initialize()
      call init_fft()

      !=============================================
      ! Set initial Uvec
      !=============================================
      CALL initial_guess
         

      if (rank == 0 .and. ((abs(viscCoefficient-1.0)>MACH_EPSILON) .or. (abs(pressureCoefficient-1.0)>MACH_EPSILON))) then
         print*, "WARNING viscCoefficient ",viscCoefficient," or pressureCoefficient ", pressureCoefficient, "not 1"
      end if
      !=========================================================
      ! Create B value and result list
      !=========================================================
      allocate( B_list(0:3) )
      allocate( optimizationResultList(1:3,0:size(B_list)) )
      do B_list_iterator=0,size(B_list)-1
         B_list(B_list_iterator) = 10.0_pr**(2.0_pr+real(B_list_iterator,pr)/5.0_pr)
      end do

      !=========================================================
      ! Loop over different values of constraint B values
      !=========================================================
      do B_list_iterator = 0,size(B_list)-1
         constraintB = B_list(B_list_iterator)
         if(rank==0) print*, "constraint B", constraintB
         !=========================================================
         ! OPTIMIZE !
         !=========================================================
         call maxdLqdt

      END DO

      call fft_deallocate
!      call fftw_mpi_cleanup()          ! Newly Added in Jan 17, but not sure whether it should be here; March 20, 2017
      CALL MPI_FINALIZE (Statinfo)
      

      IF (rank==0) THEN
         print*, "end"
      END If
 
   END PROGRAM main
