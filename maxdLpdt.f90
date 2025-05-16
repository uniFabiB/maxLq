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
      ! iguess 9 = load temp &loadTempFunctionName
      ! iguess 50 = Arnold-Beltrami-Childress, ...
      !=============================================
      iguess = 50

      if(iguess==9) then
         loadTempFunctionName = "u_result_B32_0512.nc"
      end if

      bIterOffset = 0       ! should match the loaded iteration or 0 if new

      !lebesgueQlist = (/2.0, 4.0, 5.0, 7.0, 10.0/)
      lebesgueQ = 5.0_pr

      
      IF (rank==0) print*, "start"

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
      call createDirectoryIfNonExistent(HomeDir)
      call saveUsedParams()

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
      allocate( B_list(1:64+bIterOffset) )
      allocate( optimizationResultList(1:3,0:size(B_list)) )
      

      if(abs(lebesgueQ-4.0_pr)<MACH_EPSILON) then
         B_list(1) = 0.1_pr
         do B_list_iterator=2,size(B_list)
            constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/4.0_pr)
            if(constraintB>10.0_pr) constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/8.0_pr)
            if(constraintB>50.0_pr) constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/32.0_pr)
            B_list(B_list_iterator) = constraintB         
         end do
      elseif(abs(lebesgueQ-3.0_pr)<MACH_EPSILON)
         B_list(1) = 1.0_pr
         do B_list_iterator=2,size(B_list)
            constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/4.0_pr)
            if(constraintB>150.0_pr) constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/16.0_pr)
            if(constraintB>170.0_pr) constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/32.0_pr)
            B_list(B_list_iterator) = constraintB
         end do
      else
         B_list(1) = 1.0_pr
         do B_list_iterator=2,size(B_list)
            constraintB = B_list(B_list_iterator-1)*10**(1.0_pr/4.0_pr)
            B_list(B_list_iterator) = constraintB
         end do
      end if

      !=========================================================
      ! Loop over different values of constraint B values
      !=========================================================
      do B_list_iterator = 1+bIterOffset,size(B_list)+bIterOffset
         constraintB = B_list(B_list_iterator)
         write(bIterTxt, '(i2.2)') B_list_iterator
         write(bTxt, '(ES7.1)') constraintB
         if(rank==0) print*, "constraint B"//bIterTxt, constraintB
         !=========================================================
         ! OPTIMIZE !
         !=========================================================
         call maxdLqdt

      END DO

      call fft_deallocate
!      call fftw_mpi_cleanup()          ! Newly Added in Jan 17, but not sure whether it should be here; March 20, 2017
      CALL MPI_FINALIZE (Statinfo)
      

      IF (rank==0) print*, "end"
 
   END PROGRAM main
