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
 
   PROGRAM maxdEdtHeli_main

      USE global_variables
      USE data_ops
      use fftwfunction
      USE function_ops
      USE optimization
      USE mpi
      IMPLICIT NONE
      !INCLUDE "mpif.h"    ! rather use USE mpi before implicit none
      !include 'fftw3-mpi.f03'                 ! Newly added in Jan 17, but not sure whether it should be here

      !=============================================
      ! Declare parameters
      !=============================================
      INTEGER :: RESOL, ii, kk, E0index, E1index, Kindex, Nuindex, numPts_E0, tempk
      REAL(pr) :: aux
      REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: K0_vec!, E0_vec      !E0_vec is the array of the different enstrophy values
      CHARACTER(2) :: K0txt, E0txt, E1txt, IGtxt, NUtxt
      !integer, dimension (:), allocatable :: E0List = (/20/)               !E0List = (/1, 10, 20, 30, 40, 50, 60/)
      !integer, allocatable :: E0List(:)
      integer, dimension(:), allocatable :: E0List
      !integer :: E0List(*) = (/20/)
      REAL(pr), DIMENSION (60), save :: E0_vec

      
      E0List = (/20/)      !(/1, 10, 20, 30, 40, 50, 60/)
      lebesgueQlist = (/2.0, 4.0, 5.0, 7.0, 10.0/)
      !lebesgueQlist = (/3.0/)
      
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
      
      RESOL   = 64
      Kindex  = 0
      E0index = 25      ! does not get used right now, instead the E0List is used 
      E1index = 35      ! does not get used right now, instead the E0List is used
      Nuindex = 2

      IF (rank==0) THEN
         print*, "start"
      END If

      WRITE(K0txt, '(i2.2)') Kindex
      WRITE(E0txt, '(i2.2)') E0index
      WRITE(E1txt, '(i2.2)') E1index
      WRITE(NUtxt, '(i2.2)') Nuindex
      IF (Kindex==0) THEN                         ! Initialize the ConsType values for E0
         ConsType = 1
      ELSE
         ConsType = 2
      END IF
      PI = 4.0_pr*ATAN2(1.0_pr,1.0_pr)            ! Set value of the constant pi

!=============================================
! Create K0_vec and E0_vec
!============================================= 
      SELECT CASE (ConsType)
         CASE (1)
            K0 = 0.0_pr
            do tempk = 1, size(E0List)
               ii = E0List(tempk)
               !DO ii=1,E1index
               IF (ii .LE. 2) THEN
                  aux = REAL(ii-4,pr)
                  E0_vec(ii) = 1.0_pr*10**(aux)   ! 0.001, 0.01
               ELSEIF (ii .LE. 11) THEN
                  aux = REAL(ii-2,pr)
                  E0_vec(ii) = 0.1_pr*aux         ! 0.1, 0.2, ..., 0.9
               ELSEIF (ii .LE. 31) THEN
                  aux = REAL(ii-11,pr)
                  E0_vec(ii) = aux                ! 1, 2, ..., 20
               ELSEIF (ii .LE. 47) THEN
                  aux = REAL(ii-27,pr)
                  E0_vec(ii) = 5.0_pr*aux         ! 25, 30, ..., 100
               ELSEIF (ii .LE. 56) THEN
                  aux = REAL(ii-46,pr)
                  E0_vec(ii) = 100.0_pr*aux        ! 200, 300, ..., 1000
               ELSEIF (ii .LE. 60) THEN
                  aux = REAL(ii-55,pr)
                  E0_vec(ii) = 1000.0_pr*aux      ! 2000, 3000, ..., 10000
               else
                  print*, "E0index", ii, ">60, i.e. out of programmed choices"
                  error stop "E0index > 60, i.e. out of programmed choices in main program"
               END IF
            END DO

         CASE (2)
            numPts_E0 = 12*(2**Kindex-1) + 1
            ALLOCATE(K0_vec(1:4))
            !allocate(E0_vec(1:8))
            DO ii=1,4 
               K0_vec(ii) = 10.0_pr**(ii-1)
            END DO
            K0 = K0_vec(Kindex)
            DO ii=1,numPts_E0
               aux = REAL(ii-1,pr)/REAL(numPts_E0-1, pr) 
               E0_vec(ii) = 1.1_pr*(2.0_pr*PI)**2.0_pr*K0*(10.0_pr**aux)
            END DO
      END SELECT

!=========================================================
! Loop over different values of E0
!=========================================================
      !DO kk = E0index, E1index
      do tempk = 1, size(E0List)
         kk = E0List(tempk)
         WRITE(E0txt,'(I2.2)') kk
         E0_index = kk

!=============================================
! Set parameters' values manually (not read from files).
! Newly added in Jan 23
!=============================================
         K0_index = 0


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
         iguess = 5

         NU_index = 2
         lambda1  = 5.0_pr  ! Newly added on Otc 05, 2017
         lambda2  = 0.1_pr   ! I changed this velue on March 5, 2017; lambda2 is the value in Sobolev norm
         alpha0   = 1E-4


         viscCoefficient = 1.0_pr                                    ! for testing to turn off individual terms arising from Delta u
         pressureCoefficient = 1.0_pr                                ! for testing to turn off individual terms arising from nabla p

         if (rank == 0 .and. ((abs(viscCoefficient-1.0)>MACH_EPSILON) .or. (abs(pressureCoefficient-1.0)>MACH_EPSILON))) then
            print*, "WARNING viscCoefficient ",viscCoefficient," or pressureCoefficient ", pressureCoefficient, "not 1"
         end if

         CALL MPI_BARRIER(MPI_COMM_WORLD,Statinfo)                                ! Why put MPI_BARRIER here? to wait for what?
         CALL MPI_BCAST (RESOL,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)   ! Why MPI_BCAST these values as they are global?
         CALL MPI_BCAST (K0_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)
         CALL MPI_BCAST (E0_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)
         CALL MPI_BCAST (alpha0,   1, MPI_REAL8,   0, MPI_COMM_WORLD, Statinfo)
         CALL MPI_BCAST (lambda2,  1, MPI_REAL8,   0, MPI_COMM_WORLD, Statinfo)
         CALL MPI_BCAST (lambda1,  1, MPI_REAL8,   0, MPI_COMM_WORLD, Statinfo)   ! Newly added on Otc 05, 2017
         CALL MPI_BCAST (iguess,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)
         CALL MPI_BCAST (NU_index, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, Statinfo)

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
         E0 = E0_vec(E0_index)
         IF (rank==0) THEN
            print*, "E0_index = ", E0_index, "  E0 = ", E0
         END If

!=========================================================
!--Define value of viscosity
!=========================================================
         IF (NU_index==0) THEN
            visc = 0.0_pr
         ELSE
            visc = 10.0_pr**REAL(-Nuindex,pr)
         END IF
         WRITE(IGtxt,'(i2.2)') iguess

         CALL initial_guess

         if (.false.) then
            call testFFT(uvec(:,:,:,1))
         end if

   
         call maxdLqdt

!--Deallocate vectors
         DEALLOCATE(Uvec)
         DEALLOCATE(Wvec)
         DEALLOCATE(K1)
         DEALLOCATE(K2)
         DEALLOCATE(K3)
      END DO

      call fft_deallocate
!      call fftw_mpi_cleanup()          ! Newly Added in Jan 17, but not sure whether it should be here; March 20, 2017
      CALL MPI_FINALIZE (Statinfo)
      

      IF (rank==0) THEN
         print*, "end"
      END If
 
   END PROGRAM maxdEdtHeli_main
 
