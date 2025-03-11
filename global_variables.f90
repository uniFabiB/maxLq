MODULE global_variables
  use, intrinsic :: iso_c_binding   ! Newly Added Jan 17, and the following C_INTPTR_T, type(C_PTR) and so on
  IMPLICIT NONE

  INTEGER, PARAMETER :: pr = KIND (0.0d0)
  INTEGER, PARAMETER :: MAX_ITER = 1000                   !original 1000 ! Maximal iterations of maxdEdt
  INTEGER, PARAMETER :: MAX_ITER_CONSTR = 100
  REAL(pr), PARAMETER :: OPTIM_TOL = 1.0e-10_pr            !original 1.0e-8_pr 
  REAL(pr), PARAMETER :: CONSTR_TOL = 1.0e-10_pr
  REAL(pr), PARAMETER :: MACH_EPSILON = 2.0e-16_pr
  REAL(pr), PARAMETER :: J_MAX = 1.0e15_pr

  CHARACTER(len=*), parameter :: HomeDir = "./output/"
  CHARACTER(len=:), allocatable :: ConstraintDir

  logical :: verboseOptimization = .false.                 ! verbosely output stuff to terminal
  logical :: tauDebugToConsole = .true.                 ! verbosely output stuff to terminal


  LOGICAL :: kappaTest = .false.
  LOGICAL :: toDealias = .true.
  LOGICAL :: mnbra_calcSaveAllJvalues = .true.              ! calculates J(u+tau d) for "all" tau values to get an idea of the shape of J(tau) 
  LOGICAL :: save_diag_NS = .true.
  LOGICAL :: save_data_NS = .true.
  LOGICAL :: save_diag_Constr = .true.
  LOGICAL :: save_data_Constr = .true.
  LOGICAL :: save_diag_Optim = .true.
  LOGICAL :: save_data_Optim = .true.
  LOGICAL :: save_diag_lineMin = .true.
  LOGICAL :: save_data_lineMin = .true.
  LOGICAL :: add_pert = .false.
  LOGICAL :: parallel_data

  INTEGER, DIMENSION(3), SAVE :: n
  INTEGER, SAVE :: K0_index, E0_index, NU_index, iguess, ConsType
  REAL(pr), SAVE :: E0, K0, PI, visc, lambda2, lambda1, alpha0, dV, Kmax, lebesgueQ
  real(pr), dimension(:), allocatable, save :: lebesgueQlist
  real(pr), dimension(:,:), allocatable :: optimizationResultList
  real(pr), dimension(:), allocatable :: B_list
  integer :: B_list_iterator

  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: K1, K2, K3
  REAL(pr), DIMENSION (:,:,:,:), ALLOCATABLE, SAVE :: Uvec, Wvec

  real(pr) :: viscCoefficient = 1.0_pr                                          ! for debugging to turn on and off the viscocity/pressure terms
  real(pr) :: pressureCoefficient = 1.0_pr                                      ! for debugging to turn on and off the viscocity/pressure terms
  real(pr), save :: constraintB                                                 ! constraint size, i.e. ||u||_q = B
  
  !========================================================================== 
  !                            MPI VARIABLES
  !==========================================================================
  INTEGER(C_INTPTR_T), DIMENSION(3), SAVE :: C_n
  INTEGER, SAVE :: rank, Statinfo, np
  INTEGER(C_INTPTR_T), SAVE :: C_local_alloc, C_local_N, C_local_k_offset       ! Newly Added Jan 17; Modified on March 18, 2017
  INTEGER, SAVE :: local_alloc, local_N, local_k_offset                         ! Newly Added March 18, 2017

  INTEGER(C_INTPTR_T), SAVE :: local_last_start_after_trans, C_total_local_size
  INTEGER, SAVE :: total_local_size                                             ! Newly Added March 18, 2017

  INTEGER(C_INTPTR_T), SAVE :: local_nlastFixres, local_last_startFixres, local_nlast_after_transFixres
  INTEGER(C_INTPTR_T), SAVE :: local_last_start_after_transFixres, total_local_sizeFixres, local_startFixres, local_nFixres

  INTEGER(C_INTPTR_T), SAVE :: local_nlastHighres, local_last_startHighres, local_nlast_after_transHighres
  INTEGER(C_INTPTR_T), SAVE :: local_last_start_after_transHighres, total_local_sizeHighres, local_startHighres, local_nHighres


  !========================================================================== 
  !                            FFTW VARIABLES
  !==========================================================================
  complex(C_DOUBLE_COMPLEX), pointer :: Wxvec(:,:,:), Wyvec(:,:,:), Wzvec(:,:,:)
  INTEGER, PARAMETER :: FFTW_REAL_TO_COMPLEX=-1, FFTW_COMPLEX_TO_REAL=1
  INTEGER, PARAMETER :: FFTW_OUT_OF_PLACE=0, FFTW_IN_PLACE=8
  INTEGER, PARAMETER :: FFTW_USE_WISDOM=16, FFTW_THREADSAFE=128
  INTEGER, PARAMETER :: FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0
  INTEGER, PARAMETER :: FTW_SCRAMBLED_INPUT=8192, FFTW_SCRAMBLED_OUTPUT=16384
  
END MODULE
