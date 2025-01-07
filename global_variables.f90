MODULE global_variables
  use, intrinsic :: iso_c_binding   ! Newly Added Jan 17, and the following C_INTPTR_T, type(C_PTR) and so on
  IMPLICIT NONE

  INTEGER, PARAMETER :: pr = KIND (0.0d0)
  !INTEGER, PARAMETER :: MAX_ITER = 1000                !original 1000 ! Maximal iterations of maxdEdt
  INTEGER, PARAMETER :: MAX_ITER = 0                !original 1000 ! Maximal iterations of maxdEdt ! todo change
  INTEGER, PARAMETER :: MAX_ITER_CONSTR = 100
  INTEGER, PARAMETER :: KappaPoints = 16
  REAL(pr), PARAMETER :: OPTIM_TOL = 1.0e-8_pr
  REAL(pr), PARAMETER :: CONSTR_TOL = 1.0e-10_pr
  REAL(pr), PARAMETER :: MACH_EPSILON = 2.0e-16_pr 
  REAL(pr), PARAMETER :: TAU_MAX = 10.0_pr              ! Maximal step
  REAL(pr), PARAMETER :: J_MAX = 1.0e15_pr

  REAL(pr), parameter :: WEIGHT = 1.0_pr              ! Newly added on April 21, 2017, WEIGHT*R(u)+(1-WEIGHT)
                                                      ! When WEIGHT=1.0, the problem is maximizing dE/dt, i.e., R(u)

  
!  integer, parameter :: int_WEIGHT = 100              ! Newly added on April 21, 2017, WEIGHT*R(u)+(1-WEIGHT)

!  CHARACTER(len=*), parameter :: HomeDir = "/home/yund0050/project/yund0050/maxdEdtHeli_100_02_Graham/2_01_WEIGHT100_N0128"   ! Newly added on May 8, 2017, for setting directory
!  CHARACTER(len=*), parameter :: HomeDir = "/work/yund0050/maxdEdtHeli_100_11/4_01_WEIGHT100_N0256"   ! Newly added on May 8, 2017, for setting directory
!  CHARACTER(len=*), parameter :: HomeDir = "./"   ! Newly added on May 8, 2017, for setting directory
  CHARACTER(len=*), parameter :: HomeDir = "./output/"   ! Newly added on May 8, 2017, for setting directory

  LOGICAL :: kappaTest
  LOGICAL :: toDealias 
  LOGICAL :: save_diag_NS
  LOGICAL :: save_data_NS
  LOGICAL :: save_diag_Constr
  LOGICAL :: save_data_Constr
  LOGICAL :: save_diag_Optim
  LOGICAL :: save_data_Optim
  LOGICAL :: save_diag_lineMin
  LOGICAL :: save_data_lineMin
  LOGICAL :: parallel_data
  LOGICAL :: add_pert
 
  INTEGER(C_INTPTR_T), DIMENSION(3), SAVE :: C_n
  INTEGER, DIMENSION(3), SAVE :: n                 ! Newly Added March 18, 2017
  INTEGER, SAVE :: n_dim
  INTEGER, SAVE :: K0_index, E0_index, NU_index, iguess, ConsType
  REAL(pr), SAVE :: E0, K0, PI, visc, lambda2, lambda1, alpha0, dV, Kcut, Kmax, lebesgueQ
  real(pr), dimension(:), allocatable, save :: lebesgueQlist
  !--NOTE: Kcut = cut frequency used for dealiasing. 
  !-       Kmax = maximal frequency present in solution.

  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: K1, K2, K3
  REAL(pr), DIMENSION (:,:,:,:), ALLOCATABLE, SAVE :: Uvec, Wvec
  ! REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: Ux, Uy, Uz
  ! REAL(pr), DIMENSION(:,:,:), ALLOCATABLE :: Wx, Wy, Wz
  
  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: fwork, fwork_Fixres, fworkHighres


  real(pr), save :: viscCoefficient, pressureCoefficient
  
  !========================================================================== 
  !                            MPI VARIABLES
  !==========================================================================
  INTEGER, SAVE :: rank, Statinfo, np
  INTEGER(C_INTPTR_T), SAVE :: C_local_alloc, C_local_N, C_local_k_offset       ! Newly Added Jan 17; Modified on March 18, 2017
  INTEGER, SAVE :: local_alloc, local_N, local_k_offset                         ! Newly Added March 18, 2017


  ! INTEGER(C_INTPTR_T), SAVE :: local_nlast, local_last_start, local_nlast_after_trans 
  INTEGER(C_INTPTR_T), SAVE :: local_last_start_after_trans, C_total_local_size
  INTEGER, SAVE :: total_local_size                                             ! Newly Added March 18, 2017


  INTEGER(C_INTPTR_T), SAVE :: local_nlastFixres, local_last_startFixres, local_nlast_after_transFixres
  INTEGER(C_INTPTR_T), SAVE :: local_last_start_after_transFixres, total_local_sizeFixres, local_startFixres, local_nFixres

  INTEGER(C_INTPTR_T), SAVE :: local_nlastHighres, local_last_startHighres, local_nlast_after_transHighres
  INTEGER(C_INTPTR_T), SAVE :: local_last_start_after_transHighres, total_local_sizeHighres, local_startHighres, local_nHighres


  !========================================================================== 
  !                            FFTW VARIABLES
  !==========================================================================
!  type(C_PTR) :: fwdplan, bwdplan
!  type(C_PTR) :: cxdata, cydata, czdata, ddata
!  complex(C_DOUBLE_COMPLEX), pointer :: fxdata(:,:,:), fydata(:,:,:), fzdata(:,:,:)
  complex(C_DOUBLE_COMPLEX), pointer :: Wxvec(:,:,:), Wyvec(:,:,:), Wzvec(:,:,:)
!  complex(C_DOUBLE_COMPLEX), pointer :: ftmp(:,:,:), gtmp(:,:,:)





     ! INTEGER, PARAMETER :: FFTW_FORWARD=-1, FFTW_BACKWARD=1
  INTEGER, PARAMETER :: FFTW_REAL_TO_COMPLEX=-1, FFTW_COMPLEX_TO_REAL=1
     ! INTEGER, PARAMETER :: FFTW_ESTIMATE=0, FFTW_MEASURE=1
  INTEGER, PARAMETER :: FFTW_OUT_OF_PLACE=0, FFTW_IN_PLACE=8
  INTEGER, PARAMETER :: FFTW_USE_WISDOM=16, FFTW_THREADSAFE=128
  INTEGER, PARAMETER :: FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0
  INTEGER, PARAMETER :: FTW_SCRAMBLED_INPUT=8192, FFTW_SCRAMBLED_OUTPUT=16384
     ! INTEGER(8), SAVE  :: fwdplan, bwdplan, fwdplan_Fixres, bwdplan_Fixres, fwdplan_Highres, bwdplan_Highres

END MODULE
