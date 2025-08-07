MODULE global_variables
  use, intrinsic :: iso_c_binding   ! Newly Added Jan 17, and the following C_INTPTR_T, type(C_PTR) and so on
  IMPLICIT NONE

  INTEGER, PARAMETER :: pr = KIND (0.0d0)

  !math params!
  real(pr), parameter :: visc = 1.0_pr
  
  !opt params!
  INTEGER, PARAMETER :: MAX_ITER = 99999                   !original 1000 ! Maximal iterations of maxdEdt
  integer, parameter :: banachGradIterMax = 999
  REAL(pr), PARAMETER :: OPTIM_TOL = 1.0e-7_pr            !original 1.0e-8_pr
  REAL(pr), PARAMETER :: banachIterTol = 1.0e-7_pr            !original 1.0e-8_pr
  REAL(pr), PARAMETER :: MACH_EPSILON = 2.0e-16_pr
  REAL(pr), PARAMETER :: TAU_MAX = 1.0e2_pr
  integer :: resol ! automatically when load or in maxdLpdt
  real(pr), save :: lambda1 = 0.1_pr
  logical :: useBanachGradient = .false.                   ! banach gradient (true) or hilbert gradient (false)
  logical :: useOrthogonalGradient = .true.
  logical :: useConjugateGradient = .true.
  logical :: useRiemannianGeometry = .true.
  integer :: resetMomentumTermEveryXiterations = 25                 ! <1 = never

  logical :: normalizeDirection = .true.
  logical :: use_e_u_auto_for_q_less_4 = .true.           ! automatically use e_u instead of u for q<4 and otherwise don't use e_u
  logical :: use_e_u_instead_of_uqMinus4 = .true.        ! calc |u|^{q-2} |e_u cdot partial_k u|^2 instead of |u|^{q-4} |u cdot partial_k u|^2 to avoid dividing by 0 
  logical :: dealiase_if_mult_by_e_u = .true.             ! dealiase even if just multiplied by e_u (seems to be wrong )

  !data params!

  CHARACTER(len=:), allocatable :: hostName
  CHARACTER(len=*), parameter :: HomeDir = "./output/"
  CHARACTER(len=*), parameter :: ncDir = HomeDir//"ncFiles/"
  CHARACTER(len=:), allocatable :: inputDir                 ! depends on the server -> see initialize
  CHARACTER(len=:), allocatable :: ConstraintDir
  CHARACTER(len=:), allocatable :: loadTempFunctionName


  !debug params!
  logical :: verboseOptimization = .false.                  ! verbosely output stuff to terminal
  logical :: tauDebugToConsole = .true.                     ! verbosely output stuff to terminal     
  real(pr) :: checkDivergenceTolerance = 1.0e-15            ! tolerance to still be considered divergence free
  real(pr) :: checkNormalTolerance = 1.0e-12                ! tolerance to still be considered orthogonal
  real(pr) :: checkAverageTolerance = 1.0e-3                ! tolerance such that velocity < checkAverageTolerance*constraintB 
                                                                        !is still be considered average free     

  LOGICAL :: kappaTest = .true.
  LOGICAL :: toDealias = .true.
  LOGICAL :: mnbra_calcSaveAllJvalues = .false.              ! calculates J(u+tau d) for "all" tau values to get an idea of the shape of J(tau)
  integer :: save_uvecEveryXiteration = 100                  ! <1 for never
  integer :: save_scalarFieldsEveryXiteration = -1           ! <1 for never
  integer :: save_spectraEveryXiteration = 100               ! <1 for never
  logical :: normalizeSpectrumByL2Norm = .true.              ! when calculating the spectrum calc ||u||_2^2/sum(spec)*spec instead of just spec   
  integer :: dividingByZeroWarnings = 100                    ! number of warnings when calculating |u|^{-...} where u=0 
  LOGICAL :: save_diag_fields_values = .false.
  LOGICAL :: save_diag_Constr = .true.
  LOGICAL :: save_data_Constr = .true.
  LOGICAL :: save_diag_Optim = .true.
  LOGICAL :: save_data_Optim = .true.
  LOGICAL :: save_diag_lineMin = .true.
  LOGICAL :: save_data_lineMin = .true.
  LOGICAL :: add_pert = .false.
  LOGICAL :: parallel_data


  !other params!
  INTEGER, DIMENSION(3), SAVE :: n
  INTEGER, SAVE :: K0_index, E0_index, NU_index, iguess, ConsType
  REAL(pr), SAVE :: E0, K0, PI, dV, Kmax, lebesgueQ
  character(len=:), allocatable :: lebesgueQTxt
  real(pr), dimension(:), allocatable, save :: lebesgueQlist
  real(pr), dimension(:,:), allocatable :: optimizationResultList
  real(pr), dimension(:), allocatable :: B_list
  integer :: B_list_iterator
  character(3) :: bIterTxt
  character(7) :: Btxt
  integer :: bIterOffset
  integer :: optimizationIterOffset
  character(5) :: optimizationIterationTxt = "nan"

  REAL(pr), DIMENSION (:), ALLOCATABLE, SAVE :: K1, K2, K3
  REAL(pr), DIMENSION (:,:,:,:), ALLOCATABLE, SAVE :: Uvec

  real(pr) :: viscCoefficient = 1.0_pr                                          ! for debugging to turn on and off the viscocity/pressure terms
  real(pr) :: pressureCoefficient = 1.0_pr                                      ! for debugging to turn on and off the viscocity/pressure terms
  real(pr) :: BanachGradientWCoefficient = 1.0_pr                               ! for debugging to turn on and off the derivative in W^{1,s} Banach gradient
  real(pr) :: BanachGradientLCoefficient = 1.0_pr                               ! for debugging to turn on and off the 0-th order in W^{1,s} Banach gradient
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
