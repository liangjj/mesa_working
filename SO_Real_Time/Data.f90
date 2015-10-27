!
  MODULE Data
  USE accuracy
  USE input_output,     ONLY: inp, iout
  IMPLICIT NONE
  REAL(idp)                                     ::          pi    = 3.141592653589793238462643383276D0 
  REAL(idp)                                     ::          two_pi= 6.283185307179586476925286766552D0 
  REAL(idp)                                     ::          zero    =  0.D0
  REAL(idp)                                     ::          quarter = .25D0
  REAL(idp)                                     ::          half    = .50D0
  REAL(idp)                                     ::          third   = .33333333333333333333333333333D0
  REAL(idp)                                     ::          fourth  = .25000000000000000000000000000D0
  REAL(idp)                                     ::          fifth   = .20000000000000000000000000000D0
  REAL(idp)                                     ::          sixth   = .16666666666666666666666666666D0
  REAL(idp)                                     ::          seventh = .14285714285714285714285714285D0
  REAL(idp)                                     ::          eighth  = .12500000000000000000000000000D0
  REAL(idp)                                     ::          ninth   = .11111111111111111111111111111D0
  REAL(idp)                                     ::          tenth   = .10000000000000000000000000000D0
  REAL(idp)                                     ::          one     = 1.0D0
  REAL(idp)                                     ::          two     = 2.0D0
  REAL(idp)                                     ::          three   = 3.0D0
  REAL(idp)                                     ::          four    = 4.0D0
  REAL(idp)                                     ::          five    = 5.0D0
  REAL(idp)                                     ::          six     = 6.0D0
  REAL(idp)                                     ::          seven   = 7.0D0
  REAL(idp)                                     ::          eight   = 8.0D0
  REAL(idp)                                     ::          nine    = 9.0D0
  REAL(idp)                                     ::          ten     = 10.D0
  INTEGER                                       ::          int_zero   = 0
  INTEGER                                       ::          int_one    = 1
!                                                                                                                              
!                    hbar in joule-sec                                                                                         
!                                                                                                                              
  REAL(idp)                                    ::          hbar = 1.054571596D-34
!                                                                                                                              
!                    electron mass in kg                                                                                       
!                                                                                                                              
  REAL(idp)                                    ::          massau = 9.10938188D-31
                                                                                                                              
!                    bohr radius in meters                                                                                     
!                                                                                                                              
  REAL(idp)                                    ::          lenau = 5.291772083D-11
!                                                                                                                              
!                    time for an electron to make one bohr orbit in seconds                                                    
!                                                                                                                              
  REAL(idp)                                    ::          timau    = 2.418884326D-17
  REAL(idp)                                    ::          efieldau = 5.14220624D+11
  REAL(idp)                                    ::          electric_field_to_intensity = 3.509338D+16
  REAL(idp)                                    ::          peak_electric_field = .2849540283D-03
  REAL(idp)                                    ::          pmass    = 1.67262158D-27
  REAL(idp)                                    ::          massn2p  = 1.00137841887D0
  REAL(idp)                                    ::          au_in_ev = 27.211396132D0
  CHARACTER (LEN=4096)                         :: ops
  CHARACTER (LEN=80)                           :: cpass
  CHARACTER (LEN=1600)                         :: card
  COMPLEX(idp), PARAMETER                      :: eye = (0.0d0,1.0d0)
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: X     
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: D 
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: E 
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: W
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: V
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: V_t
  REAL(idp),    DIMENSION(:,:),  ALLOCATABLE   :: Z 
  REAL(idp)                                    :: AD
  REAL(idp)                                    :: ABSTOL
  COMPLEX(idp), DIMENSION(:),    ALLOCATABLE   :: HAM_U
  COMPLEX(idp), DIMENSION(:),    ALLOCATABLE   :: HAM_L
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: HAM_D
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: WORK 
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: IWORK 
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: IFAIL 
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: Vector 
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: work_d
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: pt
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: wt
  COMPLEX(idp), DIMENSION(:),    ALLOCATABLE   :: work_z
  REAL(idp),    DIMENSION(:),    ALLOCATABLE   :: rwork
  REAL(idp)                                    :: Step_Size
  REAL(idp)                                    :: E_con
  REAL(idp)                                    :: rcond
  REAL(idp)                                    :: berr
  REAL(idp)                                    :: ferr
  LOGICAL                                      :: Get_Eigenvectors    ! GET EIGENVECTORS (TRUE)  
  LOGICAL,      DIMENSION(10)                  :: Prnt =  [ .false., .false., .false., .false., .false.,   &
                                                            .false., .false., .false., .false., .false. ] 
  LOGICAL                                      :: E_flag
  INTEGER,      DIMENSION(:),    ALLOCATABLE   :: ipiv
  INTEGER,      DIMENSION(10),   PARAMETER     :: file_no = [ 8,9,10,20,30,40,50,60,70,80 ]
  CHARACTER(LEN=80),                                                                               &
                DIMENSION(2:11), PARAMETER     :: file_nm = [ 'TD_in',    'TD_out',   'CN_LEN_R', &
                                                              'CN_LEN_I', 'CN_VEL_R', 'CN_VEL_I', &
                                                              'SO_R',     'SO_I',     'PROB_out', &
                                                              'VEC_out' ]
  CHARACTER(LEN=128)                           :: td_filename
  INTEGER                                      :: Number_of_Eigenvectors 
  INTEGER                                      :: M_Size       ! SIZE OF MATRIX
  INTEGER                                      :: test_M_Size       ! SIZE OF MATRIX
  INTEGER                                      :: st
  INTEGER                                      :: IL
  INTEGER                                      :: IU
  INTEGER                                      :: VL 
  INTEGER                                      :: VU
  INTEGER                                      :: M_FOUND
  INTEGER                                      :: NSPLIT
  INTEGER                                      :: INFO
  INTEGER,      DIMENSION(:),    ALLOCATABLE   :: ISPLIT
  INTEGER,      DIMENSION(:),    ALLOCATABLE   :: IBLOCK
  INTEGER                                      :: M_Eig
  INTEGER                                      :: lwork
  INTEGER                                      :: NO_Time_Steps
  INTEGER                                      :: Quad_Size
  CHARACTER (LEN = 1)                          :: Y_N = 'N'
  CHARACTER (LEN = 80)                         :: title
  CHARACTER(LEN=80)                            :: Method
  CHARACTER(LEN=80)                            :: pulse
  CHARACTER(LEN=80)                            :: type_potential
  CHARACTER(LEN=80)                            :: scratch
  REAL(idp)                                    :: time
  REAL(idp)                                    :: half_time
  REAL(idp)                                    :: Energy
  REAL(idp)                                    :: E_0
  REAL(idp)                                    :: E_T 
  REAL(idp)                                    :: A
  REAL(idp)                                    :: A_old
  REAL(idp)                                    :: delta_t
  REAL(idp)                                    :: Omega
  REAL(idp)                                    :: phase
  REAL(idp)                                    :: Pulse_Time         ! = DELTA_T * SIZE/2
  REAL                                         :: t_i                ! INITIAL TIME (TIMER)
  REAL                                         :: t_f 
  REAL(idp)                                    :: left_end
  REAL(idp)                                    :: right_end
  REAL(idp)                                    :: t_l
  REAL(idp)                                    :: t_u
  REAL(idp)                                    :: alpha
  REAL(idp)                                    :: cutoff
  REAL(idp)                                    :: depth
  REAL(idp)                                    :: charge
  REAL(idp)                                    :: awell
  REAL(idp)                                    :: e_c
  REAL(idp)                                    :: mass
  REAL(idp)                                    :: factor
  REAL(idp)                                    :: shift
  REAL(idp), DIMENSION(2)                      :: amp
  REAL(idp), DIMENSION(2)                      :: expnt
  REAL(idp)                                    :: pre_factor
  REAL(idp)                                    :: dscale
  INTEGER                                      :: n_p
  INTEGER                                      :: n_scale
  INTEGER                                      :: nwell
!********************************************************************************
  END MODULE Data
!********************************************************************************
