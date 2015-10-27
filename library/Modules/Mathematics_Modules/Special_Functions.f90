MODULE Special_Functions
!***begin prologue     Special_Functions
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                                                            
!***keywords           associated legendre functions                                                              
!***author             schneider, b. i.(nsf)                                                                                        
!***source                                                                                                                    
!***purpose            Compute P_lm(x) and Q_lm(x) for all x                                                                     
!***description        For P_lm and Q_lm with x=(-1,1) upward recursion is stable.
!***                   For P_lm for abs(x) > 1 upward recursion is also stable
!***                   For Q_lm for abs(x) > 1 upward recursion is unstable and we use
!***                   downward recursion starting at a high value of and then renormalizing 
!***                   using the analytically known first function.  This is called Millers
!***                   algorithm and is well know in computing Bessel functions.  
!***                                                                                                  
!                                                                                                                                   
!***references                                                                                                                      
!***routines called                                                                                                                 
!***end prologue       Special_Functions                                                                               
!
  USE Data_Module
  USE Legendre
  IMPLICIT NONE
  REAL(idp), DIMENSION(:),                          &
             ALLOCATABLE                  :: x
  REAL(idp), DIMENSION(:),                          &
             ALLOCATABLE                  :: y
  INTEGER                                 :: n_points = 1
  LOGICAL                                 :: normalize = .true.
  LOGICAL                                 :: Derivative = .false.
  LOGICAL                                 :: Print_Functions = .false.
  LOGICAL                                 :: Print_Wronskian = .false.
  LOGICAL                                 :: input_values = .true.
  LOGICAL                                 :: test_wron = .false.
  REAL(idp)                               :: norm
  REAL(idp)                               :: arg
  REAL(idp)                               :: scale_factor
  REAL(idp)                               :: log_factor
  REAL(idp)                               :: wron
  REAL(idp), DIMENSION(:), ALLOCATABLE    :: Factorial
  INTEGER                                 :: s_fac
  REAL(idp)                               :: smallest = tiny(1.d0) * 1.d+04
  REAL(idp)                               :: biggest  = huge(1.d0) * 1.d-04
  REAL(idp)                               :: eps = 1.d-10
  REAL(idp)                               :: upper=1.d0
  REAL(idp)                               :: lower=-1.d0
  REAL(idp)                               :: step
  CHARACTER (LEN=8), DIMENSION(:),                  &
                     ALLOCATABLE          :: row_label
  CHARACTER (LEN=8), DIMENSION(:),                  &
                     ALLOCATABLE          :: col_label
  CHARACTER(LEN=48)                       :: title='Test Calculation'
  CHARACTER(LEN=24)                       :: Control='compute_functions'
  CHARACTER(LEN=24)                       :: recur = 'Miller'
  CHARACTER(LEN=16)                       :: Directive = 'regular'
!
  TYPE Xi
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_Small
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_Large
  END TYPE Xi
  TYPE Eta
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_1
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_2
  END TYPE Eta

  TYPE Reg_L
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_x
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: DF_x
  END TYPE Reg_L

  TYPE Reg_M
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F_x
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: DF_x
  END TYPE Reg_M

  TYPE Reg_LM
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: DF
       REAL(idp), DIMENSION(:,:,:),                 &
                  ALLOCATABLE             :: F_x
       REAL(idp), DIMENSION(:,:,:),                 &
                  ALLOCATABLE             :: DF_x
       TYPE(Xi)                           :: Xi  
       TYPE(Eta)                          :: Eta  
  END TYPE Reg_LM

  TYPE Irreg_L
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
       REAL(idp), DIMENSION(:,:,:),                 &
                  ALLOCATABLE             :: F_x
       REAL(idp), DIMENSION(:,:,:),                 &
                  ALLOCATABLE             :: DF_x
  END TYPE Irreg_L

  TYPE Irreg_M
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
       REAL(idp), DIMENSION(:,:,:),                 &
                  ALLOCATABLE             :: F_x
       REAL(idp), DIMENSION(:,:,:),                 &
                  ALLOCATABLE             :: DF_x
  END TYPE Irreg_M

  TYPE Irreg_LM
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: DF
       REAL(idp), DIMENSION(:,:,:),                 &
                  ALLOCATABLE             :: F_x
       REAL(idp), DIMENSION(:,:,:),                 &
                  ALLOCATABLE             :: DF_x
       TYPE(Xi)                           :: Xi  
       TYPE(Eta)                          :: Eta  
  END TYPE Irreg_LM

  TYPE Up
       CHARACTER(LEN=24)                  :: Dir
  END TYPE UP

  TYPE Down_A
       CHARACTER(LEN=24)                  :: Dir
  END TYPE Down_A

  TYPE Down_B
       CHARACTER(LEN=24)                  :: Dir
  END TYPE Down_B

  TYPE Down
       TYPE(Down_A)                       :: A  
       TYPE(Down_B)                       :: B  
  END TYPE Down

  TYPE CF_Legendre
       CHARACTER(LEN=24)                  :: Dir
  END TYPE CF_Legendre
 
  TYPE Coefficients
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: a
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: b
  END TYPE Coefficients
!
  TYPE Legendre_Functions
       TYPE(Reg_L)                        :: R_L
       TYPE(Reg_M)                        :: R_M
       TYPE(Reg_LM)                       :: R_LM
       TYPE(Irreg_L)                      :: I_L
       TYPE(Irreg_M)                      :: I_M
       TYPE(Irreg_LM)                     :: I_LM
       TYPE(Up)                           :: U
       TYPE(Down)                         :: D
       TYPE(Coefficients)                 :: C_Q
  END TYPE Legendre_Functions
!
  TYPE(Legendre_Functions)                :: Leg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck N_Factorial 
!***begin prologue     N_Factorial    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            N_Factorial 
!***description        N_Factorial
!***                     
!***references         none
!                      
!***routines called
!***end prologue       N_Factorial         
      Subroutine N_Factorial 
      IMPLICIT NONE
      INTEGER                  :: i
      Factorial(0) = int_one
      Factorial(int_one) = int_one
      DO i = int_two, l_max + m_max
         Factorial(i) = i * Factorial( i - int_one )
      END DO
END SUBROUTINE N_Factorial 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Special_Functions
