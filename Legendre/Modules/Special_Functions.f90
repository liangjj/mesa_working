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
  USE accuracy
  USE Data_Module
  IMPLICIT NONE
  REAL(idp), DIMENSION(:),                          &
             ALLOCATABLE                  :: x
  REAL(idp), DIMENSION(:),                          &
             ALLOCATABLE                  :: y
  INTEGER                                 :: m_max
  INTEGER                                 :: l_max
  INTEGER                                 :: n_points
  LOGICAL                                 :: normalize
  LOGICAL                                 :: Derivative
  REAL(idp)                               :: norm
  REAL(idp)                               :: arg
  REAL(idp)                               :: scale_factor
  REAL(idp)                               :: log_factor
  REAL(idp), DIMENSION(:), ALLOCATABLE    :: Factor
  INTEGER                                 :: l
  INTEGER                                 :: m
  INTEGER                                 :: m_sign
  INTEGER                                 :: s_fac
  INTEGER                                 :: top_l
  REAL(idp)                               :: tiny  = 2.25d-307
  REAL(idp)                               :: large = 1.797d+308
  REAL(idp)                               :: eps=epsilon(1.d0)
  CHARACTER (LEN=8), DIMENSION(:),                  &
                     ALLOCATABLE          :: row_label
  CHARACTER (LEN=8), DIMENSION(:),                  &
                     ALLOCATABLE          :: col_label
  CHARACTER(LEN=80)                       :: title
  CHARACTER(LEN=32)                       :: Method = 'recursion'
  CHARACTER(LEN=16)                       :: Function_Type
!
  TYPE Reg_L
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
  END TYPE Reg_L

  TYPE Reg_M
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
  END TYPE Reg_M

  TYPE Reg_LM
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: DF
  END TYPE Reg_LM

  TYPE Irreg_L
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
  END TYPE Irreg_L

  TYPE Irreg_M
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:),                     &
                  ALLOCATABLE             :: DF
  END TYPE Irreg_M

  TYPE Irreg_LM
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: F
       REAL(idp), DIMENSION(:,:),                   &
                  ALLOCATABLE             :: DF
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

  TYPE CF_PI
       CHARACTER(LEN=24)                  :: Dir
  END TYPE CF_PI

  TYPE CF_e
       CHARACTER(LEN=24)                  :: Dir
  END TYPE CF_e

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
!
  TYPE Bessel_Functions
       TYPE(Coefficients)                 :: C_J
  END TYPE Bessel_Functions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Special_Functions
