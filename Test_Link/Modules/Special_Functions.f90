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
       CHARACTER(LEN=1)                   :: Dir
  END TYPE UP

  TYPE Down_A
       CHARACTER(LEN=1)                   :: Dir
  END TYPE Down_A

  TYPE Down_B
       CHARACTER(LEN=1)                   :: Dir
  END TYPE Down_B

  TYPE Down
       TYPE(Down_A)                       :: A  
       TYPE(Down_B)                       :: B  
  END TYPE Down

  TYPE CF_Legendre
       CHARACTER(LEN=1)                   :: Dir
  END TYPE CF_Legendre

  TYPE CF_PI
       CHARACTER(LEN=1)                   :: Dir
  END TYPE CF_PI

  TYPE CF_e
       CHARACTER(LEN=1)                   :: Dir
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
                           INTERFACE C_ab
             MODULE PROCEDURE Legendre_C_ab
                       END INTERFACE C_ab 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Legendre_C_ab 
!***begin prologue     Legendre_C_ab    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of P_LM(x) functions using various recursion
!***                   relations.  The can be found at 
!***                   http://en.wikipedia.org/wiki/Associated_Legendre_function
!***                   On caution.  When the argument is outside the cut [-1.1] the
!***                   recursion given in that article has to be modified.  This can
!***                   seen by consulting the NBS handbook on mathematical functions.
!***                   Unfortunately, that reference does not contain all of the recurrances
!***                   needed.  I note some of this below.
!***references         none
!                      P_LM(z)[Q_LM(z)] are the regular[irregular] associated legendre 
!***                   functions.  z are the real values of the argument.  Between - 1 and + 1 upward
!                      recursion can be used for both the regular and irregular function.
!                      For other values of z upward recursion is fine for the regular
!                      function but downward recusion must be used for the irregular function.
!***routines called
!***end prologue       Legendre_C_ab   
      Subroutine Legendre_C_ab(C_Q)   
      IMPLICIT NONE
      TYPE(Coefficients)               :: C_Q
      ALLOCATE(C_Q%a(0:l_max,0:m_max), C_Q%b(0:l_max,0:m_max))
      DO m = 0, m_max
         DO l = 0, l_max
            C_Q%a(l,m) = ( l + l + int_one ) / ( l + m )
            C_Q%b(l,m) = ( l - m + int_one ) / ( l + m )
         END DO
      END DO
!
END SUBROUTINE Legendre_C_ab
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                                  
END MODULE Special_Functions
