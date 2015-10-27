MODULE Associated_Legendre_Functions
!***begin prologue     Associated_Legendre_Functions
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
!***end prologue       Associated_Legendre_Functions                                                                               
  USE accuracy
  USE Data_Module
  USE input_output
  USE Matrix_Print
  IMPLICIT NONE
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

  TYPE Down
       CHARACTER(LEN=1)                   :: Dir
  END TYPE Down

  TYPE Legendre_Functions
       TYPE(Reg_L)                        :: R_L
       TYPE(Reg_M)                        :: R_M
       TYPE(Reg_LM)                       :: R_LM
       TYPE(Irreg_L)                      :: I_L
       TYPE(Irreg_M)                      :: I_M
       TYPE(Irreg_LM)                     :: I_LM
       TYPE(Up)                           :: U
       TYPE(Down)                         :: D
  END TYPE Legendre_Functions
!
  REAL(idp), DIMENSION(:),                          &
             ALLOCATABLE                  :: x
  INTEGER                                 :: m_max
  INTEGER                                 :: l_max
  INTEGER                                 :: n_points
  LOGICAL                                 :: normalize
  LOGICAL                                 :: N_O
  LOGICAL                                 :: Derivative
  REAL(idp)                               :: norm
  REAL(idp)                               :: x_1
  REAL(idp)                               :: x_2
  REAL(idp)                               :: x_3
  REAL(idp)                               :: arg
  REAL(idp)                               :: scale_factor
  REAL(idp), DIMENSION(:), ALLOCATABLE    :: Factor
  INTEGER                                 :: l
  INTEGER                                 :: m
  INTEGER                                 :: m_sign
  INTEGER                                 :: s_fac
  INTEGER                                 :: top_l
  CHARACTER(LEN=80)                       :: title
                                                                                                  
                           INTERFACE Legendre
             MODULE PROCEDURE Legendre                              
                       END INTERFACE Legendre
                                                                                                  
                           INTERFACE Legendre_Recursion
             MODULE PROCEDURE Upward_Regular_Legendre_Recursion_L,                         &
                              Upward_Regular_Legendre_Recursion_LM,                        &
                              Upward_Irregular_Legendre_Recursion_LM,                      &
                              Downward_Irregular_Legendre_Recursion_LM                              
                       END INTERFACE Legendre_Recursion
!
                           INTERFACE Initialize
             MODULE PROCEDURE Initialize_Regular_L,                                        &
                              Initialize_Regular_LM,                                       &
                              Initialize_Irregular_L,                                      &
                              Initialize_Irregular_LM                              
                       END INTERFACE Initialize
!
                           INTERFACE mlpnm
             MODULE PROCEDURE mlpnm                              
                       END INTERFACE mlpnm
!
                           INTERFACE mlqnm
             MODULE PROCEDURE mlqnm                              
                       END INTERFACE mlqnm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Regular_Legendre_Functions 
!***begin prologue     Regular_Legendre_Functions  
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
!***end prologue       Regular_Legendre_Functions 
      Subroutine Legendre (  Leg, R_LM, I_LM ) 
      IMPLICIT NONE
      TYPE(Reg_LM), OPTIONAL                         :: R_LM
      TYPE(Irreg_LM), OPTIONAL                       :: I_LM
      TYPE(Legendre_Functions)                       :: Leg
      TYPE(Up)                                       :: U
      TYPE(Down)                                     :: D
      INTEGER                                        ::  i
      INTEGER                                        ::  j
      INTEGER                                        ::  lrow=1
      INTEGER                                        ::  lcol=1
      CHARACTER(LEN=3)                               ::  itoc
      CHARACTER(LEN=16)                              ::  fptoc
      CHARACTER (LEN=8), DIMENSION(:), ALLOCATABLE   :: rowlab
      CHARACTER (LEN=8), DIMENSION(:), ALLOCATABLE   :: collab
!
!
!----------------------------------------------------------------------c
!
!         This is just a phase factor for the functions (-1)**m
!         which is a matter of convention.  It is use when abs(x) < one
!
!----------------------------------------------------------------------c
!
!----------------------------------------------------------------------c
!
  ALLOCATE( Factor(0:l_max+m_max) )
  Call Factorials 
  ALLOCATE(collab(0:m_max))
  DO i=0,m_max
     collab(i) = 'm = '//itoc(i)
  END DO
  ALLOCATE(rowlab(0:l_max))
  DO i=0,l_max
     rowlab(i) = 'l = '//itoc(i)
  END DO
  IF ( PRESENT(R_LM) ) THEN
       ALLOCATE(Leg%R_LM%F(0:l_max,0:m_max))
       IF (Derivative == .TRUE.) THEN
           ALLOCATE(Leg%R_LM%DF(0:l_max,0:m_max))
       END IF
  END IF
  IF ( PRESENT(I_LM) ) THEN
       ALLOCATE(Leg%I_LM%F(0:l_max,0:m_max))
       IF (Derivative == .TRUE.) THEN
           ALLOCATE(Leg%I_LM%DF(0:l_max,0:m_max))
       END IF
  END IF
  DO i = int_one, n_points
!----------------------------------------------------------------------c
!
     arg = x(i)
     s_fac = int_one
     IF ( abs(arg) > one ) THEN
          s_fac = - int_one
     END IF
     x_1 = one - arg * arg
     x_2 =  s_fac * x_1
     x_3 = sqrt ( x_2 )   
     IF ( PRESENT(R_LM) ) THEN
!
          IF ( N_O == .TRUE. ) THEN
               Call mlpnm (Leg%R_LM%F,x=arg,n=l_max,m=m_max)
          ELSE
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!              Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!
               Call Initialize ( Leg, Leg%R_LM ) 
!
!----------------------------------------------------------------------c
!
!              Get the other L values by upward recursion
!
               Call Legendre_Recursion ( Leg, Leg%R_LM )
!----------------------------------------------------------------------c
!
!----------------------------------------------------------------------c
          END IF
          write(iout,1) arg
          title='Regular Associated Legendre Functions'
          write(iout,2) title
          Call Print_Matrix(Leg%R_LM%F, l_max + 1, m_max + 1, iout, frmt='e',collab=collab,          &
                            rowlab=rowlab )
          IF ( abs(arg) <= one) THEN
               IF (normalize) THEN
!
!                  Normalize
!
                   Call Renormalize( Leg%R_LM%F )
                   title='Normalized Regular Associated Legendre Functions'
                   write(iout,2) title
                   Call Print_Matrix(Leg%R_LM%F, l_max + 1, m_max + 1, iout, frmt='e',collab=collab, &
                                     rowlab=rowlab )
               END IF
          END IF
     END IF
     IF ( PRESENT(I_LM) ) THEN
!
          IF ( N_O == .TRUE. ) THEN
               Call mlqnm (Leg%I_LM%F,x=arg,n=l_max,m=m_max)
          ELSE
!----------------------------------------------------------------------c
!
!    Starting values for Q_LM upward recursion.
!----------------------------------------------------------------------c
!
               IF( abs(arg) < 1.0001d0 ) THEN
                   Call Initialize( Leg, Leg%I_LM )
!
!----------------------------------------------------------------------c
!             Get other values by upward recursion in L starting with
!             the M=0,1 values and then to all L and M by upward
!             recursion on both variables.
!----------------------------------------------------------------------c
                   Call Legendre_Recursion( Leg, Leg%I_LM, Leg%U )
               ELSE
!----------------------------------------------------------------------c
!           Recur downward for m = 0,1 using Millers algorithm         c
!           The algorithm starts with a large value of l, setting the  c
!           first two values to zero and one respectively.             c
!           By recurring downward one gets quantities proportional to  c
!           exact values to some degree of precision.  The precision   c
!           depends on how high you start the recursion.  That can be  c
!           estimated based on the maximum values of l and m.          c
!           A simple rescaling based on the first two values which are c
!           easy to find, gives what is required.  Then one simply     c
!           recurs upward in m for each fixed l to get what is desired c
!----------------------------------------------------------------------c
                   Call Legendre_Recursion( Leg, Leg%I_LM, Leg%D )
!                 
               END IF
!
          END IF
          title='Irregular Associated Legendre Functions'
          write(iout,1) arg
          write(iout,2) title
          Call Print_Matrix(Leg%I_LM%F, l_max + 1, m_max + 1, iout, frmt='e',collab=collab, &
                            rowlab=rowlab )
     END IF 
  END DO
1 Format(/,25x,'Argument = ',f15.8)
2 Format(/,25x,a48)
END SUBROUTINE Legendre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Initialize_Regular_L 
!***begin prologue     Initialize_Regular_L    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation starting values of regular Legendre functions
!***                   for upward L recursion.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Regular_L         
      Subroutine Initialize_Regular_L ( Leg, R_L )
      IMPLICIT NONE
      TYPE ( Legendre_Functions )                 :: Leg
      TYPE ( Reg_L )                              :: R_L
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!         Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!       Use: 
!              P_MM = - s_fac * ( 2*M - 1) *  sqrt ( s_fac * ( 1 - x* x ) ) * P_(M-1)(M-1)
!       To step up in M after initializing at one.  
!
  Leg%R_L%F(m) = one
!
!              Overwrite until the starting value is obtained
!
  DO l = int_one, m
     Leg%R_L%F(m) = - s_fac * ( l + l - int_one ) * x_3 * Leg%R_L%F(m)
  END DO
!
!              Calculate the second term.
!
  IF ( l_max > m ) THEN
!
!          Now calculate:
!                 P_(M+1)M
!
       Leg%R_L%F(m+1) = ( m + m + int_one) * arg * Leg%R_L%F(m)
  END IF
END SUBROUTINE Initialize_Regular_L  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Initialize_Regular_LM 
!***begin prologue     Initialize_Regular_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation starting values of regular Legendre functions
!***                   for upward L recursion.M
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Regular_LM         
      Subroutine Initialize_Regular_LM ( Leg, R_LM )
      IMPLICIT NONE
      TYPE ( Legendre_Functions )                 :: Leg
      TYPE ( Reg_LM )                             :: R_LM
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!         Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!
!
!----------------------------------------------------------------------c
!       Use: 
!              P_MM = - s_fac * ( 2*M - 1) *  sqrt ( s_fac * ( 1 - x* x ) ) * P_(M-1)(M-1)
!       To step up in M after initializing at one.  
!
  DO m = 0, m_max              
!
!               Initialize first value.
!
     Leg%R_LM%F(m,m) = one
     DO l = 1, m
        Leg%R_LM%F(m,m) = - s_fac * ( l + l - int_one ) * x_3 * Leg%R_LM%F(m,m)
     END DO
!
!               Calculate second value.
!
     IF (l_max /= m) THEN
!
!          Now calculate:
!                 P_(M+1)M
!
         Leg%R_LM%F(m+1,m) = ( m + m + int_one) * arg * Leg%R_LM%F(m,m)
     END IF
  END DO
END SUBROUTINE Initialize_Regular_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Initialize_Irregular_L 
!***begin prologue     Initialize_Irregular_L    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation starting values of regular Legendre functions
!***                   for upward L recursion.M
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Irregular_L         
      Subroutine Initialize_Irregular_L ( Leg, I_L )
      IMPLICIT NONE
      TYPE ( Legendre_Functions )                 :: Leg
      TYPE ( Irreg_L )                            :: I_L
      REAL(idp)                                   :: I_0
!
!
!         Starting values for Q_LM upward recursion.
!              Use:       
!                    Q_00 = .5 * ln ( abs( ( z + 1) /( z - 1))
!                    Q_10 = z * Q_00 - 1
!                    Q_01 = - 1 / sqrt ( s_fac * ( 1 - z * z ) )
!                    Q_11 = - sqrt ( s_fac * ( 1 - z * z ) ) * ( Q_00 + z / ( 1 - z * z ) )
!
!----------------------------------------------------------------------c
!
  IF ( m == int_zero) THEN
       Leg%I_L%F(0) = half * log ( abs ( ( arg + one )/( arg - one ) ) )                   
       Leg%I_L%F(1) = arg * Leg%I_L%F(0) - one
  END IF
  IF ( m == int_one) THEN
       I_0 = half * log ( abs ( ( arg + one )/( arg - one ) ) )                   
       Leg%I_L%F(0) =  - one / x_3
       Leg%I_L%F(1) =  - x_3 * ( I_0 + arg / x_1 )
  END IF
!
END SUBROUTINE Initialize_Irregular_L  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Initialize_Irregular_LM 
!***begin prologue     Initialize_Irregular_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation starting values of regular Legendre functions
!***                   for upward L recursion.M
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Irregular_LM         
      Subroutine Initialize_Irregular_LM ( Leg, I_LM )
      IMPLICIT NONE
      TYPE ( Legendre_Functions )                 :: Leg
      TYPE ( Irreg_LM )                           :: I_LM
!
!
!         Starting values for Q_LM upward recursion.
!              Use:       
!                    Q_00 = .5 * ln ( abs( ( z + 1) /( z - 1))
!                    Q_10 = z * Q_00 - 1
!                    Q_01 = - 1 / sqrt ( s_fac * ( 1 - z * z ) )
!                    Q_11 = - sqrt ( s_fac * ( 1 - z * z ) ) * ( Q_00 + z / ( 1 - z * z ) )
!
!----------------------------------------------------------------------c
!
  Leg%I_LM%F(0,0) = half * log ( abs ( ( arg + one )/( arg - one ) ) )                   
  Leg%I_LM%F(1,0) = arg * Leg%I_LM%F(0,0) - one
  IF ( m_max > int_zero) THEN
       Leg%I_LM%F(0,1) =  - one / x_3
       Leg%I_LM%F(1,1) =  - x_3 * ( Leg%I_LM%F(0,0) + arg / x_1 )
  END IF
!
END SUBROUTINE Initialize_Irregular_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Upward_Regular_Legendre_Recursion_L 
!***begin prologue     Upward_Regular_Legendre_Recursion_L   
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using recursion
!***                   ( L + 1 - M) Q_(L+1)M = ( 2*L + 1) z Q_LM - ( L + M ) Q_(L-1)M
!***                   Recursion started with the explicit forms of Q_MM and Q_(M+1)M
!***                   The forward recursion is stable for both the regular and irregular
!***                   functions as long as abs(z) <= one and L is not huge.  It is also stable
!***                   for the regular functions for abx(z) > one.  It is NOT stable for the
!***                   irregular functions under these conditions and backward recursion is required.
!***references         none
!                      
!***routines called
!***end prologue       Upward_Regular_Legendre_Recursion_L         
      Subroutine Upward_Regular_Legendre_Recursion_L ( Leg, R_L )
      IMPLICIT NONE
      TYPE ( Legendre_Functions )                 :: Leg
      TYPE ( Reg_L )                              :: R_L
      INTEGER                                     ::  l
!
!         The upward recursion.  Stable for all values of z
!
  DO l = m + int_two, l_max
    Leg%R_L%F(l) = ( ( l + l - int_one ) * arg * Leg%R_L%F(l - int_one)    &
                                         -                                 &
                    ( l + m - int_one ) * Leg%R_L%F(l - int_two) ) / ( l - m )
  END DO
END SUBROUTINE Upward_Regular_Legendre_Recursion_L 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Upward_Regular_Legendre_Recursion_LM 
!***begin prologue     Upward_Regular_Legendre_Recursion_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using recursion
!***                   ( L + 1 - M) Q_(L+1)M = ( 2*L + 1) z Q_LM - ( L + M ) Q_(L-1)M
!***                   Recursion started with the explicit forms of Q_MM and Q_(M+1)M
!***                   The forward recursion is stable for both the regular and irregular
!***                   functions as long as abs(z) <= one and L is not huge.  It is also stable
!***                   for the regular functions for abx(z) > one.  It is NOT stable for the
!***                   irregular functions under these conditions and backward recursion is required.
!***references         none
!                      
!***routines called
!***end prologue       Upward_Regular_Legendre_Recursion_LM          
      Subroutine Upward_Regular_Legendre_Recursion_LM ( Leg, R_LM )
      IMPLICIT NONE
      TYPE ( Legendre_Functions )                   :: Leg
      TYPE ( Reg_LM )                               :: R_LM
      INTEGER                                       ::  l
!
   DO m = int_zero, m_max
!
!         The upward recursion.  Stable for all values of z
!
      DO l = m + int_two, l_max
         Leg%R_LM%F(l,m) = ( ( l + l - int_one ) * arg * Leg%R_LM%F(l - int_one,m)    &
                                                       -                              &
                           ( l + m - int_one ) * Leg%R_LM%F(l - int_two,m) ) / ( l - m )
      END DO
   END DO
END SUBROUTINE Upward_Regular_Legendre_Recursion_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Upward_Irregular_Legendre_Recursion_LM 
!***begin prologue     Upward_Irregular_Legendre_Recursion_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using recursion
!***                   ( L + 1 - M) Q_(L+1)M = ( 2*L + 1) z Q_LM - ( L + M ) Q_(L-1)M
!***                   Recursion started with the explicit forms of Q_MM and Q_(M+1)M
!***                   The forward recursion is stable for both the regular and irregular
!***                   functions as long as abs(z) <= one and L is not huge.  It is also stable
!***                   for the regular functions for abx(z) > one.  It is NOT stable for the
!***                   irregular functions under these conditions and backward recursion is required.
!***references         none
!                      
!***routines called
!***end prologue       Upward_Irregular_Legendre_Recursion_LM          
      Subroutine Upward_Irregular_Legendre_Recursion_LM ( Leg, I_LM, U )
      IMPLICIT NONE
      TYPE ( Legendre_Functions )                   :: Leg
      TYPE ( Irreg_LM )                             :: I_LM
      TYPE ( Up )                                   :: U
      INTEGER                                       ::  l
      Leg%U%Dir='U'
!
!
!         The upward recursion.  Stable for all values of z <= one
!
!----------------------------------------------------------------------c
!             Step up to get all Q_l0 and Q_l1
!----------------------------------------------------------------------c
      
      DO l = int_two, l_max
         Leg%I_LM%F(l,0) = ( ( l + l - int_one ) * arg * Leg%I_LM%F(l - int_one,0)          &
                                                       -                                    &
                          ( l - int_one ) * Leg%I_LM%F(l - int_two,0) ) / l
      END DO
      IF ( m_max > int_zero ) THEN
!
           DO l = int_two, l_max
              Leg%I_LM%F(l,1) = ( ( l + l - int_one ) * arg * Leg%I_LM%F(l - int_one,1)     &
                                                      -                                     &
                                      l * Leg%I_LM%F(l - int_two,1) ) / ( l - int_one )
           END DO
      END IF
!
!----------------------------------------------------------------------c
!             Now for each L value, step up in M
!----------------------------------------------------------------------c
      IF ( m_max > int_one ) THEN
           DO l = int_zero, l_max
!
!                         The upward recursion in m
!
              DO m = int_two, m_max
                 Leg%I_LM%F(l,m) = ( - m - m + int_two ) * arg * Leg%I_LM%F(l,m-1) / x_3    &
                                                          -                                 &
                            s_fac * ( l + m - int_one) * ( l - m + int_two) * Leg%I_LM%F(l,m-2)
              END DO
           END DO
      END IF
END SUBROUTINE Upward_Irregular_Legendre_Recursion_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Downward_Irregular_Legendre_Recursion_LM 
!***begin prologue     Downward_Irregular_Legendre_Recursion_LM    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using backward recursion
!***                   ( L + 1 - M) T_LM = ( 2*L + 3) z T_(L+1)M - ( L - M + 2 ) T_(L+2)M
!***                    Starting at a large value of L set the last value to zero and the
!***                    next to last to one.  The recur downward which is the stable direction.
!***                    The T's are proportional to the desired Q functions.  The proportionality
!***                    constant is determined by the known value of Q00.  This allows us to compute the
!***                    Q's for m=0. 
!***references         none
!                      
!***routines called
!***end prologue       Downward_Irregular_Legendre_Recursion_LM         
      Subroutine Downward_Irregular_Legendre_Recursion_LM ( Leg, I_LM, D ) 
      IMPLICIT NONE
      TYPE ( Legendre_Functions )                   :: Leg
      TYPE ( Irreg_LM )                             :: I_LM
      TYPE ( Down )                                 :: D
      REAL(idp)                                     :: I_2 
      REAL(idp)                                     :: I_1 
      REAL(idp)                                     :: I_0 
      INTEGER                                       :: l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3
      Leg%U%Dir='D'
      top_l = 40 + m_max + l_max
      IF (abs(arg) < 1.1d0) THEN
          top_l = ( 40 + m_max + l_max) * int( - one - 1.8d0 * log( arg - one ) )
      END IF
      write(iout,1) arg, top_l
      m = int_zero
      I_2 = zero
      I_1 = one
!
!     Downward recursion for m = 0
!
      n_1 = top_l + top_l + int_three
      n_2 = top_l + int_two - m
      n_3 = top_l + int_one + m
      DO l = top_l, l_max + int_one, - int_one
         I_0 = ( n_1 * arg * I_1 - n_2 * I_2 ) / n_3
         I_2 = I_1
         I_1 = I_0
         n_1 = n_1 - int_two
         n_2 = n_2 - int_one
         n_3 = n_3 - int_one
      END DO
      DO l = l_max, int_zero, - int_one
         Leg%I_LM%F(l,m) = ( n_1 * arg * I_1 - n_2 * I_2 ) / n_3 
         I_2 = I_1
         I_1 = Leg%I_LM%F(l,m)
         n_1 = n_1 - int_two
         n_2 = n_2 - int_one
         n_3 = n_3 - int_one
      END DO
!
!     Renormalize
!
      scale_factor =  ( half * log ( abs ( ( arg + one )/( arg - one ) ) ) ) / Leg%I_LM%F(0,0)
      Leg%I_LM%F(0:l_max,0) = scale_factor * Leg%I_LM%F(0:l_max,0)
!
!     Downward recursion for m = 1
!
      m = int_one
      I_2 = zero
      I_1 = one
!
!     Downward recursion 
!
      n_1 = top_l + top_l + int_three
      n_2 = top_l + int_two - m
      n_3 = top_l + int_one + m
      DO l = top_l, l_max + int_one, - int_one
         I_0 = ( n_1 * arg * I_1 - n_2 * I_2 ) / n_3
         I_2 = I_1
         I_1 = I_0
         n_1 = n_1 - int_two
         n_2 = n_2 - int_one
         n_3 = n_3 - int_one
      END DO
      DO l = l_max, int_zero, - int_one
         Leg%I_LM%F(l,m) = ( n_1 * arg * I_1 - n_2 * I_2 ) / n_3 
         I_2 = I_1
         I_1 = Leg%I_LM%F(l,m)
         n_1 = n_1 - int_two
         n_2 = n_2 - int_one
         n_3 = n_3 - int_one
      END DO
!
!     Renormalize
!
      scale_factor = ( - one / x_3 ) / Leg%I_LM%F(0,1)
      Leg%I_LM%F(0:l_max,1) = scale_factor * Leg%I_LM%F(0:l_max,1)
!
!     For each l value, step up in m
!
      DO l = 0, l_max
         n_1 = - int_two
         n_2 = l + int_one
         n_3 = l
!
!        The upward recursion in m
!
         DO m = int_two, m_max
            Leg%I_LM%F(l,m) = n_1 * arg * Leg%I_LM%F(l,m - int_one) / x_3       &
                                          -                                     &
                              s_fac * n_2 * n_3 * Leg%I_LM%F(l, m - int_two)
            n_1 = n_1 - int_two
            n_2 = n_2 + int_one
            n_3 = n_3 - int_one
         END DO
      END DO
1 Format(/,10x,'Argument = ',e15.8,1x,'Top L = ',i5)
END SUBROUTINE Downward_Irregular_Legendre_Recursion_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Factorials 
!***begin prologue     Factorials    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            Factorials 
!***description        Factorials
!***                     
!***references         none
!                      
!***routines called
!***end prologue       Factorials         
      Subroutine Factorials 
      IMPLICIT NONE
      INTEGER                  :: i
      Factor(0) = int_one
      Factor(1) = int_one
      DO i = int_two, l_max + m_max
         Factor(i) = i * Factor( i - int_one )
      END DO
END SUBROUTINE Factorials 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Renormalize 
!***begin prologue     Renormalize   
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        Normalized regular Legendre functions from their un-normalized values
!***                     
!***references         none
!                      
!***routines called
!***end prologue       Renormalizen        
      Subroutine Renormalize ( F_lm )
      IMPLICIT NONE
      REAL(idp), DIMENSION(0:l_max,0:m_max)            ::  F_lm
      INTEGER                                          :: l_fac
      DO m = 0, m_max
         l_fac = m + m + int_one
         DO l = m , l_max
            norm = sqrt ( half * l_fac * Factor ( l - m ) / Factor ( l + m ) )
            l_fac = l_fac +int_two
            F_lm(l,m) = norm * F_lm(l,m)
         END DO
      END DO
END SUBROUTINE Renormalize 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine mlpnm(pnm,dpnm,x,n,m)
  USE input_output
  INTEGER                                   :: i
  INTEGER                                   :: j
  INTEGER                                   :: k
  INTEGER                                   :: n
  INTEGER                                   :: m
  INTEGER                                   :: ls
  REAL(idp)                                 :: x
  REAL(idp)                                 :: xq
  REAL(idp)                                 :: xs
  REAL(idp), dimension(0:n,0:m)             :: pnm
  REAL(idp), dimension(0:n,0:m),  optional  :: dpnm
!       ==========================================================
!       purpose: this program computes the associated legendre 
!                functions pnm(x) and their derivatives pnm'(x) 
!       input :  x --- argument of pnm(x)
!                n --- degree of pnm(x),  n = 0,1,2,...
!                m --- degree of pnm(x),  m = 0,1,2,...
!       output:  pnm(n,m) --- pnm(x)
!                dpnm(n,m) --- dpnm'(x)
!       example: x = 0.50
!          pnm(x):
!          m\n        1            2            3            4
!         --------------------------------------------------------
!           0      .500000     -.125000     -.437500     -.289063
!           1     -.866025    -1.299038     -.324760     1.353165
!           2      .000000     2.250000     5.625000     4.218750
!           3      .000000      .000000    -9.742786   -34.099750
!           4      .000000      .000000      .000000    59.062500
!
!          dpnm(x):
!          m\n        1            2            3            4
!         --------------------------------------------------------
!           0     1.000000     1.500000      .375000    -1.562500
!           1      .577350    -1.732051    -6.278684    -5.773503
!           2      .000000    -3.000000     3.750000    33.750000
!           3      .000000      .000000    19.485572      .000000
!           4      .000000      .000000      .000000  -157.500000
!       ==========================================================
!
  pnm(:,:)=0.0d0
  pnm(0,0) = 1.d0
  IF (abs(x) == 1.0d0) THEN
      DO i=1,n
         pnm(i,0)=x**i
      END DO
      IF ( PRESENT( dpnm) ) THEN
           dpnm(:,:)=0.0d0              
           DO i=1,n
              dpnm(i,0)=0.5d0*i*(i+1.0d0)*x**(i+1)
           END DO
           dpnm(:,1) = 1.d+300
           DO j = 1, n
              dpnm(j,2)=-0.25d0*(j+2)*(j+1)*j*(j-1)*x**(j+1)
           END DO
      END IF
  ELSE
      ls=1
      if (abs(x) > 1.0d0) THEN
          ls=-1
      END IF
      xq=sqrt(ls*(1.0d0-x*x))
      xs=ls*(1.0d0-x*x)
      DO i=1,m
         pnm(i,i)=-ls*(2.0d0*i-1.0d0)*xq*pnm(i-1,i-1)
      END DO
      DO  i=0,m
          pnm(i+1,i)=(2.0d0*i+1.0d0)*x*pnm(i,i)
      END DO
      DO  i=0,m
          DO  j=i+2,n
              pnm(j,i)=((2.0d0*j-1.0d0)*x*pnm(j-1,i)-    &
                      (i+j-1.0d0)*pnm(j-2,i))/(j-i)
          END DO
      END DO
      IF ( PRESENT(dpnm) ) THEN
           dpnm(:,:)=0.d0
           DO  j=1,n
               dpnm(j,0)=ls*j*(pnm(j-1,0)-x*pnm(j,0))/xs
           END DO
           DO  i=1,m
               DO j=i,n
                  dpnm(j,i)=ls*i*x*pnm(j,i)/xs+(j+i)     &
                                *(j-i+1.0d0)/xq*pnm(j,i-1)
               END DO
           END DO
      END IF
  END IF
END SUBROUTINE mlpnm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Subroutine mlqnm(qnm,dqnm,x,n,m)
  USE input_output
  INTEGER                                   :: i
  INTEGER                                   :: j
  INTEGER                                   :: k
  INTEGER                                   :: n
  INTEGER                                   :: m
  INTEGER                                   :: ls
  INTEGER                                   :: km
  REAL(idp)                                 :: x
  REAL(idp)                                 :: xq
  REAL(idp)                                 :: xs
  REAL(idp), dimension(0:n,0:m)             :: qnm
  REAL(idp), dimension(0:n,0:m),  optional  :: dqnm
  REAL(idp)                                 :: q0
  REAL(idp)                                 :: q1
  REAL(idp)                                 :: q10
  REAL(idp)                                 :: qf
  REAL(idp)                                 :: qf0
  REAL(idp)                                 :: qf1
  REAL(idp)                                 :: qf2
!
!       ===============================================================
!       purpose: this program computes the associated legendre  
!                functions qnm(x) and their derivatives qnm'(x) using
!       input :  x --- argument of qnm(x) 
!                n --- degree of qnm(x) ( n = 0,1,2,úúú )
!                m --- order of qnm(x)  ( m = 0,1,2,úúú )
!       output:  qnm(n,m) --- qnm(x)
!                dqnm(n,m) --- qnm'(x)
!       examples:
!
!       qnm(x):  x = 0.5
!       n\m      0           1           2           3           4
!       ---------------------------------------------------------------
!        0     .549306   -1.154701    1.333333   -5.388603   26.666667
!        1    -.725347   -1.053063    2.666667   -6.158403   32.000000
!        2    -.818663     .729806    4.069272  -12.316806   42.666667
!        3    -.198655    2.491853    -.493486  -23.778868   85.333333
!        4     .440175    1.934087  -11.036781   -9.325204  186.818394
!
!       qnm'(x): x = 0.5
!       n\m      0           1           2           3           4
!       ---------------------------------------------------------------
!        0    1.333333    -.769800    4.444444  -20.014809  145.777778
!        1    1.215973   -2.377159    3.555556  -24.633611  156.444444
!        2    -.842707   -5.185328    8.796526  -24.633611  199.111111
!        3   -2.877344   -1.091406   28.115454  -50.976710  227.555556
!        4   -2.233291   11.454786   25.483527 -197.068892  412.039838
!
!       qnm(x): x = 2.0
!       n\m      0           1           2           3           4
!       ---------------------------------------------------------------
!        0     .549306    -.577350    1.333333   -5.003702   26.666667
!        1     .098612    -.203274     .666667   -3.079201   18.666667
!        2     .021184    -.064946     .277089   -1.539601   10.666667
!        3     .004871    -.019817     .104220    -.679543    5.333333
!        4     .001161    -.005887     .036816    -.276005    2.427640
!
!       qnm'(x): x = 2.0
!       n\m      0           1           2           3           4
!       ---------------------------------------------------------------
!        0    -.333333     .384900   -1.111111    5.388603  -36.444444
!        1    -.117361     .249384    -.888889    4.618802  -32.000000
!        2    -.037496     .116680    -.519437    3.079201  -23.111111
!        3    -.011442     .046960    -.253375    1.720114  -14.222222
!        4    -.003399     .017331    -.110263     .849589   -7.748516
!       ===============================================================
!
  IF (abs(x) == 1.0d0) then
      qnm(:,:)=1.0d+300
      IF ( PRESENT (dqnm) ) THEN
           dqnm(:,:)=1.0d+300
      END IF
  ELSE
      ls=1
      IF (abs(x) > 1.0d0) THEN
          ls=-1
      END IF
      xs=ls*(1.0d0-x*x)
      xq=sqrt(xs)
      q0=0.5d0*log(dabs((x+1.0d0)/(x-1.0d0)))
      IF (abs(x) < 1.0001d0) then
          qnm(0,0)=q0
          qnm(1,0)=x*q0-1.0d0
          qnm(0,1)=-1.0d0/xq
          qnm(1,1)=-xq*(q0+x/(1.0d0-x*x))
          DO i=0,1
             DO j=2,n
                qnm(j,i)=((2.0d0*j-1.0d0)*x*qnm(j-1,i)           &
                        -(j+i-1.0d0)*qnm(j-2,i))/(j-i)
             END DO
          END DO
          DO j=0,n
             DO  i=2,m
                 qnm(j,i)=-2.0d0*(i-1.0d0)*x/xq*qnm(j,i-1)-ls*   &
                          (j+i-1.0d0)*(j-i+2.0d0)*qnm(j,i-2)
             END DO
          END DO
      ELSE
          IF (abs(x) > 1.1d0) THEN
              km=40+m+n
          ELSE
              km=(40+m+n)*int(-1.0-1.8*log(x-1.0))
          END IF
          qf2=0.0d0
          qf1=1.0d0
          DO k=km,0,-1
             qf0=((2*k+3.0d0)*x*qf1-(k+2.0d0)*qf2)/(k+1.0d0)
             IF (k <= n) THEN
                   qnm(k,0)=qf0
             END IF 
             qf2=qf1
             qf1=qf0
          END DO
          qnm(0:n,0)=q0*qnm(0:n,0)/qf0
          qf2=0.0d0
          qf1=1.0d0
          DO k=km,0,-1
             qf0=((2*k+3.0d0)*x*qf1-(k+1.0d0)*qf2)/(k+2.0d0)
             IF (k <= n) THEN
                 qnm(k,1)=qf0
             END IF
             qf2=qf1
             qf1=qf0
          END DO
          q10=-1.0d0/xq
          qnm(0:n,1)=q10*qnm(0:n,1)/qf0
          DO j=0,n
             q0=qnm(j,0)
             q1=qnm(j,1)
             DO i=0,m-2
                qf=-2.0d0*(i+1)*x/xq*q1+(j-i)*(j+i+1.0d0)*q0
                qnm(j,i+2)=qf
                q0=q1
                q1=qf
             END DO
          END DO
      ENDif
      IF ( PRESENT (dqnm) ) THEN
           dqnm(0,0)=ls/xs
           DO j=1,n
              dqnm(j,0)=ls*j*(qnm(j-1,0)-x*qnm(j,0))/xs
           END DO
           DO j=0,n
              DO i=1,m
                 dqnm(j,i)=ls*i*x/xs*qnm(j,i)+(i+j)*(j-i+1.0d0)/xq*qnm(j,i-1)
              END DO
          END DO
      END IF
   END IF
END SUBROUTINE mlqnm
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Associated_Legendre_Functions
