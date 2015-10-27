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
!***                   downward recursion starting eiher at a high value of l and then 
!***                   renormalizing using the analytically known first function or using
!***                   the continued fraction approach advocated by Segura and Gil. 
!***                   This first method, known as Millers algorithm suffers from the
!***                   defect that one has to have a good estimate of the starting value of
!***                   l which is not easy to estimate.  The second approach uses a continued
!***                   fraction, the known wronskian of P_lm and Q_lm and upward recursion 
!***                   for P_lm to get the starting values of Q_lm.  Downward recursion in l
!***                   gets the Q_lm for m = (0,1) and then upward recursion in m gets the
!***                   others.  This is accurate and self contained.  Note that the continued
!***                   fraction does converge slowly when the argument is near 1 but the Miller
!***                   algorithm is even worse since one does not have a clue where to start to
!***                   get accurate values.
!***
!***references
!***routines called
!***end prologue       Associated_Legendre_Functions
!
!                          Needed Modules
!
  USE accuracy
  USE Data_Module
  USE input_output
  USE Matrix_Print
  USE Special_Functions
  USE Lentz_Thompson
  IMPLICIT NONE
                                                                                                  
                           INTERFACE Legendre
             MODULE PROCEDURE Legendre                              
                       END INTERFACE Legendre
                                                                                                  
                           INTERFACE Legendre_Recursion
             MODULE PROCEDURE Upward_Regular_Legendre_Recursion_L,                         &
                              Upward_Regular_Legendre_Recursion_LM,                        &
                              Upward_Irregular_Legendre_Recursion_LM,                      &
                              Downward_Irregular_Legendre_Recursion_LM_A,                  &  
                              Downward_Irregular_Legendre_Recursion_LM_B
                       END INTERFACE Legendre_Recursion
!
                           INTERFACE Initialize
             MODULE PROCEDURE Initialize_Regular_L,                                        &
                              Initialize_Regular_LM,                                       &
                              Initialize_Irregular_L,                                      &
                              Initialize_Irregular_LM                              
                       END INTERFACE Initialize
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!               
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Legendre_Functions 
!***begin prologue     Legendre_Functions  
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of P_LM(x) functions using various recursion
!***                   relations.  
!***
!***references         The needed relations can be found at the following two references;
!***                   1. http://en.wikipedia.org/wiki/Associated_Legendre_function
!***                   2. Abramowitz and Stegun, Handbook of Mathematical Functions
!***                           UNITED STATES DEPARTMENT OF COMMERCE
!***                                       With
!***                   Formulas, Graphs, and Mathematical Tables
!***                   Edited by Milton Abramowitz and Irene A. Stegun
!***                   National Bureau of Standards, Applied Mathematics Series - 55
!***                   Issued June 1964, Tenth Printing, December 1972, with corrections
!***                   For sale by the Superintendent of Documents, U.S. Government Printing Office 
!***                   Washington, D.C. 20402 - Price $11.35 domestic postpaid, or $10.50 GPO Bookstore 
!***
!***                   Some comments.  The recursion relationships when the argument is inside or 
!***                   ouside the cut [-1,1] are slightly different.  Reference 2. does not contain 
!***                   all of the recurrances needed.  This is noted in the subroutines.
!                      P_LM(z)[Q_LM(z)] are the regular[irregular] associated legendre 
!***                   functions.  z are the real values of the argument.  Between - 1 and + 1 upward
!                      recursion can be used for both the regular and irregular function.
!                      For other values of z upward recursion is fine for the regular
!                      function but downward recusion must be used for the irregular function.
!***routines called
!***end prologue       Legendre_Functions 
      Subroutine Legendre ( R_LM, I_LM ) 
      IMPLICIT NONE
      TYPE(Reg_LM), OPTIONAL                         :: R_LM
      TYPE(Irreg_LM), OPTIONAL                       :: I_LM
      TYPE(Up)                                       :: U
      TYPE(Down)                                     :: D
      TYPE(Down_A)                                   :: A
      TYPE(Down_B)                                   :: B
      INTEGER                                        ::  i
      INTEGER                                        ::  j
      INTEGER                                        ::  lrow=1
      INTEGER                                        ::  lcol=1
      CHARACTER(LEN=3)                               ::  itoc
      CHARACTER(LEN=16)                              ::  fptoc
!
!
!----------------------------------------------------------------------c
!
!
!----------------------------------------------------------------------c
!
!----------------------------------------------------------------------c
!
!         Compute factorials
!
  ALLOCATE( Factor(0:l_max+m_max) )
  Call Factorials 
!
!         Set print labels
!
  ALLOCATE(col_label(0:m_max))
  DO i=0,m_max
     col_label(i) = 'm = '//itoc(i)
  END DO
  ALLOCATE(row_label(0:l_max))
  DO i=0,l_max
     row_label(i) = 'l = '//itoc(i)
  END DO
!
!         Allocate arrays
!
  IF ( PRESENT(R_LM) ) THEN
       ALLOCATE(Leg%R_LM%F(0:l_max,0:m_max))
       IF (Derivative) THEN
           ALLOCATE(Leg%R_LM%DF(0:l_max,0:m_max))
       END IF
  END IF
  IF ( PRESENT(I_LM) ) THEN
       ALLOCATE(Leg%I_LM%F(0:l_max,0:m_max))
       IF (Derivative) THEN
           ALLOCATE(Leg%I_LM%DF(0:l_max,0:m_max))
       END IF
       IF (Leg%D%B%Dir == 'Continued_Fraction' ) THEN
           ALLOCATE(Leg%R_L%F(0:l_max))
       END IF
  END IF
!
!        Allocate some space for storage of often used variables
!
  ALLOCATE(y(1:20))
  DO i = int_one, n_points
!----------------------------------------------------------------------c
!
!    The definition and calculation of Legendre functions depends if the
!    argument is ( -1.,1.) or outside that range.  Using s_fac allows a uniform
!    treatment.
     arg = x(i)
     s_fac = int_one
     IF ( abs(arg) > one ) THEN
          s_fac = - int_one
     END IF
     y(1) = one - arg * arg
     y(2) =  s_fac * y(1)
     y(3) = sqrt ( y(2) )   
     y(4) = y(3)   
     y(5) = arg * arg
     y(6) = y(5) * arg
     y(7) = y(6) * arg
     y(8) = y(7) * arg 
     y(9) = y(8) * arg 
     y(10) = y(9) * arg 
     y(11) = y(10) * arg 
     y(12) = y(11) * arg 
     IF ( PRESENT(R_LM) ) THEN
!
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!              Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!
          Call Initialize ( Leg%R_LM ) 
!
!----------------------------------------------------------------------c
!
!              Get the other L values by upward recursion
!
          Call Legendre_Recursion ( Leg%R_LM )
!----------------------------------------------------------------------c
!
!----------------------------------------------------------------------c
          write(iout,1) arg
          title='Regular Associated Legendre Functions'
          write(iout,2) title
          Call Print_Matrix(Leg%R_LM%F, l_max + 1, m_max + 1, iout, frmt='e',   &
                            collab=col_label, rowlab=row_label )
          IF ( abs(arg) <= one) THEN
               IF (normalize) THEN
!
!                  Normalize
!
                   Call Renormalize( Leg%R_LM%F )
                   title='Normalized Regular Associated Legendre Functions'
                   write(iout,2) title
                   Call Print_Matrix(Leg%R_LM%F, l_max + 1, m_max + 1, iout, frmt='e',  &
                                     collab=col_label, rowlab=row_label )
               END IF
          END IF
     END IF
     IF ( PRESENT(I_LM) ) THEN
!
          log_factor = log ( abs ( ( arg + one )/( arg - one ) ) )
!----------------------------------------------------------------------c
!
!         Starting values for Q_LM upward recursion.
!----------------------------------------------------------------------c
!
          IF( abs(arg) < one ) THEN
              write(iout,*) 'upward recursion'
              Call Initialize( Leg%I_LM )
!
!----------------------------------------------------------------------c
!             Get other values by upward recursion in L starting with
!             the M=0,1 values and then to all L and M by upward
!             recursion on both variables.
!----------------------------------------------------------------------c
              Call Legendre_Recursion( Leg%I_LM, Leg%U )
          ELSE
!----------------------------------------------------------------------c
!           Recur downward for m = 0,1 using either Millers algorithm  c
!           or the continued fraction.
!----------------------------------------------------------------------c
!
               write(iout,*) 'downward recursion'
               IF ( Leg%D%A%Dir == 'Miller' ) THEN
                    Call Legendre_Recursion(  Leg%I_LM, Leg%D, Leg%D%A )
               ELSE IF (Leg%D%B%Dir == 'Continued_Fraction' ) THEN
                    Call Legendre_Recursion(  Leg%I_LM, Leg%D, Leg%D%B )
               END IF
!                 
!
          END IF
          title='Irregular Associated Legendre Functions'
          write(iout,1) arg
          write(iout,2) title
          Call Print_Matrix(Leg%I_LM%F, l_max + 1, m_max + 1, iout, frmt='e',   &
                            collab=col_label, rowlab=row_label )
     END IF 
  END DO
  DEALLOCATE( y )
  DEALLOCATE( Factor )
  DEALLOCATE(col_label)
  DEALLOCATE(row_label)
  IF ( PRESENT(R_LM) ) THEN
       DEALLOCATE(Leg%R_LM%F)
       IF (Derivative) THEN
           DEALLOCATE(Leg%R_LM%DF)
       END IF
  END IF
  IF ( PRESENT(I_LM) ) THEN
       DEALLOCATE(Leg%I_LM%F)
       IF (Derivative) THEN
           DEALLOCATE(Leg%I_LM%DF)
       END IF
  END IF
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
!***                   for upward L recursion when only one M value needed.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Regular_L         
      Subroutine Initialize_Regular_L ( R_L )
      IMPLICIT NONE
      TYPE ( Reg_L )                              :: R_L
      INTEGER                                     :: n_1
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
  n_1 = int_one  
  DO l = int_one, m
     Leg%R_L%F(m) = - s_fac * n_1 * y(3) * Leg%R_L%F(m)
     n_1 = n_1 + int_two
  END DO
!
!              Calculate the second P term.
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
!***                   for upward L recursion when multiple M values are needed.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Regular_LM         
      Subroutine Initialize_Regular_LM ( R_LM )
      IMPLICIT NONE
      TYPE ( Reg_LM )                             :: R_LM
      INTEGER                                     :: n_1
      INTEGER                                     :: n_2

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
  n_1 = int_one 
  DO m = int_zero, m_max              
!
!               Initialize first value.
!
     Leg%R_LM%F(m,m) = one
     n_2 = int_one 
     DO l = int_one, m
        Leg%R_LM%F(m,m) = - s_fac * n_2 * y(3) * Leg%R_LM%F(m,m)
        n_2 = n_2 + int_two
     END DO
!
!               Calculate second value.
!
     IF (l_max /= m) THEN
!
!          Now calculate:
!                 P_(M+1)M
!
         Leg%R_LM%F(m+1,m) = n_1 * arg * Leg%R_LM%F(m,m)
         n_1 = n_1 + int_two
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
!***                   for upward L recursion when only one M is needed.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Irregular_L         
      Subroutine Initialize_Irregular_L ( I_L )
      IMPLICIT NONE
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
       Leg%I_L%F(0) = half * log_factor
       Leg%I_L%F(1) = arg * Leg%I_L%F(0) - one
  END IF
  IF ( m == int_one) THEN
       I_0 = half * log_factor
       Leg%I_L%F(0) =  - one / y(3)
       Leg%I_L%F(1) =  - y(3) * ( I_0 + arg / y(1) )
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
!***                   for upward L recursion when multiple M values needed.
!***references         none
!                      
!***routines called
!***end prologue       Initialize_Irregular_LM         
      Subroutine Initialize_Irregular_LM ( I_LM )
      IMPLICIT NONE
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
  Leg%I_LM%F(0,0) = half * log_factor
  Leg%I_LM%F(1,0) = arg * Leg%I_LM%F(0,0) - one
  IF ( m_max > int_zero) THEN
       Leg%I_LM%F(0,1) =  - one / y(3)
       Leg%I_LM%F(1,1) =  - y(3) * ( Leg%I_LM%F(0,0) + arg / y(1) )
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
!***description        calculation of Q_LM(z) functions using forward recursion in L.
!***                   ( L + 1 - M) Q_(L+1)M = ( 2*L + 1) z Q_LM - ( L + M ) Q_(L-1)M
!***                   Recursion started with the explicit forms of Q_MM and Q_(M+1)M
!***                   The forward recursion is stable for both the regular and irregular
!***                   functions as long as abs(z) <= one and L is not huge.  It is also stable
!***                   for the regular functions for abx(z) > one.  It is NOT stable for the
!***                   irregular functions under these conditions and backward recursion is 
!***                   required.
!***references         none
!                      
!***routines called
!***end prologue       Upward_Regular_Legendre_Recursion_L         
      Subroutine Upward_Regular_Legendre_Recursion_L ( R_L )
      IMPLICIT NONE
      TYPE ( Reg_L )                              :: R_L
      INTEGER                                     ::  l
      INTEGER                                     :: n_1
      INTEGER                                     :: n_2
      INTEGER                                     :: n_3

!
!         The upward recursion.  Stable for all values of z
!
  n_1 = m + m + int_three 
  n_2 = m + m + int_one
  n_3 = int_two
  DO l = m + int_two, l_max
    Leg%R_L%F(l) = ( n_1 * arg * Leg%R_L%F(l - int_one)     &
                                         -                  &
                     n_2 * Leg%R_L%F(l - int_two) ) / n_3  
    n_1 = n_1 + int_two
    n_2 = n_2 + int_one
    n_3 = n_3 + int_one
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
!***description        calculation of Q_LM(z) functions using forward recursion in L.
!***                   ( L + 1 - M) Q_(L+1)M = ( 2*L + 1) z Q_LM - ( L + M ) Q_(L-1)M
!***                   Recursion started with the explicit forms of Q_MM and Q_(M+1)M
!***                   The forward recursion is stable for both the regular and irregular
!***                   functions as long as abs(z) <= one and L is not huge.  It is also stable
!***                   for the regular functions for abx(z) > one.  It is NOT stable for the
!***                   irregular functions under these conditions and backward recursion is 
!***                   required.
!***references         none
!                      
!***routines called
!***end prologue       Upward_Regular_Legendre_Recursion_LM          
      Subroutine Upward_Regular_Legendre_Recursion_LM ( R_LM )
      IMPLICIT NONE
      TYPE ( Reg_LM )                               :: R_LM
      INTEGER                                       ::  l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3
!
   DO m = int_zero, m_max
!
!         The upward recursion.  Stable for all values of z
!
      n_1 = m + m + int_three 
      n_2 = m + m + int_one
      n_3 = int_two
      DO l = m + int_two, l_max
         Leg%R_LM%F(l,m) = ( n_1 * arg * Leg%R_LM%F(l - int_one,m)    &
                                       -                              &
                             n_2 * Leg%R_LM%F(l - int_two,m) ) / n_3
         n_1 = n_1 + int_two
         n_2 = n_2 + int_one
         n_3 = n_3 + int_one
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
!***description        calculation of Q_LM(z) functions using forward recursion in L.
!***                   This is used when abs(arg) is < one.
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
      Subroutine Upward_Irregular_Legendre_Recursion_LM ( I_LM, U )
      IMPLICIT NONE
      TYPE ( Irreg_LM )                             :: I_LM
      TYPE ( Up )                                   :: U
      INTEGER                                       ::  l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3
      Leg%U%Dir='U'
!
!
!         The upward recursion.  Stable for all values of z <= one
!
!----------------------------------------------------------------------c
!             Step up to get all Q_l0 and Q_l1
!----------------------------------------------------------------------c

      n_1 = int_three 
      n_2 = int_one
      n_3 = int_two      
      DO l = int_two, l_max
         Leg%I_LM%F(l,0) = ( n_1 * arg * Leg%I_LM%F(l - int_one,0)          &
                                                       -                                    &
                             n_2 * Leg%I_LM%F(l - int_two,0) ) / n_3
         n_1 = n_1 + int_two
         n_2 = n_2 + int_one
         n_3 = n_3 + int_one
      END DO
      IF ( m_max > int_zero ) THEN
!
           n_1 = int_three 
           n_2 = int_two
           n_3 = int_one      
           DO l = int_two, l_max
              Leg%I_LM%F(l,1) = ( n_1 * arg * Leg%I_LM%F(l - int_one,1)     &
                                                      -                                     &
                                      l * Leg%I_LM%F(l - int_two,1) ) / ( l - int_one )
              n_1 = n_1 + int_two
              n_2 = n_2 + int_one
              n_3 = n_3 + int_one
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
               n_1 = - int_two 
               n_2 = l + int_one
               n_3 = l      
              DO m = int_two, m_max
                 Leg%I_LM%F(l,m) = ( - m - m + int_two ) * arg * Leg%I_LM%F(l,m-1) / y(3)    &
                                                          -                                 &
                            s_fac * ( l + m - int_one) * ( l - m + int_two) * Leg%I_LM%F(l,m-2)
                 n_1 = n_1 - int_two
                 n_2 = n_2 + int_one
                 n_3 = n_3 - int_one
              END DO
           END DO
      END IF
END SUBROUTINE Upward_Irregular_Legendre_Recursion_LM  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Downward_Irregular_Legendre_Recursion_LM_A 
!***begin prologue     Downward_Irregular_Legendre_Recursion_LM_A    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using backward recursion
!***                   This is used when abs(arg) is > one.
!***                   ( L + 1 - M) T_LM = ( 2*L + 3) z T_(L+1)M - ( L - M + 2 ) T_(L+2)M
!***                    Starting at a large value of L set the last value to zero and the
!***                    next to last to one.  The recur downward which is the stable direction.
!***                    The T's are proportional to the desired Q functions.  
!***                    The proportionality constant is determined by the known value of Q00.
!***                    This allows us to compute the Q's for m=0. The process is repeated 
!***                    for Q_01
!***references         none
!                      
!***routines called
!***end prologue       Downward_Irregular_Legendre_Recursion_LM_A         
      Subroutine Downward_Irregular_Legendre_Recursion_LM_A ( I_LM, D, A )
      IMPLICIT NONE
      TYPE ( Irreg_LM )                             :: I_LM
      TYPE ( Down )                                 :: D
      TYPE ( Down_A )                               :: A
      REAL(idp)                                     :: I_2 
      REAL(idp)                                     :: I_1 
      REAL(idp)                                     :: I_0 
      INTEGER                                       :: l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3
      Leg%U%Dir='D'
      top_l = 40 + m_max + l_max
!
!     This is the hocus pokus I was talking about in the introduction.  It was part of
!     the original Specfun code and I frankly hated it.  There was no reason given why
!     the top_l changed when arg was less than 1.1 nor was the code in Specfun even correct
!     when the argument was negative.  So, I played around with the highest value to get 
!     something that "worked" but it was not satisfactory.  The continued fraction method is
!     far superior in this regard.
!
      IF (abs(arg) < 1.1d0 ) THEN
!          top_l = ( 40 + m_max + l_max) * int( - one - 1.8d0 * log( abs(arg) - one ) )
          top_l = ( 40 + m_max + l_max) * int( - one - 10.d0 * log( abs(arg) - one ) )
      END IF
!      top_l=max(5000,top_l)
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
      scale_factor =  half * log_factor  / Leg%I_LM%F(0,0)
      Leg%I_LM%F(0:l_max,0) = scale_factor * Leg%I_LM%F(0:l_max,0)
      IF ( m_max > int_zero) THEN
!
!         Downward recursion for m = 1
!
          m = int_one
          I_2 = zero
          I_1 = one
!
!         Downward recursion 
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
!         Renormalize
!
          scale_factor = ( - one / y(3) ) / Leg%I_LM%F(0,1)
          Leg%I_LM%F(0:l_max,1) = scale_factor * Leg%I_LM%F(0:l_max,1)
      END IF
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
            Leg%I_LM%F(l,m) = n_1 * arg * Leg%I_LM%F(l,m - int_one) / y(3)       &
                                           -                                     &
                              s_fac * n_2 * n_3 * Leg%I_LM%F(l, m - int_two)
            n_1 = n_1 - int_two
            n_2 = n_2 + int_one
            n_3 = n_3 - int_one
         END DO
      END DO
1 Format(/,10x,'Argument = ',e15.8,1x,'Top L = ',i5)
END SUBROUTINE Downward_Irregular_Legendre_Recursion_LM_A  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Downward_Irregular_Legendre_Recursion_LM_B 
!***begin prologue     Downward_Irregular_Legendre_Recursion_LM_B    
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of Q_LM(z) functions using a novel idea of Segura and Gil;
!***                   This is used when abs(arg) is < one.
!***                   1. continued fraction for Q_L0/Q_(L-1)0 and Q_L1/Q_(L-1)1
!***                   2. upward recursion for P_L0 and P_L1
!***                   3. the wronskian  P_L0 * Q_(L-1)0 - P_(L-1)0 * Q_L0  = 1 / L 
!***                      wronskian  P_L1 * Q_(L-1)1 - P_(L-1)1 * Q_L1  =  - 1 / L 
!***                      to compute the two  highest values of Q_L0 and Q_L1.  
!***                   4. downward recursion for all Q_L0 and Q_L1.  
!***                   5. upward recursion in M for all QLM.
!***references         Gil and Segura CPC 1998
!                      
!***routines called
!***end prologue       Downward_Irregular_Legendre_Recursion_LM_B         
      Subroutine Downward_Irregular_Legendre_Recursion_LM_B ( I_LM, D, B ) 
      IMPLICIT NONE
      TYPE ( Irreg_LM )                             :: I_LM
      TYPE ( Reg_L )                                :: R_L
      TYPE ( Down )                                 :: D
      TYPE ( Down_B )                               :: B
      TYPE(CF_Legendre)                             :: CFL
      REAL(idp)                                     :: I_2 
      REAL(idp)                                     :: I_1 
      REAL(idp)                                     :: I_0 
      REAL(idp)                                     :: cf 
      INTEGER                                       :: l
      INTEGER                                       :: n_1
      INTEGER                                       :: n_2
      INTEGER                                       :: n_3

      m = int_zero
!
!                      Recur up for P_L0
!
      Call initialize(Leg%R_L)
      Call Legendre_Recursion ( Leg%R_L )
!
!                      Get continued fraction
!
      Call Continued_Fraction_Legendre(CFL,cf,arg,l_max,m)       
      Leg%I_LM%F(l_max-int_one,m) = 1.d0                            &
                                         /                          &
                      ( l_max * ( Leg%R_L%F(l_max) - Leg%R_L%F(l_max-int_one) * cf ) )      
      Leg%I_LM%F(l_max,m) = cf * Leg%I_LM%F(l_max-int_one,m)            
!
!     Downward recursion for m = 0
!
      n_1 = l_max + l_max - int_one
      n_2 = l_max 
      n_3 = l_max - int_one
      DO l = l_max, int_two, - int_one
         Leg%I_LM%F(l-int_two,m) = ( n_1 * arg * Leg%I_LM%F(l-int_one,m)         &
                                                      -                          &
                                            n_2 * Leg%I_LM%F(l,m) ) / n_3 
         n_1 = n_1 - int_two
         n_2 = n_2 - int_one
         n_3 = n_3 - int_one
      END DO
      IF ( m_max > int_zero) THEN
           m = int_one
!
!                      Recur up for P_L1
!
           Call initialize(Leg%R_L)
           Call Legendre_Recursion ( Leg%R_L )
!
!                      Get continued fraction
!
           Call Continued_Fraction_Legendre(CFL,cf,arg,l_max,m)       
           Leg%I_LM%F(l_max-int_one,m) = - l_max / ( Leg%R_L%F(l_max) - Leg%R_L%F(l_max-int_one) * cf )   
           Leg%I_LM%F(l_max,m) = cf * Leg%I_LM%F(l_max-int_one,m)            
!           write(iout,*) Leg%I_LM%F(l_max,m), Leg%I_LM%F(l_max-int_one,m)
!
!                      Downward recursion for m = 1
!
           n_1 = l_max + l_max - int_one
           n_2 = l_max - int_one
           n_3 = l_max  
           DO l = l_max, int_two, - int_one
              Leg%I_LM%F(l-int_two,m) = ( n_1 * arg * Leg%I_LM%F(l-int_one,m)         &
                                                    -                                 &
                                          n_2 * Leg%I_LM%F(l,m) ) / n_3 
              n_1 = n_1 - int_two
              n_2 = n_2 - int_one
              n_3 = n_3 - int_one
           END DO
      END IF

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
            Leg%I_LM%F(l,m) = n_1 * arg * Leg%I_LM%F(l,m - int_one) / y(3)       &
                                           -                                     &
                              s_fac * n_2 * n_3 * Leg%I_LM%F(l, m - int_two)
            n_1 = n_1 - int_two
            n_2 = n_2 + int_one
            n_3 = n_3 - int_one
         END DO
      END DO
1 Format(/,5x,'m = ',i3,1x,'Continued Fraction = ',d15.8)
2 Format(/,5x,'Starting Values of Legendre Functions')
3 Format(/,5x,5e15.8)
END SUBROUTINE Downward_Irregular_Legendre_Recursion_LM_B  
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
!***description        Normalized regular Legendre functions from their un-normalized values.
!***                   This is only used on the cut (-1,+1).
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
END MODULE Associated_Legendre_Functions
