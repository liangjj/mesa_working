MODULE Regular_Associated_Legendre_Functions
!***begin prologue     Regular_Associated_Legendre_Functions
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                  
!***keywords           associated legendre functions                                         
!***author             schneider, b. i.(nsf)                                                 
!***source                                                                                    
!***purpose            Compute P_lm(x) all x                                 
!***description        See subroutine Info in driver codes for description.
!***references
!***routines called
!***end prologue       Regular_Associated_Legendre_Functions
!
!                          Needed Modules
!
  USE accuracy
  USE Data_Module
  USE input_output
  USE Matrix_Print
  USE Special_Functions
  IMPLICIT NONE
                                                                                                  
                           INTERFACE Legendre
             MODULE PROCEDURE Legendre                              
                       END INTERFACE Legendre
!                                                                                                  
                           INTERFACE Legendre_Recursion
             MODULE PROCEDURE Upward_Regular_Legendre_Recursion_L,                         &
                              Upward_Regular_Legendre_Recursion_LM
                       END INTERFACE Legendre_Recursion
!
                           INTERFACE Initialize
             MODULE PROCEDURE Initialize_Regular_L,                                        &
                              Initialize_Regular_LM                                        
                       END INTERFACE Initialize
!
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
!***description        calculation of P_LM(x) functions using recursion
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
      Subroutine Legendre ( R_LM ) 
      IMPLICIT NONE
      TYPE(Reg_LM)                                   :: R_LM
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

!         Set print labels
!
  ALLOCATE(col_label(int_zero:m_max))
  DO i=int_zero,m_max
     col_label(i) = 'm = '//itoc(i)
  END DO
  ALLOCATE(row_label(int_zero:l_max))
  DO i=int_zero,l_max
     row_label(i) = 'l = '//itoc(i)
  END DO
!
!        Allocate some space for storage of often used variables
!
  ALLOCATE(y(int_one:int_twenty))
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
     y(int_one) = one - arg * arg
     y(int_two) =  s_fac * y(int_one)
     y(int_three) = sqrt ( y(int_two) )   
     y(int_four) = y(int_three)   
     y(int_five) = arg * arg
     y(int_six) = y(int_five) * arg
     y(int_seven) = y(int_six) * arg
     y(int_eight) = y(int_seven) * arg 
     y(int_nine) = y(int_eight) * arg 
     y(int_ten) = y(int_nine) * arg 
     y(int_eleven) = y(int_ten) * arg 
     y(int_twelve) = y(int_eleven) * arg 
!
!
     Call Legendre_Recursion ( Leg%R_LM )
     IF (Print_Functions) THEN
         write(iout,1) arg
         title='Regular Associated Legendre Functions'
         write(iout,2) title
         Call Print_Matrix(Leg%R_LM%F, l_max + int_one, m_max + int_one, iout, frmt='e',          &
                           collab=col_label, rowlab=row_label )
     END IF
     IF (normalize) THEN
!
!        Normalize
!
         Call Renormalize( Leg%R_LM%F )
         IF (Print_Functions) THEN
             title='Normalized Regular Associated Legendre Functions'
             write(iout,2) title
             Call Print_Matrix(Leg%R_LM%F, l_max + int_one, m_max + int_one, iout, frmt='e',  &
                               collab=col_label, rowlab=row_label )
         END IF
     END IF
  END DO
  DEALLOCATE( y )
  DEALLOCATE(col_label)
  DEALLOCATE(row_label)
1 Format(/,25x,'Argument = ',f15.8)
2 Format(/,25x,a48)
3 Format(/,25x,'Cannot Compute Irregular Function for Argument One')
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
!              Overwrite starting value until the result is obtained.
!
  n_1 = int_one  
  DO l = int_one, m
     Leg%R_L%F(m) = - s_fac * n_1 * y(int_three) * Leg%R_L%F(m)
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
        Leg%R_LM%F(m,m) = - s_fac * n_2 * y(int_three) * Leg%R_LM%F(m,m)
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
         Leg%R_LM%F(m+int_one,m) = n_1 * arg * Leg%R_LM%F(m,m)
         n_1 = n_1 + int_two
     END IF
  END DO
END SUBROUTINE Initialize_Regular_LM  
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
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!              Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!
!
  Call Initialize ( Leg%R_L ) 
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              Get the other L values by upward recursion in l
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!
!         The upward recursion.  Stable for all values of z
!         This is just the standard l recursion in textbooks.
!
  n_1 = m + m + int_three 
  n_2 = m + m + int_one
  n_3 = int_two
  DO l = m + int_two, l_max
    Leg%R_L%F(l) = ( n_1 * arg * Leg%R_L%F(l - int_one)                &
                                         -                             &
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
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              
!              Starting values for P_LM.
!              
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!
!
  Call Initialize ( Leg%R_LM ) 
!
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!              Get the other L values by upward recursion
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!
!
   DO m = int_zero, m_max
!
!         The upward recursion.  Stable for all values of z
!         This is just the standardv reecursion in textbooks.
!
      n_1 = m + m + int_three 
      n_2 = m + m + int_one
      n_3 = int_two
      DO l = m + int_two, l_max
         Leg%R_LM%F(l,m) = ( n_1 * arg * Leg%R_LM%F(l - int_one,m)             &
                                       -                                       &
                             n_2 * Leg%R_LM%F(l - int_two,m) ) / n_3
         n_1 = n_1 + int_two
         n_2 = n_2 + int_one
         n_3 = n_3 + int_one
      END DO
   END DO
END SUBROUTINE Upward_Regular_Legendre_Recursion_LM  
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
      Factor(int_one) = int_one
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
      REAL(idp), DIMENSION(0:l_max,int_zero:m_max)            ::  F_lm
      INTEGER                                          :: l_fac
      DO m = int_zero, m_max
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
END MODULE Regular_Associated_Legendre_Functions
