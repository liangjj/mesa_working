MODULE Associated_Legendre_Functions
!***begin prologue     Associated_Legendre_Functions
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                                                            
!***keywords           associated legendre functions                                                              
!***author             schneider, b. i.(nsf)                                                                                        
!***source                                                                                                                    
!***purpose            upward recursion for P_lm(x)                                                                      
!***description                                                                            
!***                                                                                                  
!                                                                                                                                   
!***references                                                                                                                      
!***routines called                                                                                                                 
!***end prologue       Associated_Legendre_Functions                                                                               
  USE accuracy
  IMPLICIT NONE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Legendre 
!***begin prologue     Legendre 
!***date written       091101   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           legendre functions
!***author             schneider, barry (lanl)
!***source             
!***purpose            legendre functions
!***description        calculation of P_LM(x) functions
!***references         none
!                      P_LM(X) are the associated legendre functions L=m to L=l_max    
!                      x are the values of cos(theta)
!                      dfct are the factorials from 0 to maxfac
!                      ddfct are the double factorials from 0 to maxfac    
!***routines called
!***end prologue       Legendre
      Subroutine Legendre 
      IMPLICIT NONE
      REAL(idp), DIMENSION(:,m:l_max)         ::  P_LM
      REAL(idp), DIMENSION(:)                 ::  x
      REAL(idp), DIMENSION(0:max_fac)         ::  dfct
      REAL(idp), DIMENSION(0:max_fac)         ::  ddfct
      INTEGER                                 ::  i
      INTEGER                                 ::  n_1
      INTEGER                                 ::  n_2
      INTEGER                                 ::  n_3
!----------------------------------------------------------------------c
!           start recursion with P_m,m and P_m+1,m                     c
!                      and recur upward                                c
!----------------------------------------------------------------------c
      P_LM(:,m:l_max) = zero
      IF ( m == int_zero ) 
         P_LM(:,m)=one
      ELSE
         P_LM(:,m) = ddfct(m) * (one - x(:)*x(:))**.5d+00*m
      END IF
      IF (l_max /= m) THEN
          P_LM(:,m+1) = ( m + m + int_one) * x(:) * P_LM(:,m)
          IF (l_max /= ( m + int_one)) THEN
              n_1 = int_two
              n_2 = m + m + int_three
              n_3 = n_2 - int_two
              DO i = m+int_two,l_max
                 P_LM(:,i) = ( n_2 * x(:) * P_LM(:,i - int_one) - n_3 * P_LM(:,i - int_two) ) / n_1
                 n_1 = n_1 + int_one
                 n_2 = n_2 + int_two
                 n_3 = n_3 + int_one
              END DO
          END IF
      END IF
      END SUBROUTINE Legendre
!----------------------------------------------------------------------c
!----------------------------------------------------------------------c
!*deck fact
!***begin prologue     fact
!***date written       880721   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           fact, link 2702, factorials
!***author             schneider, barry (lanl)
!***source             m2702
!***purpose            factorials
!***description        calculation of factorials
!***references         none
!***routines called
!***end prologue       fact
      Subroutine fact
      IMPLICIT NONE
      REAL(idp), DIMENSION(0:max_fac)       :: dfct
      REAL(idp), DIMENSION(0:max_fac)       :: ddfct
      INTEGER                               :: i 
!----------------------------------------------------------------------c
!               calculate factorials                                   c
!----------------------------------------------------------------------c
      dfct(0:1) = one
      IF (max_fac > int_one) THEN
          DO  i = int_two,max_fac
             dfct(i) = i * dfct(i-1)
          END DO
      END IF
!----------------------------------------------------------------------c
!           calculate (2*m-1) double factorial                         c
!----------------------------------------------------------------------c
      ddfct(0:1) = one
      ddfct(2) = three
      IF (max_fac > two) THEN
         DO i = int_three,max_fac
            ddfct(i)=(i+i-1)*ddfct(i-1)
         END DO
      END IF
     END SUBROUTINE fact
END MODULE Associated_Legendre_Functions
