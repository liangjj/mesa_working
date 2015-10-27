MODULE Lentz_Thompson
!***begin prologue     Lentz_Thompson
!***date written       021231   (yymmdd)
!***revision date               (yymmdd)                                                                                            
!***keywords           Lentz_Thompson                                                              
!***author             schneider, b. i.(nsf)                                                                                        
!***source                                                                                                                    
!***purpose            Compute continued fractions using Lentz-Thompson algoritm                                                                     
!***description        
!***                   
!***                   
!***                   
!***                   
!***                   
!***                                                                                                  
!                                                                                                                                   
!***references                                                                                                                      
!***routines called                                                                                                                 
!***end prologue       Lentz_Thompson                                                                               
!
!                          Needed Modules
!
  USE accuracy
  USE Data_Module
  USE input_output
  USE Matrix_Print
  USE Special_Functions
  IMPLICIT NONE
!
!
                           INTERFACE Continued_Fractions
             MODULE PROCEDURE Continued_Fraction_Legendre
                       END INTERFACE Continued_Fractions                                                                    
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                      
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!*deck Continued_Fraction_Legendre
!***begin prologue     Continued_Fraction_Legendre  
!***date written       100206   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           continued fractions
!***author             schneider, barry (NSF)
!***source             
!***purpose            Compute continued fraction
!***description        
!***                   
!***references         none
!***routines called
!***end prologue       Continued_Fraction_Legendre  
  Subroutine Continued_Fraction_Legendre(CFL,f,x,nu,mu) 
  IMPLICIT NONE
  TYPE(CF_Legendre)                 :: CFL
  REAL(idp)                         :: x 
  INTEGER                           :: nu 
  INTEGER                           :: mu 
  INTEGER                           :: n_0 
  REAL(idp)                         :: a 
  REAL(idp)                         :: b 
  REAL(idp)                         :: f        
  REAL(idp)                         :: C        
  REAL(idp)                         :: D        
  REAL(idp)                         :: Del        
  REAL(idp)                         :: test        
  INTEGER                           :: count         
  INTEGER                           :: iwrite         
!
  f = tiny
  C = f
  D = zero
  test = one
  a = one 
  n_0 = nu
  b = zero
!  Write(iout,1) f, C, D, test, a, b, n_0
  count = 0
!  iwrite = 0
  DO While (test > eps )              
     count = count + int_one
!     iwrite = iwrite + int_one
     b = ( ( n_0 + n_0 + one ) * x ) / ( n_0 + mu )
     D = b + a * D
     IF ( D == zero ) THEN
          D = tiny
     END IF
     C = b + a / C
     IF ( C == zero ) THEN
          C = tiny
     END IF
     D = 1.d0 / D
     Del = C * D
     f = f * Del
     test = abs ( Del - one )
     a = - ( n_0 - mu + one ) / ( n_0 + mu )
     n_0 = n_0 + 1
!     IF ( iwrite == 50 ) THEN
!          iwrite = 0
!          Write(iout,2) count, f, C, D, test, a, b, n_0
!     END IF
  END DO
  Write(iout,3) count, test, f 
1 Format(5x, 'Initial Values',/,10x,                                    &
         'f    = ',d20.12,1x,'C = ',d15.8,1x,'D = ',d15.8,/,10x,        &
         'test = ',d15.8,1x,'a  = ',d15.8,1x,'b = ',d15.8,1x,'n_0 = ',i10)
2 Format(5x, 'Loop Values',5x,'Count = ',i5/,10x,                       &
         'f    = ',d20.12,1x,'C = ',d15.8,1x,'D = ',d15.8,/,10x,        &
         'test = ',d15.8,1x,'a  = ',d15.8,1x,'b = ',d15.8,1x,'n_0 = ',i10)
3 Format(1x,'Iteration Count = ',i5,1x,'Convergence = ',d20.12,1x,     &
            'Final Value of Continued Fraction = ',d20.12)
END SUBROUTINE Continued_Fraction_Legendre
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE Lentz_Thompson
