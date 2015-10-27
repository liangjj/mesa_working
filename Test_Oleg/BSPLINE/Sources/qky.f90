   
!======================================================================
      Real(8) FUNCTION QKy (i1,j1,i2,j2,k)
!======================================================================
!
!                 k                          
!     Evaluates  Q (i1, j1; i2, j2) = 
!
!              k                     k
!             W (i1, j1; i2, j2) -  W (j1, i1; j2, i2) 
!
!----------------------------------------------------------------------


    IMPLICIT NONE
    
    Integer(4), Intent(in) :: i1,j1,i2,j2,k
    
    Real(8), External :: WKy
       
    Qky =  WKy(i1,j1,i2,j2,k) - WKy(j1,i1,j2,i2,k)
    
    END FUNCTION QKy

      
      
      
      