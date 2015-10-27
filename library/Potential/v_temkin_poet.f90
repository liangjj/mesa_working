!deck v_temkin_poet.f
!***begin prologue     v_temkin_poet
!***date written       000619   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           
!***author             schneider, barry (nsf)
!***source             
!***purpose            potential matrix elements of temkin-poet potential.
!***                   Grid in each dimension assumed the same.
!***references         
!
!***routines called    
!***end prologue       v_temkin_poet
      subroutine v_temkin_poet(pt,v_tp,n,prnt)
      USE input_output
      USE grid_global
      implicit none
      integer                     :: n, i, j
      real*8, dimension(n)        :: pt
      real*8, dimension(n,n)      :: v_tp
      logical                     :: prnt
      REAL*8                      :: r_fix 
      v_tp = 0.d0
      write(iout,1) 
      DO i=1,n
!
!            Loop for r_2 >= r_1
!           
         DO j=1,i
            v_tp(j,i) = 1.d0/pt(j)
         END DO
!
!            Loop for r_2 < r_1
!           
         r_fix = 1.d0/pt(i)
         DO j=i+1,n
            v_tp(j,i) = r_fix
         END DO                  
      END DO
      v_tp = - v_tp 
      IF(prnt) THEN
         call prntfmn('temkin-poet potential',v_tp,n,n,n,n,iout,'e')
      END IF
      return
 1    format(/,1x,'Temkin-Poet Potential')
END SUBROUTINE v_temkin_poet       



