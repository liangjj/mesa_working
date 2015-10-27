!=======================================================================
  SUBROUTINE reform(str1, str2) 
!=======================================================================
!   Reformat the configuration string from the free-format str1 to the 
!   fixed 8(1x,a3,1x,i4,1x) format
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: str1
    CHARACTER(LEN=80), INTENT(OUT) :: str2

    INTEGER :: i, is, js
!
    
    i = 0
    str2 = '  '
    is = 0
    DO
     js = index(str1(is+1:),'(')
     if (js .ne. 0) then
       if (js .gt. 8) then
         write(*,*) 'ERROR: too many shells for this code'
         STOP
       end if
       i = i+5
       str2(i-js+1:i) = str1(is+1:is+js)
       is = is + js
       js = index(str1(is+1:),')')
       i = i+5
       str2(i-js+1:i) = str1(is+1:is+js)
       is = is + js
     ELSE
       EXIT
     END IF
    END DO
  
  END SUBROUTINE reform
