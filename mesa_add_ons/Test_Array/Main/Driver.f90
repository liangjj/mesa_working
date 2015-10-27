!deck DRIVER
!**begin prologue     DRIVER
!**date written       060711   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           
!**author             schneider, barry (nsf)
!**source
!**purpose            DRIVER test code
!**description        
!**                   
!**references
!**routines called        
!                       ----                    --------
!
!**end prologue       DRIVER
  PROGRAM DRIVER
  USE input_output
  USE Arrays
  IMPLICIT NONE
  LOGICAL                                  :: dollar
  INTEGER                                  :: intkey
  INTEGER                                  :: i
  CHARACTER(LEN=80)                        :: cpass
  CHARACTER(LEN=3200)                      :: card
  REAL(idp), ALLOCATABLE, DIMENSION(:)     :: test
!
  input = inp
  output = iout
!
!  Open the input and output files
!
  OPEN(input,file='Input',status='old')
  OPEN(output,file='Output',status='unknown')
!
  ALLOCATE(test(1:20))
  DO i=1,20
     test(i) = i
  END DO
  IF ( dollar('$data',card,cpass,inp) )THEN
       n_l = intkey(card,'lower_bound',1,' ')
       n_u = intkey(card,'upper_bound',1,' ')
  END IF
  Call Set_Array
  Call Make_Array(a,test,n_l,n_u)
  Call Make_Array(a,test(n_l:n_u),n_l,n_u)
  CLOSE(input)
  CLOSE(output)
  stop
END PROGRAM DRIVER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
