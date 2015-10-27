program coulg
IMPLICIT NONE
INTEGER          :: n, i
REAL*8           :: eig_in(100), eig_diff(100), ex
READ(5,*) n
DO i=1,n
   READ(5,*) eig_in(i)
END DO
STOP
DO i=1,n
   ex=-1.d0/(2.d0*i*i)
   eig_diff(i) = abs( (eig_in(i) - ex)/ex)
END DO
DO i=1,n
   WRITE(6,*) i, eig_diff(i)
END DO
STOP
END PROGRAM COULG
