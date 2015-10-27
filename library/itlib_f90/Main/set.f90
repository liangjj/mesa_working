!deck set.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:04:30
 
!***begin prologue     set
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           set locations
!***author             schneider, barry (nsf)
!***source
!***routines called
!***end prologue      set

SUBROUTINE set(ham,v,eigv,s,eig,n)

INTEGER, INTENT(OUT)                     :: ham
INTEGER, INTENT(OUT)                     :: v
INTEGER, INTENT(OUT)                     :: eigv
INTEGER, INTENT(OUT)                     :: s
INTEGER, INTENT(OUT)                     :: eig
INTEGER, INTENT(IN)                      :: n
IMPLICIT INTEGER (a-z)
COMMON/io/inp, iout

ham=1
v=ham+n*n
eigv=v+n
s=eigv+n*n
eig=s+2*n
RETURN
END SUBROUTINE set
