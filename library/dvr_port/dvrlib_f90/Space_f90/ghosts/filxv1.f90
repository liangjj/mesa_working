!deck filxv1.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2003-01-12  Time: 12:43:15
 
!***begin prologue     filxv1
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           potential
!***author             schneider, barry (nsf)
!***source
!***purpose            fill grid and potential matrix
!***
!***description        fill the physical grid and one body potential
!***                   matrix arrays in the dvr representation.
!***
!***references

!***routines called
!***end prologue       filxv1

SUBROUTINE filxv1(pham,px,pv,n,word)

INTEGER*8, INTENT(IN OUT)                :: pham
INTEGER*8, INTENT(IN OUT)                :: px
INTEGER*8, INTENT(IN OUT)                :: pv
INTEGER, INTENT(IN)                      :: n
INTEGER, INTENT(IN OUT)                  :: word(2)
IMPLICIT INTEGER (a-z)
#ifdef decpointer

#END IF decpointer
#ifdef sgipointer

#END IF sgipointer
REAL*8 ham, x, v

COMMON/io/inp, iout
pointer (pham,ham(1))
pointer (px,x(1))
pointer (pv,v(1))

!     locate the needed arrays

h=1
vphy=h+n*n
h0=vphy+n
srf=h0+n*n
q0=srf+2
q1=q0+1
need=wpadti(1+n)
CALL getmem(need,px,word(1),'x',0)
need=wpadti(1+n)
CALL getmem(need,pv,word(2),'v',0)
CALL copy(ham(q1),x,n)
CALL copy(ham(vphy),v,n)
RETURN
END SUBROUTINE filxv1
