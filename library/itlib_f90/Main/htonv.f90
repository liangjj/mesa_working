!deck htonv.f
 
! Code converted using TO_F90 by Alan Miller
! Date: 2004-03-13  Time: 16:03:19
 
!***begin prologue     htonv
!***date written       960723   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           matrix vector
!***author             schneider, barry (nsf)
!***source             tproplib
!***purpose            matrix vector multiply for time plus space
!***                   dvr hamiltonian.
!***description        the vector vecin is transformed to the vector
!***                   vecout by operating with the hamiltonian.

!                      n = n3*n2*n1*nt
!                      nc = number of channels
!                      ht(nt,nt) = matrix representation of d
!                                                           --
!                                                           dt
!                      h1(n1,n1) = matrix representation of T(1) + V0(1)
!                      h2(n2,n2) = matrix representation of T(2) + V0(2)
!                      h3(n3,n3) = matrix representation of T(3) + V0(3)
!                      vecin(n,nc,nvc)  = vecin(n3,n2,n1,nt,nvc)
!                      vecout(n,nc,nvc) = vecout(n3,n2,n1,nt,nvc)
!                      v123t(n,nc,nc)   = interaction potential, diagonal in coordinate space
!                                         possible, off-diagonal in channels

!                            Structure of Subroutine

!                      [ -H  - d/dt ] [ VinReal ]     [ VoutReal ]
!                      [            ] [         ]  =  [          ]
!                      [ d/dt   -H  ] [ VinImag ]     [ VoutImag ]

!***references
!***routines called
!***end prologue       htonv

SUBROUTINE htonv(h1,h2,h3,ht,v123t,vecin,vecout,n,nc, n1,n2,n3,nt,nvc,dim)

REAL*8, INTENT(IN OUT)                   :: h1(n1,n1)
REAL*8, INTENT(IN OUT)                   :: h2(n2,n2)
REAL*8, INTENT(IN OUT)                   :: h3(n3,n3)
REAL*8, INTENT(IN OUT)                   :: ht(nt,nt)
REAL*8, INTENT(IN OUT)                   :: v123t(n,nc,nc)
REAL*8, INTENT(IN OUT)                   :: vecin(n,nc,2,nvc)
REAL*8, INTENT(IN OUT)                   :: vecout(n,nc,2,nvc)
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN OUT)                  :: nc
INTEGER, INTENT(IN OUT)                  :: n1
INTEGER, INTENT(IN OUT)                  :: n2
INTEGER, INTENT(IN OUT)                  :: n3
INTEGER, INTENT(IN OUT)                  :: nt
INTEGER, INTENT(IN OUT)                  :: nvc
INTEGER, INTENT(IN)                      :: dim
IMPLICIT INTEGER (a-z)



COMMON/io/inp, iout

CALL rzero(vecout,n*nc*2*nvc)
CALL vonv(vecout,v123t,vecin,n,nc,2*nvc)

!     now we need special cases

IF(dim == 1) THEN
  CALL h1tonv(h1,ht,vecout,vecin,n1,nt,nc,nvc)
ELSE IF(dim == 2) THEN
  CALL h2tonv(h1,h2,ht,vecout,vecin,n1,n2,nt,nc,nvc)
ELSE IF(dim == 3) THEN
  CALL h3tonv(h1,h2,h3,ht,vecout,vecin,n1,n2,n3,nt,nc,nvc)
END IF
RETURN
END SUBROUTINE htonv



