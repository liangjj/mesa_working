!deck blmat.f
!***begin prologue     blmat
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            lobatto blmat matrix elements
!***
!***description        computes the blmat matrix elements
!
!                      M = - 1./( 2 * mass ) [ fa(xleft) * fb'(xleft)
!                          - fa(xright)  * fb'(xright) ]
!***references

!***routines called
!***end prologue       blmat

  SUBROUTINE blmat(blr,far,dfbr,ptr,coord,nr,region)
  USE dvr_global,   ONLY  : mass, parity, iout
  USE dvr_prnt
  IMPLICIT NONE
  INTEGER                                :: nr
  REAL*8, DIMENSION(nr,nr)               :: blr, far, dfbr
  REAL*8, DIMENSION(nr)                  :: ptr
  CHARACTER (LEN=*)                      :: coord
  INTEGER                                :: region
  REAL*8                                 :: scale
  CHARACTER (LEN=80)                     :: title
  CHARACTER (LEN=3)                      :: itoc
  INTEGER                                :: i
  blr=0.d0
  scale = - .5D0/mass
  IF(coord /= 'rho') THEN
     blr(nr,:) = blr(nr,:) - far(nr,nr)*dfbr(nr,:)
     blr(1,:) = blr(1,:) + far(1,1)*dfbr(1,:)
  ELSE
     IF(parity == 'even') THEN
        blr(nr,:) = blr(nr,:) - 2.d0*ptr(nr)*ptr(nr)*far(nr,nr)*dfbr(nr,:)
        blr(1,:) = blr(1,:) + 2.d0*ptr(1)*ptr(1)*far(1,1)*dfbr(1,:)
     ELSE IF(parity == 'odd'.OR.parity == 'none') THEN
        blr(nr,:) = blr(nr,:) - ptr(nr)*far(nr,nr)*dfbr(nr,:)
        blr(1,:) = blr(1,:) + ptr(1)*far(1,1)*dfbr(1,:)
     END IF
  END IF
  blr=scale*blr
  IF(prn(3)) THEN
     title='unnormalized Blroch matrix for sector = '//itoc(region)
     CALL prntrm(title,blr,nr,nr,nr,nr,iout)
  END IF
END SUBROUTINE blmat



