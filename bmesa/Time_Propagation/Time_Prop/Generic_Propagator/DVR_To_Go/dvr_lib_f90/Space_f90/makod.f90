!deck makod.f
!***begin prologue     makod
!***date written       951229   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords
!***author             schneider, barry (nsf)
!***source
!***purpose            matrix elements between a bridge function in
!***                   region i - 1 and any function in region i.
!***                   if function in region i is not a bridge function
!***                   then there are no surface terms.  otherwise we need
!***                   to add the bloch operator contribution.

!***                   note that normalization factors are passed from 0:nfun
!***                   in this routine.  this allows us to use norm(0), which
!***                   is the normalization of the bridge function in region
!***                   i - 1 to make the needed integral.
!***references
!***routines called
!***end prologue       makod
!
  SUBROUTINE makod(hmat,kemat,pmat,prmat,norm,nr,nfun,prn)
  USE input_output
  IMPLICIT NONE
  INTEGER                                :: nr, nfun
  INTEGER                                :: i
  REAL*8, DIMENSION(*)                   :: hmat, pmat
  REAL*8, DIMENSION(nr,*)                :: kemat, prmat
  REAL*8, DIMENSION(0:nfun)              :: norm
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
!
   DO  i=1,nfun
       hmat(i) = hmat(i) + norm(i) * norm(0) * kemat(i,1)
       pmat(i) = pmat(i) + norm(i) * norm(0) * prmat(i,1)
   END DO
   IF(prn) THEN
      title='matrix elements between region (i-1) and region i'
      CALL prntrm(title,hmat,nfun,1,nfun,1,iout)
   END IF
END SUBROUTINE makod



