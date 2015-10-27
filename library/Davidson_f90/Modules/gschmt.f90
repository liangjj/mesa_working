!deck gschmt.f
SUBROUTINE gschmt(v,thresh,n,nstart,nfin,nout,schmdt,prnt)
!***begin prologue     gschmt
!***date written       960801  (yymmdd)
!***revision date              (yymmdd)

!***keywords           gram-schmidt, orthogonalization
!***author             barry schneider(nsf)
!***source
!***purpose            gram-schmidt orthogonalization.
!***description        a set of non-orthonormal vectors are input and
!                      orthonormalized using a gram-schmidt process that
!                      checks for linear dependencies. the routine outputs
!                      the orthonormal functions in the same locations as the
!                      input functions.  the number of output functions is
!                      less than or equal to the number of input functions.

!                          v(n,*) = input vectors
!                          thresh = acceptance tolerance on overlap
!                          n      = dimension of vector space
!                          nstart = beginning vector
!                          nfin   = ending vector
!                          nout   = number vectors outputted
!                          schmdt = perform either one or two
!                                   orthonormalizations on set
!***references
!***routines called    saxpy(clams), sdot(clams), sscal(clams)
!***end prologue       gschmt
REAL*8, INTENT(IN)                       :: v(n,*)
REAL*8, INTENT(IN)                       :: thresh
INTEGER, INTENT(IN OUT)                  :: n
INTEGER, INTENT(IN)                      :: nstart
INTEGER, INTENT(IN)                      :: nfin
INTEGER, INTENT(OUT)                     :: nout
LOGICAL, INTENT(IN)                      :: schmdt
LOGICAL, INTENT(IN)                      :: prnt
IMPLICIT INTEGER (a-z)


REAL*8 sdot, norm, ovrlap

CHARACTER (LEN=80) :: title

COMMON /io/inp, iout

IF(prnt) THEN
  title='trial vectors before orthonormalization'
  CALL prntrm(title,v,n,nfin,n,nfin,iout)
END IF
ntrial=nstart-1
nout=0
loop10:  DO  i=nstart,nfin
  norm=SQRT(sdot(n,v(1,i),1,v(1,i),1))
  IF(norm < thresh) CYCLE loop10
  CALL sscal(n,1.d0/norm,v(1,i),1)
  IF(ntrial /= 0) THEN
    DO  j=1,ntrial
      ovrlap=sdot(n,v(1,j),1,v(1,i),1)
      CALL saxpy(n,-ovrlap,v(1,j),1,v(1,i),1)
      norm=SQRT(sdot(n,v(1,i),1,v(1,i),1))
      IF (norm < thresh) CYCLE loop10
    END DO
    CALL sscal(n,1.0D+00/norm,v(1,i),1)
  END IF
  ntrial=ntrial+1
  nout=nout+1
  CALL copy(v(1,i),v(1,ntrial),n)
END DO loop10
IF(schmdt) THEN
  DO  i=nstart,ntrial
    norm=SQRT(sdot(n,v(1,i),1,v(1,i),1))
    CALL sscal(n,1.d0/norm,v(1,i),1)
    DO  j=1,i-1
      ovrlap=sdot(n,v(1,j),1,v(1,i),1)
      CALL saxpy(n,-ovrlap,v(1,j),1,v(1,i),1)
      norm=SQRT(sdot(n,v(1,i),1,v(1,i),1))
    END DO
    CALL sscal(n,1.0D+00/norm,v(1,i),1)
  END DO
END IF
IF(prnt) THEN
  WRITE(iout,1) nstart, nfin
  WRITE(iout,2) nout
  title='trial vectors after orthonormalization'
  CALL prntrm(title,v,n,nout,n,nout,iout)
END IF
1    FORMAT(/,5X,'orthogonalizing vectors = ',i5,' to vectors = ',i5)
2    FORMAT(/,5X,'number of output vectors = ',i5)


RETURN
END SUBROUTINE gschmt
