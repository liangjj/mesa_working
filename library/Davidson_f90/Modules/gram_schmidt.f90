!deck gram_schmidt
  SUBROUTINE gram_schmidt(v_a,v_b,thresh,size,n_a,n_b,n_out,schmidt,prnt)
!***begin prologue     gram_schmidt
!***date written       960801  (yymmdd)
!***revision date              (yymmdd)

!***keywords           gram-schmidt, orthogonalization
!***author             barry schneider(nsf)
!***source
!***purpose            gram-schmidt orthogonalization.
!***description        a set of non-orthonormal vectors, v_b, are orthogonalized
!                      to another set of vectors, v_a, using a gram-schmidt process
!                      that checks for linear dependencies. the set of vectors,
!                      v_a, are already assumed to be internally orthonormal.

!                          v_a(n,*) = input vectors
!                          v_b(n,*) = input as non-orthogonal vectors and output
!                                    as orthonormal set.
!                          thresh  = acceptance tolerance on overlap
!                          n       = dimension of vector space
!                          nstart  = beginning vector
!                          nfin    = ending vector
!                          nout    = number vectors outputted
!                          schmdt  = perform either one or two
!                                    orthonormalizations on set
!***references
!***routines called    saxpy(clams), sdot(clams), sscal(clams)

!***end prologue       gram_schmidt
  USE input_output
  IMPLICIT NONE
  INTEGER                                  :: size
  INTEGER                                  :: n_a
  INTEGER                                  :: n_b
  INTEGER                                  :: n_out
  INTEGER                                  :: n_times
  INTEGER                                  :: upper
  INTEGER                                  :: trips
  INTEGER                                  :: i
  INTEGER                                  :: j
  REAL*8, DIMENSION(size,n_a)              :: v_a
  REAL*8, DIMENSION(size,n_b)              :: v_b
  REAL*8                                   :: thresh
  LOGICAL                                  :: schmidt
  LOGICAL                                  :: prnt
  REAL*8                                   :: overlap
  REAL*8                                   :: ddot
  REAL*8                                   :: norm
!
  n_times=1
  IF(schmidt) THEN
     n_times=2
  END IF
  upper=n_b
  DO trips=1,n_times
     n_out=0
     DO i=1,upper
        DO j=1,n_a
           overlap = ddot(size,v_a(1,j),1,v_b(1,i),1)
           CALL saxpy(size,-overlap,v_a(1,j),1,v_b(1,i),1)
        END DO
        norm=SQRT(ddot(size,v_b(1,i),1,v_b(1,i),1))
        IF(norm > thresh) THEN
           CALL sscal(size,1.0D+00/norm,v_b(1,i),1)
           n_out=n_out+1
           v_b(:,i) = v_b(:,n_out)
        END IF
    END DO
    upper=n_out
  END DO
  IF(prnt) THEN
     WRITE(iout,1) n_a, n_b
     WRITE(iout,2) n_out
  END IF
1    FORMAT(/,1X,'schmidt orthogonalization of two sets of vectors', /,1X,  &
    'set a already assumed to be orthonormal',/,1X,  &
    'number of vectors in set a = ',i4,/,1X, 'number of vectors in set b = ',i4)
2    FORMAT(/,1X,'number of set b vectors after orthogonalization '  &
    'to set a vectors = ',i4)
END SUBROUTINE gram_schmidt
