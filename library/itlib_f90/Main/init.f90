!deck init.f
!***begin prologue     init
!***date written       010829   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           initialize, davidson
!***purpose            driver for initialization of davidson vectors
!***                   and matrices.
!***description
!***author             schneider, barry (nsf)
!***source
!***
!***references

!***routines called
!***end prologue       init

  SUBROUTINE init(type)
  USE io
  USE dvd_prnt
  USE dvd_global
  IMPLICIT NONE
  CHARACTER(LEN=*)                         :: type
  CHARACTER (LEN=4)                        :: itoc
  INTEGER                                  :: end
!
!     get the additional set of vectors by orthogonalizing the trials
!     to the set of roots which have been converged.  this will ensure
!     that the subspace is orthogonal to previously converged vectors.

  nout = nin
  IF(orth) THEN
     IF(rtdone /= 0) THEN
        CALL invec(resid,type,rtdone)
        CALL abschm(resid,vec(1,begin),thresh,n_dvd,rtdone,nin,    &
                    nout,.true.,.false.)
     END IF
  END IF
  END = begin + nout - 1
  IF(nout /= 0) THEN
     CALL gschmt(vec,thresh,n_dvd,begin,END,nout,.true.,log_dvd(2))
  END IF
  IF(nout == 0) THEN
     CALL lnkerr('quit davidson. no more trial vectors '// 'possible')
  END IF
  size = begin + nout - 1

!     initialize the effect of the hamiltonian on these vectors.

  title='h on initial vectors'
  CALL hvdvr
!     initialize the small hamiltonian matrix.

  CALL hsmall(.false.)
END SUBROUTINE init

