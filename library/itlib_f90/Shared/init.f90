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

  SUBROUTINE init(drctv,nin,begin,size,rtdone,code,orth)
  USE io
  USE dvr_prnt
  USE dvr_global
  IMPLICIT NONE
  INTEGER, INTENT(IN)                      :: nin
  INTEGER, INTENT(IN)                      :: begin
  INTEGER, INTENT(OUT)                     :: size
  INTEGER, INTENT(IN)                      :: rtdone
  CHARACTER (LEN=*), INTENT(IN OUT)        :: code
  LOGICAL, INTENT(IN)                      :: orth
  CHARACTER (LEN=80)                       :: title
  CHARACTER (LEN=4)                        :: itoc
  INTEGER                                  :: nout
  INTEGER                                  :: end
!
!     get the additional set of vectors by orthogonalizing the trials
!     to the set of roots which have been converged.  this will ensure
!     that the subspace is orthogonal to previously converged vectors.

  nout = nin
  IF(orth) THEN
     IF(rtdone /= 0) THEN
        CALL invec(resid,code,n_dvd,rtdone,log_dvd(1))
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
  CALL hvdvr(hx,hy,hz,diag,vec(1,begin),hvec(1,begin),  &
             n_dvd,nx,ny,nz,nout,dim,prnt(3))

!     initialize the small hamiltonian matrix.

  CALL hsmall(b,bwrk,vec,hvec,n_dvd,begin,size,maxvec,drctv,.false.)
END SUBROUTINE init

