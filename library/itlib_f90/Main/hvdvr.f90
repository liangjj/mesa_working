!deck hvdvr.f
!***begin prologue     hvdvr
!***date written       010828   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           matrix vector
!***author             schneider, barry (nsf)
!***source             tproplib
!***purpose            matrix-vector, multiply, dvr.
!***description        multiply a hamiltonian matrix, which is the sum
!***                   of up to 3 one-body matrices, onto a vector.
!**
!***                          definitions
!                      n_dvd = nphy(3)*nphy(2)*nphy(1)
!                      h(ni,ni) = matrix representation of T(i) + V0(i)
!                                 with zero diagonals
!                      diag     = full diagonal
!                      it is assumed that;
!                      H = h1 + h2 + h3 + diag
!***references
!***routines called
!***end prologue       hvdvr

  SUBROUTINE hvdvr
  USE io
  USE dvd_prnt
  USE dvd_global
  IMPLICIT NONE
  INTEGER                                :: i
  hvec(1:n_dvd,begin:size) = 0.d0
  IF(dim == 1) THEN
     call h1v(vec(1,begin),hvec(1,begin),nout)
  ELSE IF(dim == 2) THEN
     CALL h2v(vec(1,begin),hvec(1,begin),nout)
  ELSE IF(dim == 3) THEN
     CALL h3v(vec(1,begin),hvec(1,begin),nout)
  END IF
!     take care of the diagonal
  DO  i=1,n_dvd
      hvec (i,begin:size) = hvec(i,begin:size) + diag(i) * hvec(i,begin:size)
  END DO
  IF(log_dvd(3)) THEN
     title='hamiltonian on vectors'
     CALL prntrm(title,hvec(1,begin),n_dvd,nout,n_dvd,nout,iout)
  END IF
END SUBROUTINE hvdvr



