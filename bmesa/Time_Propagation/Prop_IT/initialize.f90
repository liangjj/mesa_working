!***********************************************************************
                           MODULE initialization
                           INTERFACE init
                    MODULE PROCEDURE init_d, init_z
                       END INTERFACE init
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck init_z.f
!**begin prologue     init_z
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           initialize, arnoldi
!**purpose            driver for initialization of arnoldi vectors
!**                   and matrices.
!**description
!**author             schneider, barry (nsf)
!**source
!**references
!**routines called
!**end prologue       init_z
  SUBROUTINE init_z(vector,h_vector)
  USE arnoldi_global_rt
  USE small_matrices
  USE h_on_vector
  USE finite_element_matrix_multiply
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d,end)         :: vector, h_vector
  INTEGER                                :: nout, i
  CHARACTER (LEN=4)                      :: itoc
  IF(end < begin) THEN
     WRITE(iout,1)
     CALL lnkerr('bad call to schmidt routine')
  END IF
  CALL cschmt(vector,thresh,n3d,begin,end,nout,.true.,log_prp(2))
  IF(nout == 0) then
!
!    Try generating a random vector and reorthogonalize
!
     DO i=begin,end
          vector(:,i) = 0.d0
          vector(begin,begin) = 1.d0
     END DO
     CALL cschmt(vector,thresh,n3d,begin,end,nout,.true.,  &
                 log_prp(2))     
  END IF
  IF(nout == 0) THEN
     CALL lnkerr('quit arnoldi. no more trial vectors '    &
                // 'possible')
  END IF

  size = begin + nout - 1

!     initialize the effect of the hamiltonian on these vectors.

!     title='h on initial vectors'
  CALL h_v_z(vector(:,begin:size),h_vector(:,begin:size),nout)
!  CALL finite_element_m_v_z(vector(:,begin:size),h_vector(:,begin:size),nout)
!     initialize the small hamiltonian matrix.

  CALL h_small(vector,h_vector,.false.)
1 FORMAT(/,1X,'number of new vectors less than number of old'  &
              ' vectors:quit')
END SUBROUTINE init_z
!deck init_d.f
!**begin prologue     init_d
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           initialize, arnoldi
!**purpose            driver for initialization of arnoldi vectors
!**                   and matrices.
!**description
!**author             schneider, barry (nsf)
!**source
!**references
!**routines called
!**end prologue       init_d
  SUBROUTINE init_d(vector,h_vector)
  USE arnoldi_global_it
  USE h_on_vector
  USE small_matrices
  USE finite_element_matrix_multiply
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d,end)             :: vector, h_vector
  INTEGER                                :: nout, i
  CHARACTER (LEN=4)                      :: itoc
!  write(iout,*) 'iter=',iter,'begin=',begin,'end=',end, &
!                'nout=',nout
  IF(end < begin) THEN
     WRITE(iout,1)
     CALL lnkerr('bad call to schmidt routine')
  END IF
  CALL dschmt(vector,thresh,n3d,begin,end,nout,.true.,log_prp(2))
  IF(nout == 0) then
!
!    Try generating a random vector and reorthogonalize
!
     DO i=begin,end
          vector(:,i) = 0.d0
          vector(begin,begin) = 1.d0
     END DO
     CALL dschmt(vector,thresh,n3d,begin,end,nout,.true.,  &
                 log_prp(2))     
  END IF
IF(nout == 0) THEN
     CALL lnkerr('quit arnoldi. no more trial vectors '    &
                // 'possible')
  END IF
  size = begin + nout - 1
!  write(iout,*) 'iter=',iter,'begin=',begin,'end=',end,   &
!                'nout=',nout
!     initialize the effect of the hamiltonian on these vectors.

!     title='h on initial vectors'
!  CALL h_v_d(vector(:,begin:size),h_vector(:,begin:size),nout)
  CALL finite_element_m_v(vector(:,begin:size),h_vector(:,begin:size),nout)
!     initialize the small hamiltonian matrix.
  CALL h_small(vector,h_vector,.false.)
1 FORMAT(/,1X,'number of new vectors less than number of old'  &
              ' vectors:quit')
END SUBROUTINE init_d
END MODULE initialization
