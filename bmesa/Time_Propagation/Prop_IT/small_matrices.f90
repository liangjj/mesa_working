!***********************************************************************
                           MODULE small_matrices
                           INTERFACE h_small
                    MODULE PROCEDURE h_small_d, h_small_z
                       END INTERFACE h_small
!***********************************************************************
                           CONTAINS
!***********************************************************************
!***********************************************************************
!***********************************************************************
!deck h_small_z.f
!**begin prologue     h_small_z
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           small arnoldi matrix
!**author             schneider, barry (nsf)
!**source
!**references
!**routines called
!**end prologue       h_small_z
  SUBROUTINE h_small_z(v_in,v_out,prnt)
  USE arnoldi_global_rt
  IMPLICIT NONE
  COMPLEX*16, DIMENSION(n3d,end)         :: v_in, v_out
  LOGICAL                                :: prnt
  COMPLEX*16                             :: cdotc
  INTEGER                                :: i, j 
! The matrix elements of the "small" matrix are computed as,
!\begin{equation}
!  b_{i,j} = \sum_{k=1,n} \sum_{l=1,n} V^{*}_{k,i} H_{k,l} V_{l,j}
!\end{equation}

  IF(drct == 'initialize') THEN
     DO  i=1,end
         DO  j=1,i
             b(i,j) = cdotc(n3d,v_in(1,i),1,v_out(1,j),1)
             b(j,i) = CONJG(b(i,j))
             bwrk(i,j) = b(i,j)
             bwrk(j,i) = b(j,i)
         END DO
     END DO
  ELSE IF(drct == 'fill') THEN
     DO  i=1,end
         DO  j=begin,end
             b(i,j) = cdotc(n3d,v_in(1,i),1,v_out(1,j),1)
             b(j,i) = CONJG(b(i,j))
         END DO
     END DO
     DO  i=1,END
         DO  j=1,i
             bwrk(i,j) = b(i,j)
             bwrk(j,i) = b(j,i)
         END DO
     END DO
  END IF
  IF(prnt) THEN
     title='small matrix'
     CALL prntcm(title,bwrk,end,maxvec,end,maxvec,iout)
  END IF
END SUBROUTINE h_small_z
!deck h_small_d.f
!**begin prologue     h_small_d
!**date written       010829   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           small arnoldi matrix
!**author             schneider, barry (nsf)
!**source
!**references
!**routines called
!**end prologue       h_small_d
  SUBROUTINE h_small_d(v_in,v_out,prnt)
  USE arnoldi_global_it
  IMPLICIT NONE
  REAL*8, DIMENSION(n3d,end)         :: v_in, v_out
  LOGICAL                            :: prnt
  REAL*8                             :: ddot
  INTEGER                            :: i, j 
! The matrix elements of the "small" matrix are computed as,
!\begin{equation}
!  b_{i,j} = \sum_{k=1,n} \sum_{l=1,n} V^{*}_{k,i} H_{k,l} V_{l,j}
!\end{equation}
  IF(drct == 'initialize') THEN
     DO  i=1,end
         DO  j=1,i
             b(i,j) = ddot(n3d,v_in(1,i),1,v_out(1,j),1)
             b(j,i) = b(i,j)
             bwrk(i,j) = b(i,j)
             bwrk(j,i) = b(j,i)
         END DO
     END DO
  ELSE IF(drct == 'fill') THEN
     DO  i=1,end
         DO  j=begin,end
             b(i,j) = ddot(n3d,v_in(1,i),1,v_out(1,j),1)
             b(j,i) = b(i,j)
         END DO
     END DO
     DO  i=1,END
         DO  j=1,i
             bwrk(i,j) = b(i,j)
             bwrk(j,i) = b(j,i)
         END DO
     END DO
  END IF
  IF(prnt) THEN
     title='small matrix'
     CALL prntrm(title,bwrk,end,maxvec,end,maxvec,iout)
  END IF
END SUBROUTINE h_small_d
END MODULE small_matrices
