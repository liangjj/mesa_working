!*deck exp_off_diag_v_m_3_d
!***begin prologue     exp_off_diag_v_m_3_d     
!***date written       031122   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time, dvr, real-space, split-operator, propagation
!***author             schneider, b, i. (nsf)
!***source             FEDVR_Prop
!***purpose            Time dependent Schroedinger equation using split
!***                   operator with finite difference or DVR
!
!***description        This routine simply calls
!***                   exp_off_diag_v_m_2_d in a loop over a dummy
!***                   index, k
!
!***references
!***routines called
!***end prologue       exp_off_diag_v_m_3_d
!
  SUBROUTINE exp_off_diag_v_m_3_d(real_v,      &
                                 imag_v,      &
                                 real_v_scr,  &
                                 imag_v_scr,  &
                                 ni,nj,nk,index)
  IMPLICIT NONE
  INTEGER                             :: ni, nj, nk, index
  REAL*8, DIMENSION(ni,nj,nk)         :: real_v, imag_v
  REAL*8, DIMENSION(ni,nj,nk)         :: real_v_scr, imag_v_scr
  INTEGER                             :: k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              
  DO k=1,nk
     CALL exp_off_diag_v_m_2_d(real_v(1,1,k),          &
                               imag_v(1,1,k),          &
                               real_v_scr(1,1,k),      &
                               imag_v_scr(1,1,k),      &
                               ni,nj,index) 
  END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
END SUBROUTINE exp_off_diag_v_m_3_d
