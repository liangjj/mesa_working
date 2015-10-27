!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE h_m_v
                         USE arnoldi_global
                         USE dvrprop_global
                         USE dvr_shared
                         USE dvr_global
                         USE packed_dvr_matrix_vector_multiply
                         USE fd_m_v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!deck h_m_v_1d_d
!***begin prologue     h_m_v_1d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_1d_d
  SUBROUTINE h_m_v_1d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                    :: nv
  REAL*8, DIMENSION(nphy(1),nv)              :: v_in          
  REAL*8, DIMENSION(nphy(1),nv)              :: v_out          
  REAL*8                                     :: hij
  INTEGER                                    :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_1d_d(v_in,v_out,nv)
  ELSE
     IF(row(1) == 2) THEN
        CALL ds3_bmm_d(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE IF(row(1) == 3) THEN
        CALL ds5_bmm_d(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE IF(row(1) == 4) THEN
        CALL ds7_bmm_d(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
END SUBROUTINE h_m_v_1d_d
!***********************************************************************
!deck h_m_v_1d_z
!***begin prologue     h_m_v_1d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_1d_z
  SUBROUTINE h_m_v_1d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                        :: nv
  COMPLEX*16, DIMENSION(nphy(1),nv)              :: v_in          
  COMPLEX*16, DIMENSION(nphy(1),nv)              :: v_out          
  REAL*8                                         :: hij
  INTEGER                                        :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_1d_z(v_in,v_out,nv)
  ELSE
     IF(row(1) == 2) THEN
        CALL ds3_bmm_z(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE IF(row(1) == 3) THEN
        CALL ds5_bmm_z(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE IF(row(1) == 4) THEN
        CALL ds7_bmm_z(v_in,v_out,1,1,nv,nphy(1),1)
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
END SUBROUTINE h_m_v_1d_z
!***********************************************************************
!deck h_m_v_2d_d
!***begin prologue     h_m_v_2d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_2d_d
  SUBROUTINE h_m_v_2d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                    :: nv
  REAL*8, DIMENSION(nphy(2),nphy(1),nv)      :: v_in          
  REAL*8, DIMENSION(nphy(2),nphy(1),nv)      :: v_out          
  REAL*8                                     :: hij
  INTEGER                                    :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_2d_d(v_in,v_out,nv)
  ELSE
     IF(row(2) == 2) THEN
        CALL ds3_bmm_d(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('error')
        END IF
     ELSE IF(row(2) == 3) THEN
        CALL ds5_bmm_d(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(2) == 4) THEN
        CALL ds7_bmm_d(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
END SUBROUTINE h_m_v_2d_d
!***********************************************************************
!deck h_m_v_2d_z
!***begin prologue     h_m_v_2d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_2d_z
  SUBROUTINE h_m_v_2d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                        :: nv
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)      :: v_in          
  COMPLEX*16, DIMENSION(nphy(2),nphy(1),nv)      :: v_out          
  REAL*8                                         :: hij
  INTEGER                                        :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_2d_z(v_in,v_out,nv)
  ELSE
     IF(row(2) == 2) THEN
        CALL ds3_bmm_z(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('error')
        END IF
     ELSE IF(row(2) == 3) THEN
        CALL ds5_bmm_z(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(2) == 4) THEN
        CALL ds7_bmm_z(v_in,v_out,1,nv,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
  END SUBROUTINE h_m_v_2d_z
!***********************************************************************
!deck h_m_v_3d_d
!***begin prologue     h_m_v_3d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_3d_d
  SUBROUTINE h_m_v_3d_d(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                            :: nv
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nv)      :: v_in      
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nv)      :: v_out     
  REAL*8                                             :: hij
  INTEGER                                            :: i, j, nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_3d_d(v_in,v_out,nv)
  ELSE
     IF(row(3) == 2) THEN
        CALL ds3_bmm_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 3) THEN
        CALL ds5_bmm_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 4) THEN
        CALL ds7_bmm_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_d(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('quit')
     END IF
  END IF
END SUBROUTINE h_m_v_3d_d
!***********************************************************************
!deck h_m_v_3d_z
!***begin prologue     h_m_v_3d_z
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       h_m_v_3d_z
  SUBROUTINE h_m_v_3d_z(v_in,v_out,nv)
  IMPLICIT NONE
  INTEGER                                                :: nv
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1),nv)      :: v_in  
  COMPLEX*16, DIMENSION(nphy(3),nphy(2),nphy(1),nv)      :: v_out  
  REAL*8                                                 :: hij
  INTEGER                                                :: i, j 
  INTEGER                                                :: nel 
  IF(typke /= 'fd') THEN
     CALL dvr_m_v_3d_z(v_in,v_out,nv)
  ELSE
     IF(row(3) == 2) THEN
        CALL ds3_bmm_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 3) THEN
        CALL ds5_bmm_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 4) THEN
        CALL ds7_bmm_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,nv,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt_z(v_in,v_out,1,nv,nphy(1),nphy(2)*nphy(3),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('quit')
     END IF
  END IF
  END SUBROUTINE h_m_v_3d_z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             END MODULE h_m_v
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
