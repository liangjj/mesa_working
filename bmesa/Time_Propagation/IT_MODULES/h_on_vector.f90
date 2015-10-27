!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE h_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!**begin prologue     h_on_vector
!**date written       040902   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            module for hamiltonian times vector.
!**description        this module contains all the routines necessary to
!**                   
!**references
!**                   
!**routines called    
!**end prologue       h_on_vector
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             INTERFACE h_v
             MODULE PROCEDURE h_v_1d_d, h_v_2d_d, h_v_3d_d 
                             END INTERFACE h_v
!
!
                             INTERFACE dvr_m_v        
             MODULE PROCEDURE dvr_m_v_1d_d, dvr_m_v_2d_d, dvr_m_v_3d_d
                             END INTERFACE dvr_m_v
!
!
                             INTERFACE dvr_mat_mul
             MODULE PROCEDURE dvr_mat_mul_1d_d, dvr_mat_mul_2d_d,    &
                              dvr_mat_mul_3d_d
                             END INTERFACE dvr_mat_mul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             CONTAINS
!deck h_v_1d_d
!**begin prologue     h_v_1d_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector in one dimension
!**
!**references
!**routines called
!**end prologue       h_v_1d_d
  SUBROUTINE h_v_1d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),nvec)     :: v_in          
  REAL*8, DIMENSION(nphy(1),nvec)     :: v_out          
  INTEGER                             :: i 
  v_out = 0.d0
  IF(typke /= 'fd') THEN
!
     CALL dvr_m_v(v_in,v_out)
!
!\center  For the banded matrices, we perform these as special cases.
!\center  Remember that the last few columns must be treated specially.
  ELSE
     IF(row(1) == 2) THEN
        CALL ds3_bmm(v_in,v_out,1,1,nvec,nphy(1),1)
     ELSE IF(row(1) == 3) THEN
        CALL ds5_bmm(v_in,v_out,1,1,nvec,nphy(1),1)
     ELSE IF(row(1) == 4) THEN
        CALL ds7_bmm(v_in,v_out,1,1,nvec,nphy(1),1)
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
!
! take care of the diagonal v
!
  call d_mul(v_out,v_tot,v_in)
  IF(log_prp(3)) THEN
     title='hamiltonian on vectors'
     CALL prntcm(title,v_out,n3d,nvec,n3d,nvec,iout)
  END IF
END SUBROUTINE h_v_1d_d
!deck h_v_2d_d.f
!**begin prologue     h_v_2d_d
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            space hamiltonian times space*time vector
!**                   in two dimensions
!**references
!**routines called
!**end prologue       h_v_2d_d
  SUBROUTINE h_v_2d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)         :: v_in
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)         :: v_out          
!
!
  v_out = 0.d0
  IF(typke /= 'fd') THEN
     call dvr_m_v(v_in,v_out)
  ELSE
     IF(row(2) == 2) THEN
        CALL ds3_bmm(v_in,v_out,1,nvec,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('error')
        END IF
     ELSE IF(row(2) == 3) THEN
        CALL ds5_bmm(v_in,v_out,1,nvec,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(2) == 4) THEN
        CALL ds7_bmm(v_in,v_out,1,nvec,nphy(1),nphy(2),2)
        IF(row(1) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('error in band size')
     END IF
  END IF
!
! take care of the diagonal v
!
  call d_mul(v_out,v_tot,v_in)
  IF(log_prp(3)) THEN
     title='hamiltonian on vectors'
     CALL prntcm(title,v_out,n3d,nvec,n3d,nvec,iout)
  END IF
END SUBROUTINE h_v_2d_d
!deck h_v_3d_d.f
!***begin prologue     h_v_3d_d
!***date written       011126   (yymmdd)
!***revision date      yymmdd   (yymmdd)
!***keywords           time development
!***author             schneider, barry (nsf)
!***source
!***purpose            space hamiltonian times space*time vector
!***
!***references
!***routines called
!***end prologue       h_v_3d_d
  SUBROUTINE h_v_3d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)  :: v_in        
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)  :: v_out        
!
  v_out = 0.d0
  IF(typke /= 'fd') THEN
     call dvr_m_v(v_in,v_out)
  ELSE
     IF(row(3) == 2) THEN
        CALL ds3_bmm(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(3)*nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(3)*nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 3) THEN
        CALL ds5_bmm(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(2)*nphy(3),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(3)*nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(3)*nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE IF(row(3) == 4) THEN
        CALL ds7_bmm(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),3)
        IF(row(2) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE IF(row(2) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,nvec,nphy(1),nphy(2),nphy(3),2)
        ELSE
           CALL lnkerr('quit')
        END IF
        IF(row(1) == 2) THEN
           CALL ds3_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(3)*nphy(2),1)
        ELSE IF(row(1) == 3) THEN
           CALL ds5_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(3)*nphy(2),1)
        ELSE IF(row(1) == 4) THEN
           CALL ds7_bmmt(v_in,v_out,1,nvec,nphy(1),nphy(3)*nphy(2),1)
        ELSE
           CALL lnkerr('quit')
        END IF
     ELSE
        CALL lnkerr('quit')
     END IF
  END IF
!
! take care of the diagonal v
!
  call d_mul(v_out,v_tot,v_in)
  IF(log_prp(3)) THEN
     title='hamiltonian on vectors'
     CALL prntcm(title,v_out,n3d,nvec,n3d,nvec,iout)
  END IF
END SUBROUTINE h_v_3d_d
!deck dvr_m_v_1d_d.f
!***begin prologue     dvr_m_v_1d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_1d_d
  SUBROUTINE dvr_m_v_1d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),nvec)     :: v_in          
  REAL*8, DIMENSION(nphy(1),nvec)     :: v_out          
  REAL*8                              :: hij
  INTEGER                             :: i, j, nel 
!
!    First add in the diagonal contribution.
!
  DO i=1,nphy(1)
     v_out(i,:) = v_out(i,:) +  buf(1)%d(i) * v_in(i,:)
  END DO
!
! Now loop over the non-zero, off-diagonal eleemnts.
!
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(i,:) = v_out(i,:) + hij * v_in(j,:)
      v_out(j,:) = v_out(j,:) + hij * v_in(i,:)
  END DO
END SUBROUTINE dvr_m_v_1d_d
!deck dvr_m_v_2d_d.f
!***begin prologue     dvr_m_v_2d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_2d_d
  SUBROUTINE dvr_m_v_2d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)     :: v_in, v_out
  REAL*8                                      :: hij
  INTEGER                                     :: i, nel, j
  DO  i=1,nphy(2)
      v_out(i,:,:) = v_out(i,:,:) + buf(2)%d(i) * v_in(i,:,:)
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(i,:,:) = v_out(i,:,:) + hij * v_in(j,:,:)
      v_out(j,:,:) = v_out(j,:,:) + hij * v_in(i,:,:)
  END DO
  DO  i=1,nphy(1)
      v_out(:,i,:) = v_out(:,i,:) + buf(1)%d(i) * v_in(:,i,:)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(:,i,:) = v_out(:,i,:) + hij * v_in(:,j,:)
      v_out(:,j,:) = v_out(:,j,:) + hij * v_in(:,i,:)
  END DO
END SUBROUTINE dvr_m_v_2d_d
!deck dvr_m_v_3d_d.f
!***begin prologue     dvr_m_v_3d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_m_v_3d_d
  SUBROUTINE dvr_m_v_3d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)         :: v_in, v_out
  REAL*8                                                  :: hij
  INTEGER                                                 :: i, nel, j
  do i=1,nphy(3)
     v_out(i,:,:,:) = v_out(i,:,:,:) + buf(3)%d(i) * v_in(i,:,:,:)
  end do
  do nel=1,nonz(3)
      i=buf(3)%hibuf(1,nel)
      j=buf(3)%hibuf(2,nel)
      hij=buf(3)%hbuf(nel)
      v_out(i,:,:,:) = v_out(i,:,:,:) + hij * v_in(j,:,:,:)
      v_out(j,:,:,:) = v_out(j,:,:,:) + hij * v_in(i,:,:,:)
  END DO
  DO  i=1,nphy(2)
      v_out(:,i,:,:) = v_out(:,i,:,:) + buf(2)%d(i) * v_in(:,i,:,:)
  END DO
  DO  nel=1,nonz(2)
      i=buf(2)%hibuf(1,nel)
      j=buf(2)%hibuf(2,nel)
      hij=buf(2)%hbuf(nel)
      v_out(:,i,:,:) = v_out(:,i,:,:) + hij * v_in(:,j,:,:)
      v_out(:,j,:,:) = v_out(:,j,:,:) + hij * v_in(:,i,:,:)
  END DO
  DO  i=1,nphy(1)
      v_out(:,:,i,:) = v_out(:,:,i,:) + buf(1)%d(i) * v_in(:,:,i,:)
  END DO
  DO  nel=1,nonz(1)
      i=buf(1)%hibuf(1,nel)
      j=buf(1)%hibuf(2,nel)
      hij=buf(1)%hbuf(nel)
      v_out(:,:,i,:) = v_out(:,:,i,:) + hij * v_in(:,:,j,:)
      v_out(:,:,j,:) = v_out(:,:,j,:) + hij * v_in(:,:,i,:)
  END DO
END SUBROUTINE dvr_m_v_3d_d 
!deck dvr_mat_mul_1d_d.f
!***begin prologue     dvr_mat_mul_1d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_mat_mul_1d_d
  SUBROUTINE dvr_mat_mul_1d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(1),nvec)     :: v_in          
  REAL*8, DIMENSION(nphy(1),nvec)     :: v_out          
  REAL*8                              :: hij
  INTEGER                             :: i, j, nel 
END SUBROUTINE dvr_mat_mul_1d_d
!deck dvr_mat_mul_2d_d.f
!***begin prologue     dvr_mat_mul_2d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_mat_mul_2d_d
  SUBROUTINE dvr_mat_mul_2d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)     :: v_in          
  REAL*8, DIMENSION(nphy(2),nphy(1),nvec)     :: v_out          
  REAL*8                              :: hij
  INTEGER                             :: i, j, nel 
END SUBROUTINE dvr_mat_mul_2d_d
!deck dvr_mat_mul_3d_d.f
!***begin prologue     dvr_mat_mul_3d_d
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           matrix vector multiply, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            multiply a packed DVR hamiltonian on a vector
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       dvr_mat_mul_3d_d
  SUBROUTINE dvr_mat_mul_3d_d(v_in,v_out)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)     :: v_in          
  REAL*8, DIMENSION(nphy(3),nphy(2),nphy(1),nvec)     :: v_out          
  REAL*8                              :: hij
  INTEGER                             :: i, j, nel 
END SUBROUTINE dvr_mat_mul_3d_d
!deck d_mul
!**begin prologue     d_mul
!**date written       010828   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           m_donal matrix vector
!**author             schneider, barry (nsf)
!**source             
!**purpose            m_donal matrix vector multiply.
!**description        
!**                   
!**references
!**routines called
!**end prologue       d_mul

SUBROUTINE d_mul(v_out,m_d,v_in)
  USE dvrprop_global_it
  USE dvr_global
  USE dvr_shared
  IMPLICIT NONE
  INTEGER                             :: i
  REAL*8, DIMENSION(n3d,nvec)         :: v_in, v_out
  REAL*8, DIMENSION(n3d)              :: m_d
!
  DO i=1,n3d
     v_out(i,:) = v_out(i,:) + m_d(i) * v_in(i,:)
  END DO
END  SUBROUTINE d_mul
! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{DS3_BMM: Multiply a Symmetric, Tri-Diagonal Matrix on a Matrix}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck ds3_bmm
!**begin prologue     ds3_bmm
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**
!**description        the two non-zero elements of the symmetric
!**                   tridiagonal matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmm
  SUBROUTINE ds3_bmm(v_in,v_out,n1,n2,n3,n4,dim)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n4-1
      
!\begin{eqnarray}
! v_{out}(i,j,k) &=&  v_{out}(i,j,k) + m(i,i+1) v_{in}(i+1,j,k)
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!            &=&  v_{out}(i,j,k)     + m(i+1,i) v_{in}(i+1,j,k) \nonumber
!                                    + m(i,i)   v_{in}(i,j,k)
!                                    + m(i,i-1) v_{in}(i-1,j,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                            &
                                   +                             &
                       grid(dim)%ke(2,i) * v_in(i+1,:,:,:)       &
                                   +                             &
                       grid(dim)%ke(1,i) * v_in(i,:,:,:)         &
                                   +                             &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)
  END DO

!\center Do the special case of the first and last row.

      v_out(1,:,:,:)  = v_out(1,:,:,:)                           & 
                                   +                             &
                        grid(dim)%ke(2,1) * v_in(2,:,:,:)        &
                                   +                             &
                        grid(dim)%ke(1,1) * v_in(1,:,:,:)
      v_out(n4,:,:,:) = v_out(n4,:,:,:)                          &
                                   +                             &
                        grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)  & 
                                   +                             &
                        grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)
END SUBROUTINE ds3_bmm
!deck ds3_bmmt.f
!**begin prologue     ds3_bmmt
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**
!**description        the two non-zero elements of the tridiagonal matrix
!**                   are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds3_bmmt
  SUBROUTINE ds3_bmmt(v_in,v_out,n1,n2,n3,n4,dim)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8,     DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!\center In the tridiagonal case all but the first and last rows
!\center "feel" all three elements of the matrix

  DO  i=2,n3-1
      
!\begin{eqnarray}
! v_{out}(j,i,k) &=&  v_{out}(j,i,k) + m(i,i+1) v_{in}(j,i+1,k)
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!              &=&  v_{out}(j,i,k)   + m(i+1,i) v_{in}(j,i+1,k) \nonumber
!                                    + m(i,i)   v_{in}(j,i,k)
!                                    + m(i,i-1) v_{in}(j,i-1,k)
!\end{eqnarray}
!    Since the element (i-1) is not explictly stored in column i, we
!    use the off-diagonal element from column (i-1)
      
      v_out(:,i,:,:) = v_out(:,i,:,:)                              &
                                         +                         &
                       grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)       &
                                         +                         &
                       grid(dim)%ke(1,i)   * v_in(:,i,:,:)         &
                                         +                         &
                       grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special case of the first and last row.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                 &
                                         +                         &
                    grid(dim)%ke(2,1)    * v_in(:,2,:,:)           &
                                         +                         &
                    grid(dim)%ke(1,1)    * v_in(:,1,:,:)
  v_out(:,n3,:,:) = v_out(:,n3,:,:)                                &
                                         +                         &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)        &
                                         +                         &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)
END SUBROUTINE ds3_bmmt
! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{DS5BMM: Symmetric, Penta-Diagonal Matrix a Matrix}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck ds5_bmm
!**begin prologue     ds5_bmm
!**date written       011126   (yymmdd)
!**revision date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmm
  SUBROUTINE ds5_bmm(v_in,v_out,n1,n2,n3,n4,dim)
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8,     DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO  i=3,n4-2
!\begin{eqnarray}
! v_{out}(i,j,k) =  v_{out}(i,j,k) + m(i,i+2) v_{in}(i+2,j,k)
!                                  + m(i,i+1) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!                                  + m(i,i-1) v_{in}(i-1,j,k) \nonumber \\
!                =  v_{out}(i,j,k) + m(i+2,i) v_{in}(i+2,j,k)
!                                  + m(i+1,i) v_{in}(i+1,j,k)
!                                  + m(i,i)   v_{in}(i,j,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(i-1,j,k)
!                                  + m(i,i-2) v_{in}(i-2,j,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)
  END DO
!
!\center Do the special cases of the first two and last two rows.
!
  v_out(1,:,:,:)   =  v_out(1,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(1,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,1) * v_in(3,:,:,:)
  v_out(2,:,:,:)   =  v_out(2,:,:,:)                                    &
                                      +                                 &
                      grid(dim)%ke(2,1) * v_in(1,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(1,2) * v_in(2,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(2,2) * v_in(3,:,:,:)                 &
                                      +                                 &
                      grid(dim)%ke(3,2) * v_in(4,:,:,:)
  v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                 &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)           &
                                      +                                 & 
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)           &
                                      +                                 &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)
  v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                   &
                                      +                                 &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)             &
                                      +                                 &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)           & 
                                      +                                 &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)
END SUBROUTINE ds5_bmm
!deck ds5_bmmt.f
!**begin prologue     ds5_bmmt
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the three non-zero elements of the pentadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i)   &=& hb(1,i) \nonumber \\
!              m(i+1,i) &=& hb(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& hb(3,i) = m(i,i+2)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds5_bmmt
  SUBROUTINE ds5_bmmt(v_in,v_out,n1,n2,n3,n4,dim)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8,     DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!\center In the pentadiagonal case all but the first two and last two rows
!\center "feel" all three elements of the matrix

  DO i=3,n3-2
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k)  \nonumber \\
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k)
!                                  + m(i,i)   v_{in}(j,i,k) \nonumber \\
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!\end{eqnarray}
!    Since the elements (i-1) and (i-2) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1) and (i-2)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                              &  
                                  +                               &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)           &
                                  +                               &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)         &
                                  +                               &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)
  END DO

!\center Do the special cases of the first two and last two rows.

  v_out(:,1,:,:)  = v_out(:,1,:,:)                                &
                                  +                               &
                  grid(dim)%ke(1,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,1) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,1) * v_in(:,3,:,:)
  v_out(:,2,:,:)  = v_out(:,2,:,:)                                &
                                  +                               & 
                  grid(dim)%ke(2,1) * v_in(:,1,:,:)               &
                                  +                               &
                  grid(dim)%ke(1,2) * v_in(:,2,:,:)               &
                                  +                               &
                  grid(dim)%ke(2,2) * v_in(:,3,:,:)               &
                                  +                               &
                  grid(dim)%ke(3,2) * v_in(:,4,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                           &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                             &
                                  +                               &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)         &
                                  +                               &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)       &
                                  +                               &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)
END SUBROUTINE ds5_bmmt
! \documentclass{article}
! \usepackage{graphicx}
! \setkeys{Gin}{width=\linewidth}
! \title{DS7_BMM: Symmetric, Hepta-Diagonal Matrix on a Matrix}
! \author{Barry I. Schneider}
! \date{}
! \def \<{\langle}
! \def \>{\rangle}
! \begin{document}
! \maketitle
!deck ds7_bmm.f
!**begin prologue     ds7_bmm
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as

!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%(4,i) = m(i,i+3)
!\end{eqnarray}

!**references
!**routines called
!**end prologue       ds7_bmm
  SUBROUTINE ds7_bmm(v_in,v_out,n1,n2,n3,n4,dim)
  USE arnoldi_global
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8, DIMENSION(n4,n3,n2,n1)         :: v_in, v_out
  INTEGER                                :: i
!
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix
!
  DO  i=4,n4-3
      
!\begin{eqnarray}
! v_{out}(i,:,:) =  v_{out}(i,:,:) + m(i,i+3) v_{in}(i+3,:,:)
!                                  + m(i,i+2) v_{in}(i+2,:,:)
!                                  + m(i,i+1) v_{in}(i+1,:,:)  \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:) \nonumber \\
!                =  v_{out}(i,:,:) + m(i+3,i) v_{in}(i+3,:,:)
!                                  + m(i+2,i) v_{in}(i+2,:,:)
!                                  + m(i+1,i) v_{in}(i+1,:,:) \nonumber \\
!                                  + m(i,i)   v_{in}(i,:,:)
!                                  + m(i,i-1) v_{in}(i-1,:,:)
!                                  + m(i,i-2) v_{in}(i-2,:,:)
!                                  + m(i,i-3) v_{in}(i-3,:,:)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)

      v_out(i,:,:,:) = v_out(i,:,:,:)                                    &
                                      +                                  &
                       grid(dim)%ke(4,i)   * v_in(i+3,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i)   * v_in(i+2,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(2,i)   * v_in(i+1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(1,i)   * v_in(i,:,:,:)               &
                                      +                                  &
                       grid(dim)%ke(2,i-1) * v_in(i-1,:,:,:)             &
                                      +                                  &
                       grid(dim)%ke(3,i-2) * v_in(i-2,:,:,:)             &
                                      +                                  &   
                       grid(dim)%ke(4,i-3) * v_in(i-3,:,:,:)               
  END DO

!\center Do the special cases of the first three and last three rows.

    v_out(1,:,:,:) = v_out(1,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(1,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,1) * v_in(4,:,:,:)
    v_out(2,:,:,:) = v_out(2,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(2,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(1,2) * v_in(2,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,2) * v_in(4,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(4,2) * v_in(5,:,:,:)
    v_out(3,:,:,:) = v_out(3,:,:,:)                                      &
                                      +                                  &
                   grid(dim)%ke(3,1) * v_in(1,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,2) * v_in(2,:,:,:)                     & 
                                      +                                  &
                   grid(dim)%ke(1,3) * v_in(3,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(2,3) * v_in(4,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(3,3) * v_in(5,:,:,:)                     &
                                      +                                  &
                   grid(dim)%ke(4,3) * v_in(6,:,:,:)
    v_out(n4-2,:,:,:) = v_out(n4-2,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(1,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-3) * v_in(n4-3,:,:,:)            & 
                                      +                                  &
                      grid(dim)%ke(3,n4-4) * v_in(n4-4,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-5) * v_in(n4-5,:,:,:)
    v_out(n4-1,:,:,:) = v_out(n4-1,:,:,:)                                &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(1,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(2,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-3) * v_in(n4-3,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-4) * v_in(n4-4,:,:,:)
    v_out(n4,:,:,:)   = v_out(n4,:,:,:)                                  &
                                      +                                  &
                      grid(dim)%ke(1,n4)   * v_in(n4,:,:,:)              &
                                      +                                  &
                      grid(dim)%ke(2,n4-1) * v_in(n4-1,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(3,n4-2) * v_in(n4-2,:,:,:)            &
                                      +                                  &
                      grid(dim)%ke(4,n4-3) * v_in(n4-3,:,:,:)
END SUBROUTINE ds7_bmm
!deck ds7_bmmt.f
!**begin prologue     ds7_bmmt
!**date written       011126   (yymmdd)
!**rev_insion date      yymmdd   (yymmdd)
!**keywords           time development
!**author             schneider, barry (nsf)
!**source
!**purpose            hamiltonian times space vector
!**description        the four non-zero elements of the heptadiagonal
!**                   matrix are stored by columns as
!\begin{eqnarray}
!              m(i,i) &=& grid(dim)%ke(1,i) \nonumber \\
!              m(i+1,i) &=& grid(dim)%ke(2,i) = m(i,i+1) \nonumber \\
!              m(i+2,i) &=& grid(dim)%ke(3,i) = m(i,i+2) \nonumber \\
!              m(i+3,i) &=& grid(dim)%ke(4,i) = m(i,i+3)
!\end{eqnarray}
!**references
!**routines called
!**end prologue       ds7_bmmt
  SUBROUTINE ds7_bmmt(v_in,v_out,n1,n2,n3,n4,dim)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: n1, n2, n3, n4, dim
  REAL*8,     DIMENSION(n4,n3,n2,n1)     :: v_in, v_out
  INTEGER                                :: i
!\center In the heptadiagonal case all but the first three and last three rows
!\center "feel" all four elements of the matrix

  DO i=4,n3-3
      
!\begin{eqnarray}
! v_{out}(j,i,k) =  v_{out}(j,i,k) + m(i,i+3) v_{in}(j,i+3,k)
!                                  + m(i,i+2) v_{in}(j,i+2,k)
!                                  + m(i,i+1) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k) \nonumber \\
!                =  v_{out}(j,i,k) + m(i+3,i) v_{in}(j,i+3,k)
!                                  + m(i+2,i) v_{in}(j,i+2,k)
!                                  + m(i+1,i) v_{in}(j,i+1,k) \nonumber \\
!                                  + m(i,i)   v_{in}(j,i,k)
!                                  + m(i,i-1) v_{in}(j,i-1,k)
!                                  + m(i,i-2) v_{in}(j,i-2,k)
!                                  + m(i,i-3) v_{in}(j,i-3,k)
!\end{eqnarray}
!    Since the elements (i-1), (i-2) and (i-3) are not explictly stored in column i, we
!    use the off-diagonal element from column (i-1), (i-2) and (i-3)
      
     v_out(:,i,:,:) = v_out(:,i,:,:)                                &
                                      +                             &
                    grid(dim)%ke(4,i)   * v_in(:,i+3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i)   * v_in(:,i+2,:,:)           & 
                                      +                             &
                    grid(dim)%ke(2,i)   * v_in(:,i+1,:,:)           &
                                      +                             &
                    grid(dim)%ke(1,i)   * v_in(:,i,:,:)             &
                                      +                             &
                    grid(dim)%ke(4,i-3) * v_in(:,i-3,:,:)           &
                                      +                             &
                    grid(dim)%ke(3,i-2) * v_in(:,i-2,:,:)           &
                                      +                             &
                    grid(dim)%ke(2,i-1) * v_in(:,i-1,:,:)
  END DO

!\center Do the special cases of the first three and last three rows.

  v_out(:,1,:,:) = v_out(:,1,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(1,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,1) * v_in(:,4,:,:)
  v_out(:,2,:,:) = v_out(:,2,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(2,1) * v_in(:,1,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(1,2) * v_in(:,2,:,:)                  & 
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,2) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,2) * v_in(:,5,:,:)
  v_out(:,3,:,:) = v_out(:,3,:,:)                                   &
                                   +                                &
                 grid(dim)%ke(3,1) * v_in(:,1,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,2) * v_in(:,2,:,:)                  &
                                   +                                &
                 grid(dim)%ke(1,3) * v_in(:,3,:,:)                  &
                                   +                                &
                 grid(dim)%ke(2,3) * v_in(:,4,:,:)                  &
                                   +                                &
                 grid(dim)%ke(3,3) * v_in(:,5,:,:)                  &
                                   +                                &
                 grid(dim)%ke(4,3) * v_in(:,6,:,:)
  v_out(:,n3-2,:,:) = v_out(:,n3-2,:,:)                             &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3,:,:)           &
                                   +                                &
                   grid(dim)%ke(2,n3-2) * v_in(:,n3-1,:,:)          &
                                   +                                &
                   grid(dim)%ke(1,n3-2) * v_in(:,n3-2,:,:)          &
                                   +                                &
                   grid(dim)%ke(2,n3-3) * v_in(:,n3-3,:,:)          &
                                   +                                &
                   grid(dim)%ke(3,n3-4) * v_in(:,n3-4,:,:)          &
                                   +                                &
                   grid(dim)%ke(4,n3-5) * v_in(:,n3-5,:,:)
  v_out(:,n3-1,:,:) = v_out(:,n3-1,:,:)                             &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(1,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(2,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-3) * v_in(:,n3-3,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-4) * v_in(:,n3-4,:,:)
  v_out(:,n3,:,:)   = v_out(:,n3,:,:)                               &
                                   +                                &
                    grid(dim)%ke(1,n3)   * v_in(:,n3,:,:)           &
                                   +                                &
                    grid(dim)%ke(2,n3-1) * v_in(:,n3-1,:,:)         &
                                   +                                &
                    grid(dim)%ke(3,n3-2) * v_in(:,n3-2,:,:)         &
                                   +                                &
                    grid(dim)%ke(4,n3-3) * v_in(:,n3-3,:,:)
END SUBROUTINE ds7_bmmt
!deck pack_h.f
!***begin prologue     pack_h
!***date written       960718   (yymmdd)
!***revision date               (yymmdd)
!***keywords           pack, hamiltonian
!***author             schneider, b. i.(nsf)
!***source
!***purpose            The one-body matrix multiply of the DVR Hamiltonian
!                      on a vector may be performed using a packed form of
!                      the DVR hamiltonian in which only the non-zeros are
!                      used.  An arrays of indices and non-zero elements      
!                      are created in this subroutine.
!***

!***references

!***routines called    iosys, util and mdutil
!***end prologue       pack_h
  SUBROUTINE pack_h(dim)
  USE dvrprop_global_it
  USE dvr_shared
  USE dvr_global
  IMPLICIT NONE
  INTEGER                                :: dim
  REAL*8                                 :: eps= 1.d-15
  INTEGER                                :: posble, i, j, count
  posble=nphy(dim)*(nphy(dim)+1)/2
!
!  The diagonals are stored separately.
!
  DO  i=1,nphy(dim)
      buf(dim)%d(i) = grid(dim)%ke(i,i) 
  END DO
!
! Do the off-diagonals
! 
  count = 0
  DO  i=1,nphy(dim)
      DO  j=1,i-1
          IF(ABS(grid(dim)%ke(i,j)) > eps) THEN
             count = count + 1
             buf(dim)%hbuf(count) = grid(dim)%ke(i,j)
             buf(dim)%hibuf(1,count)=i
             buf(dim)%hibuf(2,count)=j
          END IF
      END DO
  END DO
  nonz(dim)=count
  count = count + nphy(dim)
  WRITE(iout,1) count, posble
1 FORMAT (/,1X,'number of non-zero matrix elements          = ', &
                i10, &
          /,1X,'possible number of non-zero matrix elements = ',i10)
END SUBROUTINE pack_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             END MODULE h_on_vector
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
