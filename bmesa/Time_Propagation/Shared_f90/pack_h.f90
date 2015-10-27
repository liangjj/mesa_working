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
  USE arnoldi_global
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
!  title='kinetic energy'
!  call prntrm(title,grid(dim)%ke,nphy(dim),nphy(dim),nphy(dim),nphy(dim),iout)
!  write(iout,*) 'the diagonals'
!  write(iout,*) buf(dim)%d
!  write(iout,*) 'the off-diagonals'
!  do i=1,nonz(dim)
!     write(iout,*) buf(dim)%hibuf(1,i),buf(dim)%hibuf(2,i),buf(dim)%hbuf(i)
!  end do
1 FORMAT (/,1X,'number of non-zero matrix elements          = ', &
                i10, &
          /,1X,'possible number of non-zero matrix elements = ',i10)
END SUBROUTINE pack_h
