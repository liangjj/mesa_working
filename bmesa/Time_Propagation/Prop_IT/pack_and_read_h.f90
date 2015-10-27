!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             MODULE pack_and_read_h
                         USE io
                         USE dvrprop_global_it
                         USE dvr_shared
                         USE dvr_global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         CONTAINS
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
  WRITE(iout,1) dim, count, posble
1 FORMAT (/,1X,'spatial dimension                           = ',i10,/,1x,   &
               'number of non-zero matrix elements          = ',i10,/,1x,   &
               'possible number of non-zero matrix elements = ',i10)
END SUBROUTINE pack_h
!deck rd_hpak.f
!***date written       020206   (yymmdd)
!***revision date               (yymmdd)
!***keywords           pack, dvr
!***author             schneider, b. i.(nsf)
!***source
!***purpose            read and store a packed DVR hamiltonian
!***
!***references
!***routines called    iosys, util and mdutil
!***end prologue       rd_hpak
  SUBROUTINE rd_hpak(diag,hbuf,ihbuf,n,prn)
  IMPLICIT NONE
  INTEGER                                :: n, nel, i
  REAL*8, DIMENSION(n)                   :: diag
  REAL*8, DIMENSION(*)                   :: hbuf
  INTEGER, DIMENSION(2,*)                :: ihbuf
  LOGICAL                                :: prn
  CHARACTER (LEN=80)                     :: title
  CALL iosys('read real "diagonals of dvr hamiltonian" from rwf',    &
              n,diag,0,' ')
  CALL iosys('read integer "number dvr matrix elements" from rwf',1, &
              nel,0,' ')
  CALL iosys('read integer "packed dvr buffer" from rwf',2*nel,      &
              ihbuf,0,' ')
  CALL iosys('read real "packed dvr hamiltonian" from rwf',nel,      &
              hbuf,0,' ')
  IF(prn) THEN
     title='diagonals'
     CALL prntrm(title,diag,n,1,n,1,iout)
     WRITE(iout,1)
     DO  i=1,nel
         WRITE(iout,2) ihbuf(1,i), ihbuf(2,i), hbuf(i)
     END DO
  END IF
1    FORMAT(/,1X,'non-zero matrix elements',/,5X,'I',5X,'J', 5X,'Element')
2    FORMAT(i5,5X,i5,5X,e15.8)
END SUBROUTINE rd_hpak
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                             END MODULE pack_and_read_h
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
