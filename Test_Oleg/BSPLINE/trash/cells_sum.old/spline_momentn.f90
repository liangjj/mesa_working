!====================================================================
    MODULE spline_momentc
!====================================================================
!                                                          _
!   contains moments defining as <B_i|r^k|B_j> or <B_i|r^k|B_j> 
!   over an interval, where  _
!                            B = B' - B/r
!
!   These moments are used for calculating of two-electron integrals
!   according to the sell algorithm.
!
!--------------------------------------------------------------------
!
!   rkd(1:ks*ks;1:ks*ks,1:nv) - 
!   
!      the two-dimensional array of integrals <B_i B_j|...|B_i' B_j'>
!      over a triangle (or square) diagonal cell
!     
!   rkd[1,2,3,4](1:ks*ks,1:nv) - different moments defining as
!
!     <B_i|r^k|B_j>  and <B_i|r^k|B_j>  over an interval iv
!    
!   rkd, rkd1, ... differ from rkt, rkt1, ... in module spline_moments)
!   only by reduced dimensions, that increases slightly the speed of
!   calculations
!--------------------------------------------------------------------

    IMPLICIT NONE
    SAVE

    INTEGER :: knk = -100
    CHARACTER(3) :: ntype = 'aaa'

    REAL(8), ALLOCATABLE, DIMENSION(:,:,:) :: rkd
    REAL(8), ALLOCATABLE, DIMENSION(:,:) :: rkd1, rkd2, rkd3, rkd4

    END MODULE spline_momentc



!====================================================================
      SUBROUTINE allocate_momentc
!====================================================================
!
! ... allocates space for arrays in MODULE spline_moments
!
!--------------------------------------------------------------------

      USE spline_param
      USE spline_momentc
   
      INTEGER :: ierr, jk

      if(allocated(rkd)) DEALLOCATE(rkd, rkd1,rkd2,rkd3,rkd4)

      jk = ks*ks
      ALLOCATE(rkd(jk,jk,nv), rkd1(jk,nv),rkd2(jk,nv), &
                              rkd3(jk,nv),rkd4(jk,nv), STAT=ierr)
      if(ierr.ne.0) Stop ' Problem with allocation of moments'
 
      ntype='bbb'
      knk=-100

      END SUBROUTINE allocate_momentc


!====================================================================
      SUBROUTINE dealloc_momentc
!====================================================================
!
!     deallocate space in MODULE spline_moments
!
!--------------------------------------------------------------------

      USE spline_momentc
   
      if(allocated(rkd)) DEALLOCATE(rkd, rkd1,rkd2,rkd3,rkd4)
      ntype='aaa'
      knk=-100

      END SUBROUTINE dealloc_momentc


