!======================================================================
    FUNCTION mk1 (i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  M (i1, j1; i2, j2)
!
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_integrals

    IMPLICIT NONE
    INTEGER, INTENT(in) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: jp,j,ip,i
    REAL(KIND=8), DIMENSION(ns,ks) :: a,b
    REAL(KIND=8) :: mk1,mki

    ! .. define spline integrals

    Call Mmk1(k)

    ! .. form cross-products

    Call density (ns,ks,a,p(1,i1),p(1,i2),'s')
    Call density (ns,ks,b,p(1,j1),p(1,j2),'s')

    ! .. assemble the spline integrals

    mk1 = 0.d0
    do jp = 1,ks
      do j = 1,ns-jp+1
        mki = 0.d0
        do ip = 1,ks
          do i = 1,ns-ip+1
            mki = mki+a(i,ip)*rkb(i,j,ip,jp)
          end do
        end do
        mk1 = mk1 + b(j,jp)*mki
      end do
    end do

    mk1 = mk1 / 2.d0

    END FUNCTION mk1
