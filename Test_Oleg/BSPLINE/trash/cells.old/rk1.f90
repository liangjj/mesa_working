!======================================================================
    FUNCTION RK1 (i1,j1,i2,j2,k)
!======================================================================
!             k
!   Returns  R (i1, j1; i2, j2)
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
    REAL(KIND=8) :: rk1,rki

    ! .. define spline integrals

    Call MRK1(k)

    ! .. form cross-products

    Call density (ns,ks,a,p(1,i1),p(1,i2),'s')
    Call density (ns,ks,b,p(1,j1),p(1,j2),'s')

    ! .. assemble the spline integrals

    rk1 = 0.d0
    do jp = 1,ks
      do j = 2,ns-jp+1
        rki = 0.d0
        do ip = 1,ks
          do i = 2,ns-ip+1
            rki = rki+a(i,ip)*rkb(i,j,ip,jp)
          end do
        end do
        rk1 = rk1 + b(j,jp)*rki
      end do
    end do


    END FUNCTION RK1
