!======================================================================
      FUNCTION NK (I1,J1,I2,J2,K)
!======================================================================
!                 k
!     Evaluates  N (i1, j1; i2, j2) with diff. eq. method
!
!     Calls: MNK
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_integrals
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: jp,j,ip,i
    REAL(KIND=8), DIMENSION(ns,ks) :: a,b
    REAL(KIND=8) :: nk,rki

    ! .. define nk-integrals on B-spline basis

      Call MNK(k)

    ! ..    form cross-products

      Call density (ns,ks,a,p(1,i1),p(1,i2),'s')
      Call density (ns,ks,b,p(1,j1),p(1,j2),'s')

        nk = 0.d0

        do jp = 1,ks
         do j = 1,ns-jp+1

              rki = 0.d0
              do ip = 1,ks
               do i = 1,ns-ip+1
                 rki = rki + a(i,ip)*rkb(i,j,ip,jp)
               end do
              end do

              nk = nk + b(j,jp)*rki

         end do
        end do


      END FUNCTION NK