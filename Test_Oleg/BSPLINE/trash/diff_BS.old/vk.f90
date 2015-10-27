!======================================================================
      FUNCTION VK(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  V (i1, j1; i2, j2)
!
!
!     Calls:  MVK
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_integrals

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: i,j, imin,imax, jp,jj,ip,ii, iimin, iimax
    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a
    REAL(KIND=8), DIMENSION(ns,ks) :: b
    REAL(KIND=8) :: Vk,Vki

! .. define vk-integrals on B-spline basis

      Call MVK(k)

! ... form cross-products

      Call density (ns,ks,a,p(1,i1),p(1,i2),'n')
      Call density (ns,ks,b,p(1,j1),p(1,j2),'s')

! ... assemble the B-spline integrals

        vk = 0.d0

        do jp = 1,ks
         do jj = 1,ns-jp+1

            vki = 0.d0
            do ip = 1,ks+ks-1
              iimin=max( 1, 1 + ks-ip)
              iimax=min(ns,ns + ks-ip)
              do ii = iimin,iimax
               vki = vki + a(ii,ip)*rkb(ii,jj,ip,jp)
              end do
            end do

            vk = vk + b(jj,jp)*vki

         end do
        end do

      END FUNCTION VK
