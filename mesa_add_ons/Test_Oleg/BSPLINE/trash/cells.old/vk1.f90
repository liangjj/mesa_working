!======================================================================
      FUNCTION VK1(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  V (i1, j1; i2, j2)
!
!
!     Calls:  MVK1
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals
    USE spline_integrals
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: i,j, ip, jp, imin, imax
    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a
    REAL(KIND=8), DIMENSION(ns,ks) :: b
    REAL(KIND=8) :: vk1,vki

! ... define vk-integrals on B-spline basis

      Call MVK1(k)

! ... form cross-products

      Call density (ns,ks,a,pbs(1,i1),pbs(1,i2),'n')
      Call density (ns,ks,b,pbs(1,j1),pbs(1,j2),'s')

! ... assemble the B-spline integrals

      vk1 = 0.d0

      do jp = 1,ks
       do j = 1,ns-jp+1

          vki = 0.d0
          do ip = 1,ks+ks-1
            imin=max( 1, 1 + ks-ip)
            imax=min(ns,ns + ks-ip)
            do i = imin,imax
             vki = vki + a(i,ip)*rkb(i,j,ip,jp)
            end do
          end do

          vk1 = vk1 + b(j,jp)*vki

       end do
      end do

      END FUNCTION VK1
