!======================================================================
      FUNCTION WK(i1,j1,i2,j2,k)
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
    REAL(KIND=8) :: vk,vki

    ! .. define vk-integrals on B-spline basis

      Call MVK(k)

!     form cross-products for symmetric b-array

      b = 0.d0
!                                           central diagonal:
      jp = 1
      do jj =1,ns
        b(jj,jp) = p(jj,j1)*p(jj,j2)
      end do
!                                           upper  diagonals:
      do jp = 2,ks
        do jj = 1,ns-jp+1
          b(jj,jp) = p(jj,j1)*p(jj+jp-1,j2) + p(jj,j2)*p(jj+jp-1,j1)
        end do
      end do

! ..  form cross-products for non-symmetric a-array

      a = 0.d0
      do j = 1,ks+ks-1
        imin=max0( 1, 1 + ks-j)
        imax=min0(ns,ns + ks-j)
        do i = imin,imax
          a(i,j) = p(i,i1)*p(i+j-ks,i2)
        end do
      end do

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

            k = vk + b(jj,jp)*vki

         end do
        end do

      END FUNCTION WK
