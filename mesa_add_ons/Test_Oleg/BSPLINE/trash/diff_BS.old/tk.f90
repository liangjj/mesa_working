!======================================================================
      FUNCTION TK(i1,j1,i2,j2,k)
!======================================================================
!                 k
!     Evaluates  T (i1, j1; i2, j2)
!
!     Calls:  MTK
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_integrals
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: i,j, ip,jp, imin,imax, jmin, jmax
    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: tk,tki

    ! .. define Tk-integrals on B-spline basis

      Call MTK(k)

! ... form cross-products

      Call density (ns,ks,a,p(1,i1),p(1,i2),'n')
      Call density (ns,ks,b,p(1,j1),p(1,j2),'n')

! ... assemble the B-spline integrals

        tk = 0.d0

        do jp = 1,ks+ks-1
         jmin=max( 1, 1 + ks-jp)
         jmax=min(ns,ns + ks-jp)
         do j = jmin,jmax

            tki = 0.d0

            do ip = 1,ks+ks-1
              imin=max( 1, 1 + ks-ip)
              imax=min(ns,ns + ks-ip)
              do i = imin,imax

               tki = tki + a(i,ip)*rkb(i,j,ip,jp)

              end do
            end do

            tk = tk + b(j,jp)*tki

         end do
        end do

      END
