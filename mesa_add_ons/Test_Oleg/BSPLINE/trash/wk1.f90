!======================================================================
      FUNCTION wk1(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  W (i1, j1; i2, j2)
!
!
!     Calls:  Mwk1
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_integrals
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: i,j, imin,imax, jp,ip
    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a
    REAL(KIND=8), DIMENSION(ns,ks) :: b
    REAL(KIND=8) :: wk1,wki

! ... define wk-integrals on B-spline basis

      Call Mwk1(k)

! ... form cross-products 

    Call density (ns,ks,a,p(1,i1),p(1,i2),'n')
    Call density (ns,ks,b,p(1,j1),p(1,j2),'s')

! ... assemble the B-spline orbitals

      wk1 = 0.d0

      do jp = 1,ks
       do j = 1,ns-jp+1

          wki = 0.d0
          do ip = 1,ks+ks-1
            imin=max( 1, 1 + ks-ip)
            imax=min(ns,ns + ks-ip)
            do i = imin,imax
             wki = wki + a(i,ip)*rkb(i,j,ip,jp)
            end do
          end do

          wk1 = wk1 + b(j,jp)*wki

       end do
      end do

      END FUNCTION wk1
