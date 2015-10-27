!======================================================================
      FUNCTION qk1(i1,j1,i2,j2,k)
!======================================================================
!
!                 k
!     Evaluates  Q (i1, j1; i2, j2)
!
!
!     Calls:  Mqk1
!----------------------------------------------------------------------

    USE spline_param
    USE spline_orbitals, p => pbs
    USE spline_integrals
    USE spline_atomic

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i1,j1,i2,j2,k

    ! .. local variables

    INTEGER :: i,j, imin,imax, jmin,jmax, jp,ip
    REAL(KIND=8), DIMENSION(ns,2*ks-1) :: a,b
    REAL(KIND=8) :: qk1,qki

! ... define qk-integrals on B-spline basis

      Call Mqk1(k)

! ... form cross-products 

    Call density (ns,ks,a,p(1,i1),p(1,i2),'n')
    Call density (ns,ks,b,p(1,j1),p(1,j2),'n')

! ... assemble the B-spline orbitals

      qk1 = 0.d0

      do jp = 1,ks+ks-1
       jmin=max( 1, 1 + ks - jp)
       jmax=min(ns,ns + ks - jp)
       do j = jmin,jmax

          qki = 0.d0
          do ip = 1,ks+ks-1
           imin=max( 1, 1 + ks - ip)
           imax=min(ns,ns + ks - ip)
           do i = imin,imax
            qki = qki + a(i,ip)*rkb(i,j,ip,jp) 
           end do
          end do

          qk1 = qk1 + b(j,jp)*qki

       end do
      end do

      END FUNCTION qk1
