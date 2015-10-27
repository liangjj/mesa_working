!======================================================================
      REAL(8) FUNCTION wk(i1,j1,i2,j2,k,meth)
!======================================================================
!               k
!     Returns  W (i1, j1; i2, j2) base on the assembling the B-spline
!     integrals, which are supposed to be placed in the module
!     spline-integrals. If not, they are calculated by programs
!     mwk_diff or mwk_cell, depending of the parameter 'meth':
!     meth = 'd' - differential equation method
!          = 'c' - cell integration method
!----------------------------------------------------------------------

      USE spline_param
      USE spline_orbitals
      USE spline_integrals
  
      IMPLICIT NONE
      INTEGER(4), INTENT(in) :: i1,j1,i2,j2,k
      CHARACTER(1), INTENT(in), OPTIONAL :: meth
  
      ! .. local variables
  
      INTEGER(4) :: i,ip, j,jp, imin,imax
      REAL(8), DIMENSION(ns,ks+ks-1) :: a
      REAL(8), DIMENSION(ns,ks) :: b
      REAL(8) :: wki
  
      ! .. check the B-spline integrals in module spline-integrals
  
      if(k.ne.krk.or.itype.ne.'wk ') then
!       if(meth.eq.'d') then
!         Call Mwk_diff(k)
!       else
         Call Mwk_cell(k)
!       end if
      end if
  
      ! .. form cross-products
  
      Call density (ns,ks,a,pbs(1,i1),pbs(1,i2),'n')
      Call density (ns,ks,b,pbs(1,j1),pbs(1,j2),'s')
  
      ! .. assembling the B-spline integrals

      wk = 0.d0

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
        wk = wk + b(j,jp)*wki
       end do
      end do

      END FUNCTION wk
  