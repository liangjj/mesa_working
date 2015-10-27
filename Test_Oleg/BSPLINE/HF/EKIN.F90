!=======================================================================
   FUNCTION ekin(i,ii,rel) 
!=======================================================================
!
!       returns the value of the integral of
!
!         (2/r)p (y p  + x )
!               j  i i    i
!
!   integrated with respect to r.
!
!----------------------------------------------------------------------
!
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i,ii
    LOGICAL, INTENT(IN) :: rel
    REAL(KIND=8) :: ekin

    REAL(KIND=8), DIMENSION(ns,nt) :: yrm, xv, w

      call xch(i,rel,xv)
      call potl(i,rel,yrm)
      call bxvpw(nt,ks,ns,2.0d0,yrm,p(1,i),xv,w)
      ekin = SUM(p(:,ii)*xv)
      print *, 'Ekin', i,ii,ekin
   END FUNCTION ekin
