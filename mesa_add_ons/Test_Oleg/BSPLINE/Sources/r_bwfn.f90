!======================================================================
      Subroutine R_bwfn (nu)
!======================================================================
!
!     read radial orbitals in B-spline representation from file 'nu'
!
!----------------------------------------------------------------------

      USE spline_atomic
      USE spline_param
      USE spline_orbitals
      
      IMPLICIT NONE
      
      Integer(4), Intent(in) :: nu
      Character(4) :: el
      Real(8) :: zw,hw,hmw,rmw
      Integer(4) :: ksw,nsw,mw,n,l,k,i
      Integer(4), External :: Iadd_bsorb

      rewind(nu)
    1 read(nu,end=2) el,zw,hw,hmw,rmw,ksw,nsw,mw
      if(zw.ne.z) Stop ' R_bwfn:  z <> zw'
      if(abs(hw-h).gt.1.d-12) Stop ' R_bwfn:  h <> hw'
      if(abs(hmw-hmax).gt.1.d-12) Stop ' R_bwfn:  hmw <> hmax'
      if(abs(rmw-rmax).gt.1.d-12) Stop ' R_bwfn:  rmw <> rmax'
      if(ksw.ne.ks) Stop ' R_bwfn:  ksw <> ks'
      if(nsw.ne.ns) Stop ' R_bwfn:  nsw <> ns'
      Call EL4_nlk(el,n,l,k)
      i = Iadd_bsorb(n,l,k)
      read(nu) pbs(1:mw,i)
      mbs(i) = mw
      if(mw.lt.ns) pbs(mw+1:ns,i) = 0.d0 
      go to 1
    2 Continue

      End Subroutine R_bwfn


