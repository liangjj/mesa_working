*deck @(#)secmom.f	5.1  11/6/94
      subroutine secmom(maxap3,natoms,a,atmchg,eigval,eigvec)
      implicit real*8(a-h,o-z)
c
c     compute the principal second moments of charge and axes of
c     charge for the natoms whose coordinates are in a and
c     atomic numbers are in atmchg.
c
      dimension a(maxap3,3), atmchg(1), t(6), e(9), e2(18), eigval(1),
     $     eigvec(3,3)
      data zero/0.0d+00/
      save zero
c
c     call rtrace(6hsecmom,1)
      do 20 i=1,6
         t(i) = zero
 20   continue
c
      do 40 iat=1,natoms
         an = atmchg(iat)
         x  = a(iat,1)
         y  = a(iat,2)
         z  = a(iat,3)
         t(1) = t(1) + an * (y*y+z*z)
         t(3) = t(3) + an * (x*x+z*z)
         t(6) = t(6) + an * (x*x+y*y)
         t(2) = t(2) - an * x * y
         t(4) = t(4) - an * x * z
         t(5) = t(5) - an * y * z
 40   continue
      call rsp(3,3,6,t,eigval,1,eigvec,e,e2,ierr)
      return
      end
