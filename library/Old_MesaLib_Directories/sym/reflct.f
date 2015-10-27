*deck @(#)reflct.f	5.1  11/6/94
      subroutine reflct(maxap3,a,b,natoms,t,ixyz)
      implicit real*8(a-h,o-z)
c
c     the coordinates of the natoms in a are reflected in a plane
c     perpindicular to cartesian axis ixyz using the transformation
c     matrix t.  the reflected coordinates are returned in b.
c
      dimension t(3,3), a(maxap3,3), b(maxap3,3)
      data zero,one/0.0d+00,1.0d+00/
      save zero,one
c
c     call rtrace(6hreflct,1)
      do 20 i = 1, 3
         do 10 j = 1, 3
 10         t(i,j) = zero
 20         t(i,i) = one
            t(ixyz,ixyz) = -one
            call tform(maxap3,t,a,b,natoms)
            return
            end
