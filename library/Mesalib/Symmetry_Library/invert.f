*deck @(#)invert.f	5.1  11/6/94
      subroutine invert(maxap3,a,b,natoms,t)
      implicit real*8(a-h,o-z)
c
c     the coordinates of the natoms in a are inverted through the
c     origin using the transformation matrix t.  the inverted
c     coordinates are returned in b.
c
      dimension a(maxap3,3), b(maxap3,3), t(3,3)
      data zero,one/0.0d+00,1.0d+00/
      save zero,one
c
c     call rtrace(6hinvert,1)
      do 50 i=1,3
         do 50 j=1,3
 50         t(i,j) = zero
            t(1,1) = - one
            t(2,2) = - one
            t(3,3) = - one
            call tform(maxap3,t,a,b,natoms)
            return
            end
