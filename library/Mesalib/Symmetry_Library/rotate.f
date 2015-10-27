*deck @(#)rotate.f	5.1  11/6/94
      subroutine rotate(maxap3,a,b,natoms,t,ixyz,theta)
      implicit real*8(a-h,o-z)
c
c     the coordinates of the natoms in a are rotated counterclockwise
c     by an angle theta around the cartesian axis ixyz using the
c     transformation matrix t.  the rotated coordinates are returned
c     in b.
c
      dimension a(maxap3,3), b(maxap3,3), t(3,3)
      data zero,one/0.d0,1.d0/
      save zero,one
c
c     call rtrace(6hrotate,1)
      i1 = ixyz
      i2 = 1 + mod(i1,3)
      i3 = 1 + mod(i2,3)
      s = sin(theta)
      c = cos(theta)
      t(i1,i1) = one
      t(i1,i2) = zero
      t(i1,i3) = zero
      t(i2,i1) = zero
      t(i2,i2) = c
      t(i2,i3) = s
      t(i3,i1) = zero
      t(i3,i2) = -s
      t(i3,i3) = c
      call tform(maxap3,t,a,b,natoms)
      return
      end
