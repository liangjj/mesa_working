*deck @(#)putsym.f	5.1  11/6/94
      subroutine putsym(maxap3,a,b,t,v,natoms,ixyz)
      implicit real*8(a-h,o-z)
c
c     the n sets of coordinates in a are rotated so that the point
c     specified by v is placed on the ixyz axis.  the axis of rotation
c     is given by the vector product of v with the unit vector defining
c     the axis ixyz.  the angle of rotation is given by the appropriate
c     direction cosine.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), v(3), t(3,3), t1(3,3)
      data one/1.0d+00/
      save one
c
c     define the axes i1, i2, and i3 where i1=ixyz:
c     i1:  x  z  y
c     i2:  y  x  z
c     i3:  z  y  x
c     the projections on these (i1, i2, i3) are v1, v2, and v3.
c     call rtrace(6hputsym,1)
c
      i1 = ixyz
      i2 = 1 + mod(i1,3)
      i3 = 1 + mod(i2,3)
      v1 = v(i1)
      v2 = v(i2)
      v3 = v(i3)
      vnorm = sqrt(v1*v1 + v2*v2 + v3*v3)
      if (abs(abs(v1)-vnorm) .lt. toler) return
c
c     compute the direction cosines and some common factors.
c
      alph = v1 / vnorm
      beta = v2 / vnorm
      gamm = v3 / vnorm
      v2v2 = v2 * v2
      v3v3 = v3 * v3
      v2233 = one / (v2v2+v3v3)
c
c     form the transformation matrix in the i1, i2, i3 coordinates.
c     this matrix will place the point (v1,v2,v3) on the pos've i1 axis.
c
      t1(1,1) =  alph
      t1(1,2) =  beta
      t1(1,3) =  gamm
      t1(2,1) =  - t1(1,2)
      t1(3,1) =  - t1(1,3)
      t1(2,3) =  v2 * v3 * (alph-one) * v2233
      t1(3,2) =  t1(2,3)
      t1(2,2) =  (v2v2*alph + v3v3) * v2233
      t1(3,3) =  (v3v3*alph + v2v2) * v2233
c
c     transform to the original coordinate system.
c
      t(i1,i1) = t1(1,1)
      t(i1,i2) = t1(1,2)
      t(i1,i3) = t1(1,3)
      t(i2,i1) = t1(2,1)
      t(i2,i2) = t1(2,2)
      t(i2,i3) = t1(2,3)
      t(i3,i1) = t1(3,1)
      t(i3,i2) = t1(3,2)
      t(i3,i3) = t1(3,3)
c
c     carry out the rotation.
c
      call tform(maxap3,t,a,b,natoms)
      call move(maxap3,b,a,natoms)
      return
      end
