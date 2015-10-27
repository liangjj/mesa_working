*deck @(#)ornax.f	5.1  11/6/94
      integer function ornax(maxap3,a,natoms,ian)
      implicit real*8(a-h,o-z)
c
c     return the number of the cartesian axes which passes through
c     the greatest number of atoms.  if there is no one axis which
c     passes through more atoms than the other two then return the
c     number of the axes which passes through the largest number
c     of bonds.  if this is not decisive return a zero.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3), ian(1), cutoff(10), iclass(18), ix(4)
      data cutoff/0.89, 1.73, 2.91, 1.25, 2.53, 1.70, 1.66, 2.4 ,
     $     2.04, 2.53/
      data iclass/1,1,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4/
      data ix/0,1,3,6/, zero,one/0.0d+00,1.0d+00/
      save cutoff,iclass,ix,zero,one
      lind(i,j) = ((max(i,j)*(max(i,j)-1)/2)+min(i,j))
c
c     call rtrace(6hornax ,1)
      nx = 0
      ny = 0
      nz = 0
c
c     determine the number of atoms on each axis.
c
      do 20 iat=1,natoms
         xx = a(iat,1) * a(iat,1)
         yy = a(iat,2) * a(iat,2)
         zz = a(iat,3) * a(iat,3)
         dx = sqrt(yy+zz)
         dy = sqrt(xx+zz)
         dz = sqrt(xx+yy)
         if (abs(dx) .lt. toler) nx = nx + 1
         if (abs(dy) .lt. toler) ny = ny + 1
         if (abs(dz) .lt. toler) nz = nz + 1
 20   continue
c
c     is any one count larger than the other two?
c
      if (nz-ny)  40, 50, 60
 40   if (ny-nx)  70,100, 80
 50   if (nx-nz) 100,100, 70
 60   if (nz-nx)  70,100, 90
 70   ornax = 1
      return
 80   ornax = 2
      return
 90   ornax = 3
      return
c
c     determine the number of bonds cut by each cartesian axis.
c     the criteria are
c     1-- the axis and the line connecting atoms a and b must
c     intersect.
c     2-- the distance between a and b must be less than some
c     cutoff.  in order to establish reasonable cutoffs
c     the atoms have been classified
c     class     atoms
c     h         h, he
c     a1        li, be
c     a2        b, c, n, o, f, ne
c     b         na, mg, al, si, p, s, cl, ar
c     these four classes of atoms produce ten types of bonds.
c     the cutoff is 20  greater than the average standard
c     model bond length for the individual bonds of a given
c     type except for a1-a2 bonds where the cutoff is 10
c     greater than the maximum standard model bond length of
c     that type.
c
 100  continue
c
c     for each pair of bonded atoms determine which axis the bond
c     intersects.
c
      do 200 iat=1,natoms
         do 200 jat=1,iat
            ic = 4
            if(ian(iat).le.18) ic = iclass(ian(iat))
            jc = 4
            if(ian(jat).le.18) jc = iclass(ian(jat))
            dist = sqrt( (a(iat,1)-a(jat,1))**2 +
     $           (a(iat,2)-a(jat,2))**2 +
     $           (a(iat,3)-a(jat,3))**2)
            if(dist.lt.tol2.or.dist.gt.cutoff(lind(ic,jc))) goto 200
            do 190 i1=1,3
               i2 = 1 + mod(i1,3)
               i3 = 1 + mod(i2,3)
               qa2 = a(iat,i2)
               qa3 = a(iat,i3)
               qb2 = a(jat,i2)
               qb3 = a(jat,i3)
c
c     reject on-axis atoms.
c
               tst1 = sqrt(qa2*qa2+qa3*qa3)
               if (tst1 .lt. toler) goto 190
               tst2 = sqrt(qb2*qb2+qb3*qb3)
               if (tst2 .lt. toler) goto 190
c
c     does the i1 axis intersect the line defined by iat and jat?
c     (see crc math tables, 20th ed., p 365)
c
               tst1 = qa3 * (qa2-qb2)
               tst2 = qa2 * (qa3-qb3)
               if (abs(tst1-tst2) .gt. tol2) goto 190
c
c     is the point of intersection between the atoms?
c
               if (abs(qa2).lt.toler .and. abs(qb2).lt.toler) goto 150
               tst1 = sign(one,qa2) + sign(one,qb2)
               if (tst1 .eq. zero) goto 160
               goto 190
 150           tst2 = sign(one,qa3) + sign(one,qb3)
               if (tst2 .ne. zero) goto 190
c
c     increment the appropriate counter.
c
 160           goto(165,170,175),i1
 165           nx = nx + 1
               goto 190
 170           ny = ny + 1
               goto 190
 175           nz = nz + 1
c
 190        continue
 200     continue
c
c     pick the biggest count, if any, and return.
c
         if (nz-ny) 240,250,260
 240     if (ny-nx) 270,300,280
 250     if (nx-nz) 300,300,270
 260     if (nz-nx) 270,300,290
 270     ornax = 1
         return
 280     ornax = 2
         return
 290     ornax = 3
         return
 300     ornax = 0
         return
         end
