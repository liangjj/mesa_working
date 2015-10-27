*deck @(#)orcn.f	5.1  11/6/94
      subroutine orcn(maxap3,a,b,d,aset,atmchg,npop,nset,natoms,dump,
     $     savang,angmat)
      implicit real*8(a-h,o-z)
c
c     this routine orients molecules in the groups cs, cn, sn, cnh,
c     and i.  these point groups have one uniquely defined axis which
c     is coincident with the z cartesian axis.  rotation about this
c     axis, however, does not change the position of any symmetry
c     elements.  thus, this routine rotates the molecule about the
c     z-axis so as to maximize the number of pairs of heavy (non-
c     hydrogen) atoms which are parallel to the y axis.  if there is
c     only one heavy atom or if there is no orientation in which
c     more than one heavy atom pair can be aligned with y, then the
c     key atom as defined in function subroutine orkey is put in the
c     yz plane so as to give it a positive y coordinate.  if two or
c     more orientations give the same number of "bonds" parallel to y
c     then the one which maximizes the y coordinate of the key atom
c     is selected.
c
      integer orkey
      logical dump
      common/io/inp,iout
      common/tol/ toler,tol2
      dimension a(maxap3,3), b(maxap3,3), d(maxap3,3), atmchg(1),
     $     npop(1), nset(1), aset(maxap3,3), angmat(1),
     $     savang(1), t(3,3)
      data zero,one,two,flag/0.0d+00,1.0d+00,2.0d+00,100.0d+00/
      data heavy/2.0d+00/
      save zero,one,two,flag,heavy
 1000 format(17h orcn-- key atom ,i3)
 1010 format(39h orcn-- heavy atom matrix truncated at ,i3)
c
c     call rtrace(6horcn  ,1)
      numatm = natoms + 3
      halfpi = two * atan(one)
      pi     = two * halfpi
      key    = orkey(maxap3,natoms,a,atmchg,nset,npop,aset)
      if(dump) write(iout,1000) key
c
c     calculate the angles of rotation neccessary to align each heavy
c     atom pair with y.
c
      nhvy = 0
      idx  = 0
      do 60 iat=1,natoms
         if (atmchg(iat) .le. heavy) goto 60
         x1 = a(iat,1)
         y1 = a(iat,2)
         j2 = iat - 1
         nhvy = nhvy + 1
         if (nhvy .eq. 1) goto 60
         do 40 jat=1,j2
            if (atmchg(jat) .le. heavy) goto 40
            idx = idx + 1
            x = a(jat,1) - x1
            y = a(jat,2) - y1
            theta = halfpi
            if (abs(y) .gt. toler) theta = - atan(x/y)
            angmat(idx) = theta
 40      continue
 60   continue
c
c     which angle occurs most frequently?
c
 80   if (nhvy .eq. 1) goto 280
      if (nhvy .gt. 2) goto 100
      savang(1) = angmat(1)
      nsav = 1
      goto 200
 100  continue
      i2 = nhvy * (nhvy-1) / 2
      nmax = 0
      nsav = 0
      do 180 i=1,i2
         curang = angmat(i)
         if (curang .eq. flag) goto 180
         j1 = i + 1
         ncur = 1
         if(j1.gt.i2)go to 125
         do 120 j=j1,i2
            if (abs(curang-angmat(j)) .gt. tol2) goto 120
            ncur = ncur + 1
            angmat(j) = flag
 120     continue
 125     continue
         if (nmax-ncur) 140,160,180
 140     nsav = 1
         savang(1) = curang
         nmax = ncur
         goto 180
 160     nsav = nsav + 1
         savang(nsav) = curang
 180     if (nmax .eq. 1) goto 280
c
c     if nsav is one then a unique orientation has been selected.
c     rotate the molecule and return.
c
 200     if (nsav .gt. 1) goto 220
         call rotate(maxap3,a,b,numatm,t,3,savang(1))
         call move(maxap3,b,a,numatm)
         call oraxis(maxap3,a,b,natoms,atmchg,2)
         return
c
c     find which of several orientations maximizes the y-coordinate of
c     the key atom.  for zero or equal y values the x-coordinate is
c     tested.
c
 220     d(1,1) = a(key,1)
         d(1,2) = a(key,2)
         d(1,3) = a(key,3)
         call rotate(maxap3,d,b,1,t,3,savang(1))
         bestx = b(1,1)
         besty = b(1,2)
         ibest = 1
         do 260 i=2,nsav
            call rotate(maxap3,d,b,1,t,3,savang(i))
            if (abs(abs(besty)-abs(b(1,2))) .lt. toler) goto 240
            if (abs(besty).lt.toler .and. abs(b(1,2)).lt.toler) goto 240
            if (besty .gt. b(1,2)) goto 260
            ibest = i
            besty = b(1,2)
            bestx = b(1,1)
            goto 260
 240        if (abs(abs(bestx)-abs(b(1,1))) .lt. toler) goto 260
            if (bestx .gt. b(1,1)) goto 260
            ibest = i
            bestx = b(1,1)
            besty = b(1,2)
 260     continue
c
         call rotate(maxap3,a,b,numatm,t,3,savang(ibest))
         call move(maxap3,b,a,numatm)
         return
c
c     no orientation aligns more than one heavy atom pair with y.
c
 280     continue
         theta = halfpi
         x = a(key,1)
         y = a(key,2)
         if (abs(y) .gt. toler) theta = - atan(x/y)
         call rotate(maxap3,a,b,numatm,t,3,theta)
         if (b(key,2) .gt. zero) goto 300
         call rotate(maxap3,b,a,numatm,t,3,pi)
         goto 320
 300     call move(maxap3,b,a,numatm)
c
c     planar molecules contained in the xy plane have not been
c     completely specified.
c
 320     continue
         call orptst(maxap3,a,natoms,ixyz)
         if (ixyz .eq. 3) call oraxis(maxap3,a,b,natoms,atmchg,1)
         return
         end
