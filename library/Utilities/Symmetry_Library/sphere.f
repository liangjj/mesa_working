*deck @(#)sphere.f	5.1  11/6/94
      subroutine sphere(maxap3,natoms,a,b,d,atmchg,nset,
     $     npop,norder,dump,mset1,mpop,init,idx)
      implicit real*8(a-h,o-z)
c
c     this routine is called for spherical top molecules.
c     it has two primary functions
c     1--  the highest order proper rotation axis is found and its
c     value placed in norder.
c     the possibilities are
c     5 for the point groups ih, i
c     4 for the point groups oh, o
c     3 for the point groups td, t, th
c     2--  the molecule's gross orientation is fixed as follows
c     t, td, th  the three mutually perpindicular c2 axes are
c     aligned with the cartesian axes so as to max-
c     imize the z-coordinate of the key atom (defined
c     below).
c     o, oh      the three mutually perpindicular c4 axes are
c     aligned with the cartesian axes so as to max-
c     imize the z-coordiante of the key atom.
c     i, ih      on of the six c5 axes is aligned with the
c     cartesian z-axis so as to maximize the z-coord-
c     inate of the key atom.
c
      logical dump
      common/tol/ toler,tol
      common/io/inp,iout
      dimension a(maxap3,3), b(maxap3,3), d(maxap3,3), atmchg(1),
     $     npop(1), nset(1), mset1(1), mpop(1), init(1), idx(1)
      dimension centr(3), save(3), save2(3), t(3,3)
      data zero/0.0d+00/,half,one,two/0.5d+00,1.0d+00,2.0d+00/
      save zero,half,one,two
 1000 format(' sphere-- key atom',i4)
c
c     call rtrace(6hsphere,1)
      norder = 0
      numatm = natoms + 3
      halfpi = two * atan(one)
      pi     = two * halfpi
c
c     find the spherical sets.
c
      call sphset(maxap3,natoms,a,atmchg,nset,npop,d,numset,mset1,
     $     mpop,init,idx)
c
c     find the key atom -- the lowest numbered atom in the first
c     spherical set, where the sperical sets have been ordered
c     in sphset.
c
      key = natoms
      itop = npop(1)
      do 20 i=1,itop
         key = min(key,nset(i))
 20   continue
      if(dump) write(iout,1000) key
c
c     define the smallest spherical set.  axes will be searched for
c     within thes set.
c
      ioff = 0
      moff = 0
      mpop1 = natoms
      do 40 iset=1,numset
         mpop2 = min(mpop1,npop(iset))
         if (mpop2 .eq. mpop1) goto 40
         mpop1 = mpop2
         mset = iset
         moff = ioff
 40      ioff = ioff + npop(iset)
c
c     pick three atoms from the set selected above.
c
         num3 = 0
         i2 = mpop1 - 2
         j2 = mpop1 - 1
         do 160 i=1,i2
            iat = nset(moff+i)
            j1 = i + 1
            do 160 j=j1,j2
               jat = nset(moff+j)
               k1 = j + 1
               do 160 k=k1,mpop1
                  kat = nset(moff+k)
c
c     has a 4- or 5-fold axis been found in a previous triplet of atoms?
c     if yes  then branch to 120
c     else    search for 4- and 5-fold axes.
c
                  if (norder .gt. 3) goto 120
                  call move(maxap3,a,d,numatm)
                  call tstc5(maxap3,d,b,natoms,atmchg,iat,jat,kat,
     $                 centr,itst)
                  if (itst .eq. 0) goto 60
                  norder = 5
                  savez = abs(d(key,3))
                  save(1) = centr(1)
                  save(2) = centr(2)
                  save(3) = centr(3)
                  goto 160
c
 60               continue
                  call move(maxap3,a,d,numatm)
                  call tstc4(maxap3,d,b,natoms,atmchg,iat,jat,kat,
     $                 centr,itst)
                  if (itst .eq. 0) goto 70
                  norder = 4
                  call move(maxap3,d,a,numatm)
                  goto 160
c
c     no axes have been found or num3 3-fold axis have been found.
c     test for a 3-fold axis
c     if no  then get 3 new atoms
c     else   test num3
c     if zero  then set norder, save centr in save, continue
c     if  one  then save cnetr in save2, continue
c     if  two  then continue
c     continue the search
c
 70               if (num3 .eq. 2) goto 160
                  call move(maxap3,a,d,numatm)
                  call tstc3(maxap3,d,b,natoms,atmchg,iat,jat,kat,
     $                 centr,itst)
                  if (itst .eq. 0) goto 160
                  ihop = num3 + 1
                  goto(80,100), ihop
 80               norder = 3
                  num3 = 1
                  save(1) = centr(1)
                  save(2) = centr(2)
                  save(3) = centr(3)
                  goto 160
 100              num3 = 2
                  save2(1) = centr(1)
                  save2(2) = centr(2)
                  save2(3) = centr(3)
                  goto 160
c
 120              if (norder .eq. 4) goto 140
c
c     one 5-fold axis has already been found.
c     test the current atom triplet for a 5-fold axis
c     if no  then continue the search
c     else   is the z-coordinate of the key atom > savez?
c     if yes  then save it
c     else    continue the search
c
                  call move(maxap3,d,a,numatm)
                  call tstc5(maxap3,d,b,natoms,atmchg,iat,jat,kat,
     $                 centr,itst)
                  if (itst .eq. 0) goto 160
                  curz = abs(d(key,3))
                  if (abs(curz-savez) .lt. toler) goto 160
                  if (savez .gt. curz) goto 160
                  savez = curz
                  save(1) = centr(1)
                  save(2) = centr(2)
                  save(3) = centr(3)
                  goto 160
c
c     one 4-fold axis has already been found.
c     test the current atoms for a second 4-fold axis
c     if no  then continue the search
c     else   rotate about z to align the new c4 with y and branch
c     to 180
c
 140              continue
                  call move(maxap3,a,d,numatm)
                  call tstc4(maxap3,d,b,natoms,atmchg,iat,jat,kat,
     $                 centr,itst)
                  if (itst .eq. 0) goto 160
                  if (abs(centr(3)) .gt. toler) goto 160
                  x = centr(1)
                  y = centr(2)
                  theta =halfpi
                  if (abs(y) .gt. toler) theta = atan(x/y)
                  call rotate(maxap3,d,a,numatm,t,3,theta)
                  goto 180
c
 160           continue
c
c     branch on the value of the highest order axis.  return if no
c     axis was found.
c
 180           continue
               if (norder .eq. 0) return
               ihop = norder - 2
               goto(200,220,280), ihop
c
c     norder = 3, point groups t, td, th.
c     vectors coincident with two of the 3-fold axes are in save and
c     save2.  align the c2 which bisects these with z.
c
 200           continue
               centr(1) = half * ( save(1) + save2(1) )
               centr(2) = half * ( save(2) + save2(2) )
               centr(3) = half * ( save(3) + save2(3) )
               call putsym(maxap3,a,b,t,centr,numatm,3)
               b(1,1) = save(1)
               b(1,2) = save(2)
               b(1,3) = save(3)
               b(2,1) = save2(1)
               b(2,2) = save2(2)
               b(2,3) = save2(3)
               call putsym(maxap3,b,d,t,centr,2,3)
               save(1) = b(1,1)
               save(2) = b(1,2)
               save(3) = b(1,3)
               save2(1) = b(2,1)
               save2(2) = b(2,2)
               save2(3) = b(2,3)
c
c     find a second c2 and align it with y.
c
               x = half * ( save(1) - save2(2) )
               y = half * ( save(2) + save2(1) )
               theta = halfpi
               if (abs(y) .gt. toler) theta = atan(x/y)
               call rotate(maxap3,a,b,numatm,t,3,theta)
c
c     put that c2 on z which will maximize the z-coordinate of the
c     key atom.
c
               call move(maxap3,b,a,numatm)
c
c     norder = 4, point groups o, oh.
c     the three c4 axes are aligned with the cartesian axes.  align that
c     c4 with z that will maximize the z-coordinate of the key atom.
c
 220           continue
               x = a(key,1)
               y = a(key,2)
               z = a(key,3)
               cmax = abs(z)
               ixyz = 3
               do 240 i=1,2
                  tst = abs(a(key,i))
                  if (abs(tst-cmax) .lt. toler) goto 240
                  if (cmax .gt. tst) goto 240
                  ixyz = i
                  cmax = tst
 240           continue
               if (ixyz .eq. 3) goto 260
               ixyz = iabs(ixyz-2) + 1
               call rotate(maxap3,a,b,numatm,t,ixyz,halfpi)
               call move(maxap3,b,a,numatm)
 260           if (a(key,3) .gt. zero) return
               call rotate(maxap3,a,b,numatm,t,1,pi)
               call move(maxap3,b,a,numatm)
               return
c
c     norder = 5, point groups i, ih.
c     the c5 axis to be aligned with z is indicated by save.
c
 280           continue
               call putsym(maxap3,a,b,t,save,numatm,3)
               goto 260
               end
