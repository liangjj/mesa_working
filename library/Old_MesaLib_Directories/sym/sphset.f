*deck @(#)sphset.f	5.1  11/6/94
      subroutine sphset(maxap3,natoms,a,atmchg,nset,npop,aset,numset,
     $     mset,mpop,init,idx)
      implicit real*8(a-h,o-z)
c
c     a "spherical-set" of atoms is hereby defined as consisting of
c     those atoms which have the same atomic number and which are
c     equidistant from the molecules charge center.  any atoms in a
c     molecule which are equivalent by symmetry must belong to the
c     same spherical-set.
c
c     this routine searches for spherical-sets of atoms.
c     nset(i)    gives the number of the atom in each set where the
c     boundrys between sets can be determined from npop.
c     the list is sorted in terms of increasing distance
c     from the origin and secondarily in terms of
c     increasing atomic number.
c     npop(j)    is the number of atoms in set j.
c     init(j)    is the number of the first atom in set j.
c     aset(i,1)  is atmchg(i) and is also used as a flag.
c     aset(i,2)  is the distance of the i'th atom from the origin.
c
      common/tol/ toler,tol2
      dimension a(maxap3,3),atmchg(1),npop(1),nset(1),aset(maxap3,3)
      dimension mset(1),mpop(1),init(1),idx(1)
      data zero/0.0d+00/
      save zero
c
c     fill aset.
c
c     call rtrace(6hsphset,1)
      do 20 iat=1,natoms
         aset(iat,1) = atmchg(iat)
         aset(iat,2) = sqrt( a(iat,1)**2 + a(iat,2)**2 + a(iat,3)**2 )
 20   continue
c
c     flag any atoms at the origin.
c
      do 40 iat=1,natoms
         if (abs(aset(iat,2)) .lt. toler) aset(iat,1) = zero
 40   continue
c
c     fill mset and mpop.
c
      ic     = 0
      iset   = 0
      do 80 iat=1,natoms
         if (aset(iat,1) .eq. zero) goto 80
         ic   = ic + 1
         iset = iset + 1
         mpop(iset) = 1
         mset(ic)   = iat
         init(iset) = iat
         j1 = iat + 1
         if (iat .eq. natoms) goto 80
         do 60 jat=j1,natoms
            if (aset(jat,1) .eq. zero) goto 60
            if (abs(aset(jat,2)-aset(iat,2)) .gt. toler) goto 60
            ic = ic + 1
            mpop(iset) = mpop(iset) + 1
            mset(ic)   = jat
            aset(jat,1) = zero
 60      continue
 80   continue
      numset = iset
      ictop = ic
c
c     sort the list in terms of increasing distance from the origin.
c     if more than on set is at the same distance place the lower
c     atomic numbered one first.
c
c     the proper ordering is first extablished in idx.
c
      do 100 i=1,numset
 100     idx(i) = i
         if (numset .eq. 1) goto 185
c
         i = 0
 120     i = i + 1
         j = idx(i)
         iat = init(j)
         curd = aset(iat,2)
         curz = aset(iat,1)
         k1 = i + 1
         do 180 k=k1,numset
            l = idx(k)
            jat = init(l)
            if (abs(curd-aset(jat,2)) .gt. toler) goto 140
            if (curz .lt. aset(jat,1)) goto 180
            goto 160
 140        if (curd .lt. aset(jat,2)) goto 180
 160        idx(i) = l
            idx(k) = j
            iat = init(l)
            curd = aset(iat,2)
            curz = aset(iat,1)
 180     continue
         if (i .lt. numset-1) goto 120
c
c     move the data from mset and mpop to nset and npop using the
c     order stored in idx.
c
 185     ic = 0
         do 220 iset=1,numset
            jset = idx(iset)
            npop(iset) = mpop(jset)
            num = npop(iset)
            jc = 0
            if (jset .eq. 1) goto 195
            j2 = jset - 1
            do 190 j=1,j2
 190           jc = jc + mpop(j)
 195           continue
               do 200 i=1,num
                  ic = ic + 1
                  nset(ic) = mset(jc+i)
 200           continue
 220        continue
            return
            end
