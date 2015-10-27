      subroutine mcijij(dj,dkii,dkjj,dkij,nij,nii,njj,niis,njjs,ni,nj,
     $     den)
c
c***begin prologue
c***date written       871022   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c
c***keywords
c***author             lengsfield, byron (brl)
c***source             %W%   %G%
c
c***purpose
c
c***description
c
c***references
c
c***routines called    (none)
c
c***end prologue
c
c     implicit real*8(a-h,o-p,r-z),             integer*2(q)
cc
cmp   extended dummy dj,dkii,dkjj,dkij,den
cc
      dimension dj(nij,nij),dkii(nii,njjs),dkjj(njj,niis),dkij(nij,nij)
      dimension den(2)
      common /number/ zero,pt5,one,two,four,eight
c
c
c
      ic = 0
      do 20 ix = 1, nij
         do 10 jx = 1, ix
            ic = ic + 1
            dj(jx,ix) = den(ic)
            dj(ix,jx) = den(ic)
 10      continue
 20   continue
c
c      dkii
c
      kx = 0
c     jx = 0
      do 140 k = 1, nj
         lx = 0
         jx=k
         do 130 l = 1, nj
c     jx = jx + 1
            ix = 0
            do 120 i = 1, ni
               kxi = kx + i
               do 110 j = 1, i
                  ix = ix + 1
                  dkii(ix,jx) = dj(kxi,lx+j)
 110           continue
 120        continue
            lx = lx + ni
 130     jx=jx+nj
c 130 continue
         kx = kx + ni
 140  continue
c
c     dkjj
c
c     lx = 0
      do 240 i = 1, ni
         lx=i
         do 230 j = 1, ni
c     lx = lx + 1
            ki = i
            kx = 0
            do 220 k = 1, nj
               lj = j
               do 210 l = 1, k
                  kx = kx + 1
                  dkjj(kx,lx) = dj(ki,lj)
                  lj = lj + ni
 210           continue
               ki = ki + ni
 220        continue
 230     lx=lx+ni
c 230 continue
 240  continue
c
c     dkij
c
      jb = 0
      ix = 0
      do 340 j = 1, nj
         do 330 i = 1, ni
            ix = ix + 1
            kb = 0
            jx = 0
            do 320 k = 1, j
               lend = ni
               if (j .eq. k) lend = i
               do 310 l = 1, lend
                  jx = jx + 1
                  dkij(ix,jx) = dj(kb+i,jb+l)
                  dkij(jx,ix) = dj(kb+i,jb+l)
 310           continue
               kb = kb + ni
 320        continue
 330     continue
         jb = jb + ni
 340  continue
c
c     scale the coulomb block
c
      do 500 ij = 1, nij
         dj(ij,ij) = two * dj(ij,ij)
         dkij(ij,ij) = two * dkij(ij,ij)
 500  continue
c
c     dkii
c
      if(nj.eq.1)go to 521
      ix = 0
      do 520 i = 1, ni
         ix = ix + i
         jx = nj
         do 510 j = 2, nj
            j1 = j - 1
            do 505 jj=1,j1
               dkii(ix,jx+jj) = zero
 505        continue
            jx = jx + nj
 510     continue
 520  continue
 521  continue
c
c     dkjj
c
      if(ni.eq.1)go to 621
      ix = 0
      do 620 i = 1, nj
         ix = ix + i
         jx = ni
         do 610 j = 2, ni
            j1 = j - 1
            do 605 jj=1,j1
               dkjj(ix,jx+jj) = zero
 605        continue
            jx = jx + ni
 610     continue
 620  continue
 621  continue
c
      return
      end
