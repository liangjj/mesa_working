      subroutine mciijj(djii,dk,nii,njj,nij,ni,nj,den,ipair)
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
cmp   extended dummy djii,dk,den,ipair
cc
      dimension djii(nii,njj),dk(nij,nij)
      dimension den(2),ipair(2)
      common /number/ zero,pt5,one,two,four,eight
c
c        note   j > i   so exchange block is dkji
c
      ic = 0
      do 20 jx = 1, njj
         do 10 ix = 1, nii
            ic = ic + 1
            djii(ix,jx) = den(ic)
 10      continue
 20   continue
c
c     build exchange matrix
c
      ji = 0
      do 60 i = 1, ni
         ipi = ipair(i)
         do 50 j = 1, nj
            ji = ji + 1
            ipj = ipair(j)
            jjii = 0
            do 40 ii = 1, i
               ipii = ipi + ii
               do 30 jj = 1, j
                  jjii = jjii + 1
                  dk(ji,jjii) = djii(ipii,ipj+jj)
                  dk(jjii,ji) = djii(ipii,ipj+jj)
 30            continue
               if (ii .eq. i) go to 40
               if (j .eq. nj) go to 40
               j1 = j + 1
               do 25 jj = j1, nj
                  jjii = jjii + 1
                  dk(ji,jjii) = djii(ipii,ipair(jj)+j)
                  dk(jjii,ji) = djii(ipii,ipair(jj)+j)
 25            continue
 40         continue
 50      continue
 60   continue
c
c    scale    dk(jp,jq)      p>q
c
      ix = nj
      do 140 i = 2, ni
         i1 = i - 1
         iixx = 0
         do 130 ii = 1, i1
            do 120 j = 1, nj
               dk(ix+j,iixx+j) = two * dk(ix+j,iixx+j)
               dk(iixx+j,ix+j) = dk(ix+j,iixx+j)
 120        continue
            iixx = iixx + nj
 130     continue
         ix = ix + nj
 140  continue
c
c    scale    dk(pi,qi)      p>q
c
      ix = 0
      do 180 i = 1, ni
         if (nj .eq. 1) go to 165
         do 160 j = 2, nj
            j1 = j - 1
            do 150 jj = 1, j1
               dk(ix+j,ix+jj) = two * dk(ix+j,ix+jj)
               dk(ix+jj,ix+j) = dk(ix+j,ix+jj)
 150        continue
 160     continue
 165     continue
         do 170 j = 1, nj
            dk(ix+j,ix+j) = four * dk(ix+j,ix+j)
 170     continue
         ix = ix + nj
 180  continue
c
      return
      end
