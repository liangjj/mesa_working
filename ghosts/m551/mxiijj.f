*deck %W%  %G%
      subroutine mxiijj(djii,dk,nii,njj,nij,ni,nj,den,ipair)
C
C***Begin prologue
C***Date written       871022   (yymmdd)
C***Revision date      yymmdd   (yymmdd)
C
C***Keywords
C***Author             Lengsfield, Byron (BRL)
C***Source             %W%   %G%
C
C***Purpose
C
C***Description
C
C***References
C
C***Routines called    (none)
C
C***End prologue
C
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy djii,dk,den,ipair
cc
      dimension djii(nii,njj),dk(nij,nij)
      dimension den(2),ipair(2)
      common /number/ zero,pt5,one,two,four,eight
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
      ij = 0
      do 60 j = 1, nj
         ipj = ipair(j)
         do 50 i = 1, ni
            ij = ij + 1
            ipi = ipair(i)
            iijj = 0
            do 40 jj = 1, j
c     ipjj = ipj + j
               ipjj = ipj + jj
               do 30 ii = 1, i
                  iijj = iijj + 1
                  dk(ij,iijj) = djii(ipi+ii,ipjj)
                  dk(iijj,ij) = djii(ipi+ii,ipjj)
 30            continue
               if (jj .eq. j) go to 40
               if (i .eq. ni) go to 40
               i1 = i + 1
               do 25 ii = i1, ni
                  iijj = iijj + 1
                  dk(ij,iijj) = djii(ipair(ii)+i,ipjj)
                  dk(iijj,ij) = djii(ipair(ii)+i,ipjj)
 25            continue
 40         continue
 50      continue
 60   continue
c
c    scale    dk(ip,iq)      p>q
c
      jx = ni
      do 140 j = 2, nj
         j1 = j - 1
         jjxx = 0
         do 130 jj = 1, j1
            do 120 i = 1, ni
               dk(jx+i,jjxx+i) = two * dk(jx+i,jjxx+i)
               dk(jjxx+i,jx+i) = dk(jx+i,jjxx+i)
 120        continue
            jjxx = jjxx + ni
 130     continue
         jx = jx + ni
 140  continue
c
c    scale    dk(pj,qj)      p>q
c
      jx = 0
      do 180 j = 1, nj
         if (ni .eq. 1) go to 165
         do 160 i = 2, ni
            i1 = i - 1
            do 150 ii = 1, i1
               dk(jx+i,jx+ii) = two * dk(jx+i,jx+ii)
               dk(jx+ii,jx+i) = dk(jx+i,jx+ii)
 150        continue
 160     continue
 165     continue
         do 170 i = 1, ni
            dk(jx+i,jx+i) = four * dk(jx+i,jx+i)
 170     continue
         jx = jx + ni
 180  continue
c
      return
      end
