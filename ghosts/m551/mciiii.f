*deck %W%  %G%
      subroutine mciiii(dj,dkk,dk,dkt,ni,nii,den,ipair)
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
      implicit real*8(a-h,o-z)
cc
cmp   extended dummy dj,dk,dkt,dkk,den,ipair
cc
      dimension dj(nii,nii),dk(nii,nii),dkt(nii,nii),dkk(nii,2)
      dimension den(2),ipair(2)
      common /number/ zero,pt5,one,two,four,eight
      common /io/ inp,iout
c
c
c
      if (ni .eq. 1) go to 1000
c
      ix = 0
      do 50 i = 1, nii
         do 40 j = 1, i
            ix = ix + 1
            dj(i,j) = den(ix)
            dj(j,i) = den(ix)
            dkt(i,j) = zero
            dkt(j,i) = zero
 40      continue
 50   continue
c
c       write(iout,*)' dj nii ',nii
c       call vecout(dj,nii,nii)
c
c     build exchange block
c
cccc
      do 100 i = 1, ni
cccc
         ipi = ipair(i)
ccc
         do 200 k = 1, i
ccc
            ipk = ipair(k)
            ipik = ipi + k
cc
            do 300 j = 1, ni
cc
               ipj = ipair(j)
c
               if (j .gt. i) go to 315
               ipij = ipi + j
               if (k .gt. j) go to 310
               ipjk = ipj + k
               go to 320
 310           continue
               ipjk = ipk + j
               go to 320
c
 315           continue
               ipij = ipj + i
               ipjk = ipj + k
 320           continue
c
c     build exchange blocks for the iiii symmetry case
c
               do 400 l = 1, j
                  ipl = ipair(l)
                  if (l .gt. i) go to 410
                  ipil = ipi + l
                  go to 415
 410              continue
                  ipil = ipl + i
 415              continue
c
                  if (l .gt. j) go to 420
                  iplj = ipj + l
                  go to 425
 420              continue
                  iplj = ipl + j
 425              continue
                  if (l .gt. k) go to 430
                  ipkl = ipk + l
                  go to 435
 430              continue
                  ipkl = ipl + k
 435              continue
c
                  dk(ipik,iplj) = dj(ipij,ipkl)
 400           continue
 300        continue
 200     continue
 100  continue
c
c
      do 540 i=2,ni
         i1=i-1
         ipi=ipair(i)
         do 530 j=1,i1
            ipj=ipair(j)
            ij=ipi+j
            do 520 k=2,ni
               k1=k-1
               ipk=ipair(k)
               if(k.lt.j) go to 501
               jk=ipk+j
               go to 502
 501           continue
               jk=ipj+k
 502           continue
               do 510 l=1,k1
                  kl=ipk+l
                  if(l.gt.i) go to 503
                  il=ipi+l
                  go to 504
 503              continue
                  il=ipair(l)+i
 504              continue
                  dkt(ij,kl)=dj(il,jk)
 510           continue
 520        continue
 530     continue
 540  continue
c
c
c     scale the coulomb term
c
      do 550 ix = 1, nii
         dj(ix,ix) = dj(ix,ix) * two
         dk(ix,ix) = dk(ix,ix) * four
 550  continue
c
c---------------------------------------------------c
c    jiki  exchange terms
c      j>k with  k>i and k=i
c      2111 3121 3222 4111 4121 4222 4131 4232 etc.
c---------------------------------------------------c
      do 730 i = 2, ni
         ipi = ipair(i)
         i1 = i - 1
         do 720 j = 1, i1
            ipj = ipair(j)
            do 710 k = 1, j
               dk(ipi+k,ipj+k) = two * dk(ipi+k,ipj+k)
               dk(ipj+k,ipi+k) = dk(ipi+k,ipj+k)
 710        continue
 720     continue
 730  continue
c---------------------------------------------------c
c    ijik exchange terms
c      i>j and i=j,  j>k
c      3231 3331 3332  etc.
c---------------------------------------------------c
      do 830 i = 2, ni
         ipi = ipair(i)
         do 820 j = 2, i
            j1 = j - 1
            do 810 k = 1, j1
               dk(ipi+j,ipi+k) = two * dk(ipi+j,ipi+k)
               dk(ipi+k,ipi+j) = dk(ipi+j,ipi+k)
 810        continue
 820     continue
 830  continue
c-------------------c
c     dkt           c
c-------------------c
c
c--------------------------------c
c   ijjl transpose exchange terms
c        i>j ,  j >or= l
c     2111 3221 3222 etc.
c  *****  and  *****
c   ijij transpose exchange terms
c        i>j
c--------------------------------c
c
      do 940 i = 2, ni
         i1 = i - 1
         ipi = ipair(i)
         do 930 j = 1, i1
            ipij = ipi + j
            ipj = ipair(j)
            do 920 l = 1, j
               dkt(ipj+l,ipij) = two * dkt(ipj+l,ipij)
               dkt(ipij,ipj+l) = dkt(ipj+l,ipij)
 920        continue
            dkt(ipij,ipij) = two * dkt(ipij,ipij)
 930     continue
 940  continue
c-------------------------------------c
c  move  dk and dkt  into  dkk        c
c-------------------------------------c
c
      do 1040 i=1,nii
c     jx=0
         ix=0
         do 1030 j=1,ni
            kx=0
            do 1020 k=1,j
               ix=ix+1
c     dkk(i,jx+k)=dk(i,ix)
               dkk(i,kx+j)=dk(i,ix)
               kx=kx+ni
 1020       continue
c     jx=jx+ni
 1030    continue
 1040 continue
c
      do 1140 i=1,nii
         ix=1
         jx=ni
         do 1130 j=2,ni
c     kx=0
            j1=j-1
            do 1120 k=1,j1
               ix=ix+1
c     dkk(i,kx+j)=dkt(i,ix)
               dkk(i,jx+k)=dkt(i,ix)
 1120       continue
            ix=ix+1
            jx=jx+ni
 1130    continue
 1140 continue
c
      return
c
 1000 continue
      dj(1,1) = two * den(1)
      dkk(1,1) = four * den(1)
      dkt(1,1) = zero
      return
      end
