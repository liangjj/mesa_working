*deck @(#)twoout.f	5.1  11/6/94
      subroutine twoout(lab2,valtwo,nbf,m)
c
c  take a one-dimensional array of two-electron integrals stored in
c  canonical order and write them out with labels.
c  zero integrals are not listed.
c
c  arguments:
c            lab2   integral label
c           valtwo  the array to be written out
c            nbf    the number of basis functions
c             m     the number of integrals written per line (le 4)
c
      implicit real*8 (a-h,o-z)
      character*8 lab2
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      dimension valtwo(*), ii(4), jj(4), kk(4), ll(4), a(4)
c
      data xlimit/1.0d-8/
c
      if( (m .lt. 1) .or. (m .gt. 4) )  m = 4
c
      write (iw,10) lab2
   10 format(///5x,'list of the two-electron integrals'////
     1          5x,'integral label - ',a8/)
      it = 0
      ijkl = 0
c
      do 60 i=1,nbf
      do 60 j=1,i
      do 60 k=1,i
      if( k .eq. i ) then
        lmax = j
      else
        lmax = k
      endif
      do 50 l=1,lmax
c
      ijkl = ijkl + 1
      if( abs(valtwo(ijkl)) .lt. xlimit ) go to 50
c
      if( it .eq. m ) then
        write (iw,70) (ii(mm), jj(mm), kk(mm), ll(mm), a(mm), mm=1,m)
        it = 0
      endif
c
      it = it + 1
      ii(it) = i
      jj(it) = j
      kk(it) = k
      ll(it) = l
      a(it)  = valtwo(ijkl)
   50 continue
   60 continue
      write (iw,70) (ii(mm), jj(mm), kk(mm), ll(mm), a(mm), mm=1,it)
c
   70 format(4(4x,4i3,f15.8))
c
      write (iw,80) xlimit
   80 format(/5x,'integrals not listed if less than',1p,d10.1)
c
      return
      end
