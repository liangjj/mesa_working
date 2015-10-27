*deck @(#)oneout.f	5.1  11/6/94
      subroutine oneout(lab1,valone,nbf,m,n)
c
c  take a one-dimensional array of one-electron integrals stored in
c  canonical order and write them out with labels.
c  zero integrals are not listed.
c
c  arguments:
c            lab1   integral label
c           valone  the array to be written out
c            nbf    the number of basis functions
c             m     the number of integrals written per line (le 5)
c             n     the integer format for writing out the array
c
      implicit real*8 (a-h,o-z)
      character*(*) lab1
      character*20 fmt1
      common/tapes/iw,iunt1a,iunt1b,iunt2a,iunts1,iunts2
      dimension valone(*), ii(5), jj(5), s(5)
c
      data xlimit/1.0d-8/
c
      if( (m .lt. 1) .or. (m .gt. 5) ) m = 5
      if( n .gt. 3 ) m = min(m,4)
c
      write (fmt1,10) n
   10 format('(4x,5(2i',i1,',f15.8,4x))')
c
      write (iw,30) lab1
   30 format('1',4x,'list of the one-electron integrals'////
     1           5x,'integral label - ',a32/)
      it = 0
      ij = 0
c
      do 60 i=1,nbf
      do 50 j=1,i
c
      ij = ij + 1
      if( abs(valone(ij)) .lt. xlimit ) go to 50
c
      if( it .eq. m ) then
        write (iw,fmt1) (ii(mm), jj(mm), s(mm), mm=1,m)
        it = 0
      endif
c
      it = it + 1
      ii(it) = i
      jj(it) = j
      s(it)  = valone(ij)
   50 continue
   60 continue
      write (iw,fmt1) (ii(mm), jj(mm), s(mm), mm=1,it)
c
      write (iw,70) xlimit
   70 format(/5x,'integrals not listed if less than',1p,d10.1)
c
      return
      end
