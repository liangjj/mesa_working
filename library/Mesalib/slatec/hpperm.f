*deck hpperm
      subroutine hpperm (hx, n, iperm, work, ier)
c***begin prologue  hpperm
c***purpose  rearrange a given array according to a prescribed
c            permutation vector.
c***library   slatec
c***category  n8
c***type      character (spperm-s, dpperm-d, ipperm-i, hpperm-h)
c***keywords  application of permutation to data vector
c***author  mcclain, m. a., (nist)
c           rhoads, g. s., (nbs)
c***description
c
c         hpperm rearranges the data vector hx according to the
c         permutation iperm: hx(i) <--- hx(iperm(i)).  iperm could come
c         from one of the sorting routines ipsort, spsort, dpsort or
c         hpsort.
c
c     description of parameters
c         hx - input/output -- character array of values to be
c                 rearranged.
c         n - input -- number of values in character array hx.
c         iperm - input -- permutation vector.
c         work - character variable which must have a length
c                   specification at least as great as that of hx.
c         ier - output -- error indicator:
c             =  0  if no error,
c             =  1  if n is zero or negative,
c             =  2  if work array is not long enough,
c             =  3  if iperm is not a valid permutation.
c
c***references  (none)
c***routines called  xermsg
c***revision history  (yymmdd)
c   901004  date written
c   920507  modified by m. mcclain to revise prologue text and to add
c           check for length of work array.
c***end prologue  hpperm
      integer n, iperm(*), i, ier, indx, indx0, istrt
      character*(*) hx(*), work
c***first executable statement  hpperm
      ier=0
      if(n.lt.1)then
         ier=1
         call xermsg ('slatec', 'hpperm',
     +    'the number of values to be rearranged, n, is not positive.',
     +    ier, 1)
         return
      endif
      if(len(work).lt.len(hx(1)))then
         ier=2
         call xermsg ('slatec', 'hpperm',
     +    'the length of the work variable, work, is too short.',ier,1)
         return
      endif
c
c     check whether iperm is a valid permutation
c
      do 100 i=1,n
         indx=abs(iperm(i))
         if((indx.ge.1).and.(indx.le.n))then
            if(iperm(indx).gt.0)then
               iperm(indx)=-iperm(indx)
               goto 100
            endif
         endif
         ier=3
         call xermsg ('slatec', 'hpperm',
     +    'the permutation vector, iperm, is not valid.', ier, 1)
         return
  100 continue
c
c     rearrange the values of hx
c
c     use the iperm vector as a flag.
c     if iperm(i) > 0, then the i-th value is in correct location
c
      do 330 istrt = 1 , n
         if (iperm(istrt) .gt. 0) goto 330
         indx = istrt
         indx0 = indx
         work = hx(istrt)
  320    continue
         if (iperm(indx) .ge. 0) goto 325
            hx(indx) = hx(-iperm(indx))
            indx0 = indx
            iperm(indx) = -iperm(indx)
            indx = iperm(indx)
            goto 320
  325    continue
         hx(indx0) = work
  330 continue
c
      return
      end
