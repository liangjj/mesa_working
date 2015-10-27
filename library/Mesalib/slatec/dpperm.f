*deck dpperm
      subroutine dpperm (dx, n, iperm, ier)
c***begin prologue  dpperm
c***purpose  rearrange a given array according to a prescribed
c            permutation vector.
c***library   slatec
c***category  n8
c***type      double precision (spperm-s, dpperm-d, ipperm-i, hpperm-h)
c***keywords  permutation, rearrangement
c***author  mcclain, m. a., (nist)
c           rhoads, g. s., (nbs)
c***description
c
c         dpperm rearranges the data vector dx according to the
c         permutation iperm: dx(i) <--- dx(iperm(i)).  iperm could come
c         from one of the sorting routines ipsort, spsort, dpsort or
c         hpsort.
c
c     description of parameters
c         dx - input/output -- double precision array of values to be
c                   rearranged.
c         n - input -- number of values in double precision array dx.
c         iperm - input -- permutation vector.
c         ier - output -- error indicator:
c             =  0  if no error,
c             =  1  if n is zero or negative,
c             =  2  if iperm is not a valid permutation.
c
c***references  (none)
c***routines called  xermsg
c***revision history  (yymmdd)
c   901004  date written
c   920507  modified by m. mcclain to revise prologue text.
c***end prologue  dpperm
      integer n, iperm(*), i, ier, indx, indx0, istrt
      double precision dx(*), dtemp
c***first executable statement  dpperm
      ier=0
      if(n.lt.1)then
         ier=1
         call xermsg ('slatec', 'dpperm',
     +    'the number of values to be rearranged, n, is not positive.',
     +    ier, 1)
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
         ier=2
         call xermsg ('slatec', 'dpperm',
     +    'the permutation vector, iperm, is not valid.', ier, 1)
         return
  100 continue
c
c     rearrange the values of dx
c
c     use the iperm vector as a flag.
c     if iperm(i) > 0, then the i-th value is in correct location
c
      do 330 istrt = 1 , n
         if (iperm(istrt) .gt. 0) goto 330
         indx = istrt
         indx0 = indx
         dtemp = dx(istrt)
  320    continue
         if (iperm(indx) .ge. 0) goto 325
            dx(indx) = dx(-iperm(indx))
            indx0 = indx
            iperm(indx) = -iperm(indx)
            indx = iperm(indx)
            goto 320
  325    continue
         dx(indx0) = dtemp
  330 continue
c
      return
      end
