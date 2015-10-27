*deck qs2i1r
      subroutine qs2i1r (ia, ja, a, n, kflag)
c***begin prologue  qs2i1r
c***subsidiary
c***purpose  sort an integer array, moving an integer and real array.
c            this routine sorts the integer array ia and makes the same
c            interchanges in the integer array ja and the real array a.
c            the array ia may be sorted in increasing order or decreas-
c            ing order.  a slightly modified quicksort algorithm is
c            used.
c***library   slatec (slap)
c***category  n6a2a
c***type      single precision (qs2i1r-s, qs2i1d-d)
c***keywords  singleton quicksort, slap, sort, sorting
c***author  jones, r. e., (snla)
c           kahaner, d. k., (nbs)
c           seager, m. k., (llnl) seager@llnl.gov
c           wisniewski, j. a., (snla)
c***description
c     written by rondall e jones
c     modified by john a. wisniewski to use the singleton quicksort
c     algorithm. date 18 november 1976.
c
c     further modified by david k. kahaner
c     national bureau of standards
c     august, 1981
c
c     even further modification made to bring the code up to the
c     fortran 77 level and make it more readable and to carry
c     along one integer array and one real array during the sort by
c     mark k. seager
c     lawrence livermore national laboratory
c     november, 1987
c     this routine was adapted from the isort routine.
c
c     abstract
c         this routine sorts an integer array ia and makes the same
c         interchanges in the integer array ja and the real array a.
c         the array ia may be sorted in increasing order or decreasing
c         order.  a slightly modified quicksort algorithm is used.
c
c     description of parameters
c        ia - integer array of values to be sorted.
c        ja - integer array to be carried along.
c         a - real array to be carried along.
c         n - number of values in integer array ia to be sorted.
c     kflag - control parameter
c           = 1 means sort ia in increasing order.
c           =-1 means sort ia in decreasing order.
c
c***see also  ss2y
c***references  r. c. singleton, algorithm 347, an efficient algorithm
c                 for sorting with minimal storage, communications acm
c                 12:3 (1969), pp.185-7.
c***routines called  xermsg
c***revision history  (yymmdd)
c   761118  date written
c   890125  previous revision date
c   890915  made changes requested at july 1989 cml meeting.  (mks)
c   890922  numerous changes to prologue to make closer to slatec
c           standard.  (fnf)
c   890929  numerous changes to reduce sp/dp differences.  (fnf)
c   900805  changed xerror calls to calls to xermsg.  (rwc)
c   910411  prologue converted to version 4.0 format.  (bab)
c   910506  made subsidiary to ss2y and corrected reference.  (fnf)
c   920511  added complete declaration section.  (wrb)
c   920929  corrected format of reference.  (fnf)
c   921012  added e0's to f.p. constants.  (fnf)
c***end prologue  qs2i1r
cvd$r novector
cvd$r noconcur
c     .. scalar arguments ..
      integer kflag, n
c     .. array arguments ..
      real a(n)
      integer ia(n), ja(n)
c     .. local scalars ..
      real r, ta, tta
      integer i, iit, ij, it, j, jjt, jt, k, kk, l, m, nn
c     .. local arrays ..
      integer il(21), iu(21)
c     .. external subroutines ..
      external xermsg
c     .. intrinsic functions ..
      intrinsic abs, int
c***first executable statement  qs2i1r
      nn = n
      if (nn.lt.1) then
         call xermsg ('slatec', 'qs2i1r',
     $      'the number of values to be sorted was not positive.', 1, 1)
         return
      endif
      if( n.eq.1 ) return
      kk = abs(kflag)
      if ( kk.ne.1 ) then
         call xermsg ('slatec', 'qs2i1r',
     $      'the sort control parameter, k, was not 1 or -1.', 2, 1)
         return
      endif
c
c     alter array ia to get decreasing order if needed.
c
      if( kflag.lt.1 ) then
         do 20 i=1,nn
            ia(i) = -ia(i)
 20      continue
      endif
c
c     sort ia and carry ja and a along.
c     and now...just a little black magic...
      m = 1
      i = 1
      j = nn
      r = .375e0
 210  if( r.le.0.5898437e0 ) then
         r = r + 3.90625e-2
      else
         r = r-.21875e0
      endif
 225  k = i
c
c     select a central element of the array and save it in location
c     it, jt, at.
c
      ij = i + int ((j-i)*r)
      it = ia(ij)
      jt = ja(ij)
      ta = a(ij)
c
c     if first element of array is greater than it, interchange with it.
c
      if( ia(i).gt.it ) then
         ia(ij) = ia(i)
         ia(i)  = it
         it     = ia(ij)
         ja(ij) = ja(i)
         ja(i)  = jt
         jt     = ja(ij)
         a(ij)  = a(i)
         a(i)   = ta
         ta     = a(ij)
      endif
      l=j
c
c     if last element of array is less than it, swap with it.
c
      if( ia(j).lt.it ) then
         ia(ij) = ia(j)
         ia(j)  = it
         it     = ia(ij)
         ja(ij) = ja(j)
         ja(j)  = jt
         jt     = ja(ij)
         a(ij)  = a(j)
         a(j)   = ta
         ta     = a(ij)
c
c     if first element of array is greater than it, swap with it.
c
         if ( ia(i).gt.it ) then
            ia(ij) = ia(i)
            ia(i)  = it
            it     = ia(ij)
            ja(ij) = ja(i)
            ja(i)  = jt
            jt     = ja(ij)
            a(ij)  = a(i)
            a(i)   = ta
            ta     = a(ij)
         endif
      endif
c
c     find an element in the second half of the array which is
c     smaller than it.
c
  240 l=l-1
      if( ia(l).gt.it ) go to 240
c
c     find an element in the first half of the array which is
c     greater than it.
c
  245 k=k+1
      if( ia(k).lt.it ) go to 245
c
c     interchange these elements.
c
      if( k.le.l ) then
         iit   = ia(l)
         ia(l) = ia(k)
         ia(k) = iit
         jjt   = ja(l)
         ja(l) = ja(k)
         ja(k) = jjt
         tta   = a(l)
         a(l)  = a(k)
         a(k)  = tta
         goto 240
      endif
c
c     save upper and lower subscripts of the array yet to be sorted.
c
      if( l-i.gt.j-k ) then
         il(m) = i
         iu(m) = l
         i = k
         m = m+1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m+1
      endif
      go to 260
c
c     begin again on another portion of the unsorted array.
c
  255 m = m-1
      if( m.eq.0 ) go to 300
      i = il(m)
      j = iu(m)
  260 if( j-i.ge.1 ) go to 225
      if( i.eq.j ) go to 255
      if( i.eq.1 ) go to 210
      i = i-1
  265 i = i+1
      if( i.eq.j ) go to 255
      it = ia(i+1)
      jt = ja(i+1)
      ta =  a(i+1)
      if( ia(i).le.it ) go to 265
      k=i
  270 ia(k+1) = ia(k)
      ja(k+1) = ja(k)
      a(k+1)  =  a(k)
      k = k-1
      if( it.lt.ia(k) ) go to 270
      ia(k+1) = it
      ja(k+1) = jt
      a(k+1)  = ta
      go to 265
c
c     clean up, if necessary.
c
  300 if( kflag.lt.1 ) then
         do 310 i=1,nn
            ia(i) = -ia(i)
 310     continue
      endif
      return
c------------- last line of qs2i1r follows ----------------------------
      end
