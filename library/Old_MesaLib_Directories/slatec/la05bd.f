*deck la05bd
      subroutine la05bd (a, ind, ia, n, ip, iw, w, g, b, trans)
c***begin prologue  la05bd
c***subsidiary
c***purpose  subsidiary to dsplp
c***library   slatec
c***type      double precision (la05bs-s, la05bd-d)
c***author  (unknown)
c***description
c
c     this subprogram is a slight modification of a subprogram
c     from the c. 1979 aere harwell library.  the name of the
c     corresponding harwell code can be obtained by deleting
c     the final letter =d= in the names used here.
c     revised sep. 13, 1979.
c
c     royalties have been paid to aere-uk for use of their codes
c     in the package given here.  any primary usage of the harwell
c     subroutines requires a royalty agreement and payment between
c     the user and aere-uk.  any usage of the sandia written codes
c     dsplp( ) (which uses the harwell subroutines) is permitted.
c
c ip(i,1),ip(i,2) point to start of row/column i of u.
c iw(i,1),iw(i,2) are lengths of row/col i of u.
c iw(.,3),iw(.,4) hold row/col numbers in pivotal order.
c
c***see also  dsplp
c***routines called  xermsg, xsetun
c***common blocks    la05dd
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900402  added type section.  (wrb)
c   920410  corrected second dimension on iw declaration.  (wrb)
c***end prologue  la05bd
      double precision a(*), b(*), am, w(*), g, small
      logical trans
      integer ind(ia,2), iw(n,8)
      integer ip(n,2)
      common /la05dd/ small, lp, lenl, lenu, ncp, lrow, lcol
c***first executable statement  la05bd
      if (g.lt.0.d0) go to 130
      kll = ia - lenl + 1
      if (trans) go to 80
c
c     multiply vector by inverse of l
      if (lenl.le.0) go to 20
      l1 = ia + 1
      do 10 kk=1,lenl
         k = l1 - kk
         i = ind(k,1)
         if (b(i).eq.0.d0) go to 10
         j = ind(k,2)
         b(j) = b(j) + a(k)*b(i)
   10 continue
   20 do 30 i=1,n
         w(i) = b(i)
         b(i) = 0.d0
   30 continue
c
c     multiply vector by inverse of u
      n1 = n + 1
      do 70 ii=1,n
         i = n1 - ii
         i = iw(i,3)
         am = w(i)
         kp = ip(i,1)
         if (kp.gt.0) go to 50
         kp = -kp
         ip(i,1) = kp
         nz = iw(i,1)
         kl = kp - 1 + nz
         k2 = kp + 1
         do 40 k=k2,kl
            j = ind(k,2)
            am = am - a(k)*b(j)
   40    continue
   50    if (am.eq.0.) go to 70
         j = ind(kp,2)
         b(j) = am/a(kp)
         kpc = ip(j,2)
         kl = iw(j,2) + kpc - 1
         if (kl.eq.kpc) go to 70
         k2 = kpc + 1
         do 60 k=k2,kl
            i = ind(k,1)
            ip(i,1) = -abs(ip(i,1))
   60    continue
   70 continue
      go to 140
c
c     multiply vector by inverse of transpose of u
   80 do 90 i=1,n
         w(i) = b(i)
         b(i) = 0.d0
   90 continue
      do 110 ii=1,n
         i = iw(ii,4)
         am = w(i)
         if (am.eq.0.d0) go to 110
         j = iw(ii,3)
         kp = ip(j,1)
         am = am/a(kp)
         b(j) = am
         kl = iw(j,1) + kp - 1
         if (kp.eq.kl) go to 110
         k2 = kp + 1
         do 100 k=k2,kl
            i = ind(k,2)
            w(i) = w(i) - am*a(k)
  100    continue
  110 continue
c
c     multiply vector by inverse of transpose of l
      if (kll.gt.ia) return
      do 120 k=kll,ia
         j = ind(k,2)
         if (b(j).eq.0.d0) go to 120
         i = ind(k,1)
         b(i) = b(i) + a(k)*b(j)
  120 continue
      go to 140
c
  130 call xsetun(lp)
      if (lp .gt. 0) call xermsg ('slatec', 'la05bd',
     +   'earlier entry gave error return.', -8, 2)
  140 return
      end
