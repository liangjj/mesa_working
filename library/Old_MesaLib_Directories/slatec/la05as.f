*deck la05as
      subroutine la05as (a, ind, nz, ia, n, ip, iw, w, g, u)
c***begin prologue  la05as
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (la05as-s, la05ad-d)
c***author  (unknown)
c***description
c
c     this subprogram is a slight modification of a subprogram
c     from the c. 1979 aere harwell library.  the name of the
c     corresponding harwell code can be obtained by deleting
c     the final letter =s= in the names used here.
c     revisions made by r j hanson, snla, august, 1979.
c     revised sep. 13, 1979.
c
c     royalties have been paid to aere-uk for use of their codes
c     in the package given here.  any primary usage of the harwell
c     subroutines requires a royalty agreement and payment between
c     the user and aere-uk.  any usage of the sandia written codes
c     splp( ) (which uses the harwell subroutines) is permitted.
c
c ip(i,1),ip(i,2) point to the start of row/col i.
c iw(i,1),iw(i,2) hold the number of non-zeros in row/col i.
c during the main body of this subroutine the vectors iw(.,3),iw(.,5),
c     iw(.,7) are used to hold doubly linked lists of rows that have
c     not been pivotal and have equal numbers of non-zeros.
c iw(.,4),iw(.,6),iw(.,8) hold similar lists for the columns.
c iw(i,3),iw(i,4) hold first row/column to have i non-zeros
c     or zero if there are none.
c iw(i,5), iw(i,6) hold row/col number of row/col prior to row/col i
c     in its list, or zero if none.
c iw(i,7), iw(i,8) hold row/col number of row/col after row/col i
c     in its list, or zero if none.
c for rows/cols that have been pivotal iw(i,5),iw(i,6) hold negation of
c     position of row/col i in the pivotal ordering.
c
c***see also  splp
c***routines called  la05es, mc20as, r1mach, xermsg, xsetun
c***common blocks    la05ds
c***revision history  (yymmdd)
c   811215  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890605  corrected references to xerrwv.  (wrb)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900402  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c***end prologue  la05as
      integer ip(n,2)
      integer ind(ia,2), iw(n,8)
      real a(*), amax, au, am, g, u, small, w(*)
      logical first
      character*8 xern0, xern1, xern2
c
      common /la05ds/ small, lp, lenl, lenu, ncp, lrow, lcol
c eps is the relative accuracy of floating-point computation
      save eps, first
      data first /.true./
c***first executable statement  la05as
      if (first) then
         eps = 2.0e0 * r1mach(4)
      endif
      first = .false.
c
c     set the output unit number for the error processor.
c     the usage of this error processor is documented in the
c     sandia labs. tech. rept. sand78-1189, by r e jones.
      call xsetun(lp)
      if (u.gt.1.0e0) u = 1.0e0
      if (u.lt.eps) u = eps
      if (n.lt.1) go to 670
      g = 0.
      do 50 i=1,n
         w(i) = 0.
         do 40 j=1,5
            iw(i,j) = 0
   40    continue
   50 continue
c
c flush out small entries, count elements in rows and columns
      l = 1
      lenu = nz
      do 80 idummy=1,nz
         if (l.gt.lenu) go to 90
         do 60 k=l,lenu
            if (abs(a(k)).le.small) go to 70
            i = ind(k,1)
            j = ind(k,2)
            g = max(abs(a(k)),g)
            if (i.lt.1 .or. i.gt.n) go to 680
            if (j.lt.1 .or. j.gt.n) go to 680
            iw(i,1) = iw(i,1) + 1
            iw(j,2) = iw(j,2) + 1
   60    continue
         go to 90
   70    l = k
         a(l) = a(lenu)
         ind(l,1) = ind(lenu,1)
         ind(l,2) = ind(lenu,2)
         lenu = lenu - 1
   80 continue
c
   90 lenl = 0
      lrow = lenu
      lcol = lrow
c mcp is the maximum number of compresses permitted before an
c     error return results.
      mcp = max(n/10,20)
      ncp = 0
c check for null row or column and initialize ip(i,2) to point
c     just beyond where the last component of column i of a will
c     be stored.
      k = 1
      do 110 ir=1,n
         k = k + iw(ir,2)
         ip(ir,2) = k
         do 100 l=1,2
            if (iw(ir,l).le.0) go to 700
  100    continue
  110 continue
c reorder by rows
c check for double entries while using the newly constructed
c     row file to construct the column file. note that by putting
c    the entries in backwards and decreasing ip(j,2) each time it
c     is used we automatically leave it pointing to the first element.
      call mc20as(n, lenu, a, ind(1,2), ip, ind(1,1), 0)
      kl = lenu
      do 130 ii=1,n
         ir = n + 1 - ii
         kp = ip(ir,1)
         do 120 k=kp,kl
            j = ind(k,2)
            if (iw(j,5).eq.ir) go to 660
            iw(j,5) = ir
            kr = ip(j,2) - 1
            ip(j,2) = kr
            ind(kr,1) = ir
  120    continue
         kl = kp - 1
  130 continue
c
c set up linked lists of rows and cols with equal numbers of non-zeros.
      do 150 l=1,2
         do 140 i=1,n
            nz = iw(i,l)
            in = iw(nz,l+2)
            iw(nz,l+2) = i
            iw(i,l+6) = in
            iw(i,l+4) = 0
            if (in.ne.0) iw(in,l+4) = i
  140    continue
  150 continue
c
c
c start of main elimination loop.
      do 590 ipv=1,n
c find pivot. jcost is markowitz cost of cheapest pivot found so far,
c     which is in row ipp and column jp.
         jcost = n*n
c loop on length of column to be searched
         do 240 nz=1,n
            if (jcost.le.(nz-1)**2) go to 250
            j = iw(nz,4)
c search columns with nz non-zeros.
            do 190 idummy=1,n
               if (j.le.0) go to 200
               kp = ip(j,2)
               kl = kp + iw(j,2) - 1
               do 180 k=kp,kl
                  i = ind(k,1)
                  kcost = (nz-1)*(iw(i,1)-1)
                  if (kcost.ge.jcost) go to 180
                  if (nz.eq.1) go to 170
c find largest element in row of potential pivot.
                  amax = 0.
                  k1 = ip(i,1)
                  k2 = iw(i,1) + k1 - 1
                  do 160 kk=k1,k2
                     amax = max(amax,abs(a(kk)))
                     if (ind(kk,2).eq.j) kj = kk
  160             continue
c perform stability test.
                  if (abs(a(kj)).lt.amax*u) go to 180
  170             jcost = kcost
                  ipp = i
                  jp = j
                  if (jcost.le.(nz-1)**2) go to 250
  180          continue
               j = iw(j,8)
  190       continue
c search rows with nz non-zeros.
  200       i = iw(nz,3)
            do 230 idummy=1,n
               if (i.le.0) go to 240
               amax = 0.
               kp = ip(i,1)
               kl = kp + iw(i,1) - 1
c find largest element in the row
               do 210 k=kp,kl
                  amax = max(abs(a(k)),amax)
  210          continue
               au = amax*u
               do 220 k=kp,kl
c perform stability test.
                  if (abs(a(k)).lt.au) go to 220
                  j = ind(k,2)
                  kcost = (nz-1)*(iw(j,2)-1)
                  if (kcost.ge.jcost) go to 220
                  jcost = kcost
                  ipp = i
                  jp = j
                  if (jcost.le.(nz-1)**2) go to 250
  220          continue
               i = iw(i,7)
  230       continue
  240    continue
c
c pivot found.
c remove rows and columns involved in elimination from ordering vectors.
  250    kp = ip(jp,2)
         kl = iw(jp,2) + kp - 1
         do 290 l=1,2
            do 280 k=kp,kl
               i = ind(k,l)
               il = iw(i,l+4)
               in = iw(i,l+6)
               if (il.eq.0) go to 260
               iw(il,l+6) = in
               go to 270
  260          nz = iw(i,l)
               iw(nz,l+2) = in
  270          if (in.gt.0) iw(in,l+4) = il
  280       continue
            kp = ip(ipp,1)
            kl = kp + iw(ipp,1) - 1
  290    continue
c store pivot
         iw(ipp,5) = -ipv
         iw(jp,6) = -ipv
c eliminate pivotal row from column file and find pivot in row file.
         do 320 k=kp,kl
            j = ind(k,2)
            kpc = ip(j,2)
            iw(j,2) = iw(j,2) - 1
            klc = kpc + iw(j,2)
            do 300 kc=kpc,klc
               if (ipp.eq.ind(kc,1)) go to 310
  300       continue
  310       ind(kc,1) = ind(klc,1)
            ind(klc,1) = 0
            if (j.eq.jp) kr = k
  320    continue
c bring pivot to front of pivotal row.
         au = a(kr)
         a(kr) = a(kp)
         a(kp) = au
         ind(kr,2) = ind(kp,2)
         ind(kp,2) = jp
c
c perform elimination itself, looping on non-zeros in pivot column.
         nzc = iw(jp,2)
         if (nzc.eq.0) go to 550
         do 540 nc=1,nzc
            kc = ip(jp,2) + nc - 1
            ir = ind(kc,1)
c search non-pivot row for element to be eliminated.
            kr = ip(ir,1)
            krl = kr + iw(ir,1) - 1
            do 330 knp=kr,krl
               if (jp.eq.ind(knp,2)) go to 340
  330       continue
c bring element to be eliminated to front of its row.
  340       am = a(knp)
            a(knp) = a(kr)
            a(kr) = am
            ind(knp,2) = ind(kr,2)
            ind(kr,2) = jp
            am = -a(kr)/a(kp)
c compress row file unless it is certain that there is room for new row.
            if (lrow+iw(ir,1)+iw(ipp,1)+lenl.le.ia) go to 350
            if (ncp.ge.mcp .or. lenu+iw(ir,1)+iw(ipp,1)+lenl.gt.ia) go
     *       to 710
            call la05es(a, ind(1,2), ip, n, iw, ia, .true.)
            kp = ip(ipp,1)
            kr = ip(ir,1)
  350       krl = kr + iw(ir,1) - 1
            kq = kp + 1
            kpl = kp + iw(ipp,1) - 1
c place pivot row (excluding pivot itself) in w.
            if (kq.gt.kpl) go to 370
            do 360 k=kq,kpl
               j = ind(k,2)
               w(j) = a(k)
  360       continue
  370       ip(ir,1) = lrow + 1
c
c transfer modified elements.
            ind(kr,2) = 0
            kr = kr + 1
            if (kr.gt.krl) go to 430
            do 420 ks=kr,krl
               j = ind(ks,2)
               au = a(ks) + am*w(j)
               ind(ks,2) = 0
c if element is very small remove it from u.
               if (abs(au).le.small) go to 380
               g = max(g,abs(au))
               lrow = lrow + 1
               a(lrow) = au
               ind(lrow,2) = j
               go to 410
  380          lenu = lenu - 1
c remove element from col file.
               k = ip(j,2)
               kl = k + iw(j,2) - 1
               iw(j,2) = kl - k
               do 390 kk=k,kl
                  if (ind(kk,1).eq.ir) go to 400
  390          continue
  400          ind(kk,1) = ind(kl,1)
               ind(kl,1) = 0
  410          w(j) = 0.
  420       continue
c
c scan pivot row for fills.
  430       if (kq.gt.kpl) go to 520
            do 510 ks=kq,kpl
               j = ind(ks,2)
               au = am*w(j)
               if (abs(au).le.small) go to 500
               lrow = lrow + 1
               a(lrow) = au
               ind(lrow,2) = j
               lenu = lenu + 1
c
c create fill in column file.
               nz = iw(j,2)
               k = ip(j,2)
               kl = k + nz - 1
               if (nz .eq. 0) go to 460
c if possible place new element at end of present entry.
               if (kl.ne.lcol) go to 440
               if (lcol+lenl.ge.ia) go to 460
               lcol = lcol + 1
               go to 450
  440          if (ind(kl+1,1).ne.0) go to 460
  450          ind(kl+1,1) = ir
               go to 490
c new entry has to be created.
  460          if (lcol+lenl+nz+1.lt.ia) go to 470
c compress column file if there is not room for new entry.
               if (ncp.ge.mcp .or. lenu+lenl+nz+1.ge.ia) go to 710
               call la05es(a, ind, ip(1,2), n, iw(1,2), ia, .false.)
               k = ip(j,2)
               kl = k + nz - 1
c transfer old entry into new.
  470          ip(j,2) = lcol + 1
               if (kl .lt. k) go to 485
               do 480 kk=k,kl
                  lcol = lcol + 1
                  ind(lcol,1) = ind(kk,1)
                  ind(kk,1) = 0
  480          continue
  485          continue
c add new element.
               lcol = lcol + 1
               ind(lcol,1) = ir
  490          g = max(g,abs(au))
               iw(j,2) = nz + 1
  500          w(j) = 0.
  510       continue
  520       iw(ir,1) = lrow + 1 - ip(ir,1)
c
c store multiplier
            if (lenl+lcol+1.le.ia) go to 530
c compress col file if necessary.
            if (ncp.ge.mcp) go to 710
            call la05es(a, ind, ip(1,2), n, iw(1,2), ia, .false.)
  530       k = ia - lenl
            lenl = lenl + 1
            a(k) = am
            ind(k,1) = ipp
            ind(k,2) = ir
            lenu = lenu - 1
  540    continue
c
c insert rows and columns involved in elimination in linked lists
c     of equal numbers of non-zeros.
  550    k1 = ip(jp,2)
         k2 = iw(jp,2) + k1 - 1
         iw(jp,2) = 0
         do 580 l=1,2
            if (k2.lt.k1) go to 570
            do 560 k=k1,k2
               ir = ind(k,l)
               if (l.eq.1) ind(k,l) = 0
               nz = iw(ir,l)
               if (nz.le.0) go to 720
               in = iw(nz,l+2)
               iw(ir,l+6) = in
               iw(ir,l+4) = 0
               iw(nz,l+2) = ir
               if (in.ne.0) iw(in,l+4) = ir
  560       continue
  570       k1 = ip(ipp,1) + 1
            k2 = iw(ipp,1) + k1 - 2
  580    continue
  590 continue
c
c reset column file to refer to u and store row/col numbers in
c     pivotal order in iw(.,3),iw(.,4)
      do 600 i=1,n
         j = -iw(i,5)
         iw(j,3) = i
         j = -iw(i,6)
         iw(j,4) = i
         iw(i,2) = 0
  600 continue
      do 620 i=1,n
         kp = ip(i,1)
         kl = iw(i,1) + kp - 1
         do 610 k=kp,kl
            j = ind(k,2)
            iw(j,2) = iw(j,2) + 1
  610    continue
  620 continue
      k = 1
      do 630 i=1,n
         k = k + iw(i,2)
         ip(i,2) = k
  630 continue
      lcol = k - 1
      do 650 ii=1,n
         i = iw(ii,3)
         kp = ip(i,1)
         kl = iw(i,1) + kp - 1
         do 640 k=kp,kl
            j = ind(k,2)
            kn = ip(j,2) - 1
            ip(j,2) = kn
            ind(kn,1) = i
  640    continue
  650 continue
      return
c
c     the following instructions implement the failure exits.
c
  660 if (lp.gt.0) then
         write (xern1, '(i8)') ir
         write (xern2, '(i8)') j
         call xermsg ('slatec', 'la05as', 'more than one matrix ' //
     *      'entry.  here row = ' // xern1 // ' and col = ' // xern2,
     *      -4, 1)
      endif
      g = -4.
      return
c
  670 if (lp.gt.0) call xermsg ('slatec', 'la05as',
     *   'the order of the system, n, is not positive.', -1, 1)
      g = -1.0e0
      return
c
  680 if (lp.gt.0) then
         write (xern0, '(i8)') k
         write (xern1, '(i8)') i
         write (xern2, '(i8)') j
         call xermsg ('slatec', 'la05as', 'element k = ' // xern0 //
     *      ' is out of bounds.$$here row = ' // xern1 //
     *      ' and col = ' // xern2, -3, 1)
      endif
      g = -3.
      return
c
  700 if (lp.gt.0) then
         write (xern1, '(i8)') l
         call xermsg ('slatec', 'la05as', 'row or column has no ' //
     *      'elements.  here index = ' // xern1, -2, 1)
      endif
      g = -2.
      return
c
  710 if (lp.gt.0) call xermsg ('slatec', 'la05as',
     *   'lengths of arrays a(*) and ind(*,2) are too small.', -7, 1)
      g = -7.
      return
c
  720 ipv = ipv + 1
      iw(ipv,1) = ir
      do 730 i=1,n
         ii = -iw(i,l+4)
         if (ii.gt.0) iw(ii,1) = i
  730 continue
c
      if (lp.gt.0) then
         xern1 = 'rows'
         if (l.eq.2) xern1 = 'columns'
         call xermsg ('slatec', 'la05as', 'dependant ' // xern1, -5, 1)
c
  740    write (xern1, '(i8)') iw(i,1)
         xern2 = ' '
         if (i+1.le.ipv) write (xern2, '(i8)') iw(i+1,1)
         call xermsg ('slatec', 'la05as',
     *      'dependent vector indices are ' // xern1 // ' and ' //
     *      xern2, -5, 1)
         i = i + 2
         if (i.le.ipv) go to 740
      endif
      g = -5.
      return
      end
