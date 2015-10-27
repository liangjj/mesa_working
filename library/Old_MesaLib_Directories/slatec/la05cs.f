*deck la05cs
      subroutine la05cs (a, ind, ia, n, ip, iw, w, g, u, mm)
c***begin prologue  la05cs
c***subsidiary
c***purpose  subsidiary to splp
c***library   slatec
c***type      single precision (la05cs-s, la05cd-d)
c***author  (unknown)
c***description
c
c     this subprogram is a slight modification of a subprogram
c     from the c. 1979 aere harwell library.  the name of the
c     corresponding harwell code can be obtained by deleting
c     the final letter =s= in the names used here.
c     revised sep. 13, 1979.
c
c     royalties have been paid to aere-uk for use of their codes
c     in the package given here.  any primary usage of the harwell
c     subroutines requires a royalty agreement and payment between
c     the user and aere-uk.  any usage of the sandia written codes
c     splp( ) (which uses the harwell subroutines) is permitted.
c
c***see also  splp
c***routines called  la05es, xermsg, xsetun
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
c   920410  corrected second dimension on iw declaration.  (wrb)
c   920422  changed upper limit on do from last to last-1.  (wrb)
c***end prologue  la05cs
      real a(*), g, u, am, w(*), small, au
      integer ind(ia,2), iw(n,8)
      integer ip(n,2)
      character*8 xern1
c
      common /la05ds/ small, lp, lenl, lenu, ncp, lrow, lcol
c***first executable statement  la05cs
      call xsetun(lp)
      if (g.lt.0.0e0) go to 620
      jm = mm
c mcp limits the value of ncp permitted before an error return results.
      mcp = ncp + 20
c remove old column
      lenu = lenu - iw(jm,2)
      kp = ip(jm,2)
      im = ind(kp,1)
      kl = kp + iw(jm,2) - 1
      iw(jm,2) = 0
      do 30 k=kp,kl
         i = ind(k,1)
         ind(k,1) = 0
         kr = ip(i,1)
         nz = iw(i,1) - 1
         iw(i,1) = nz
         krl = kr + nz
         do 10 km=kr,krl
            if (ind(km,2).eq.jm) go to 20
   10    continue
   20    a(km) = a(krl)
         ind(km,2) = ind(krl,2)
         ind(krl,2) = 0
   30 continue
c
c insert new column
      do 110 ii=1,n
         i = iw(ii,3)
         if (i.eq.im) m = ii
         if (abs(w(i)).le.small) go to 100
         lenu = lenu + 1
         last = ii
         if (lcol+lenl.lt.ia) go to 40
c compress column file if necessary.
         if (ncp.ge.mcp .or. lenl+lenu.ge.ia) go to 610
         call la05es(a, ind, ip(1,2), n, iw(1,2), ia, .false.)
   40    lcol = lcol + 1
         nz = iw(jm,2)
         if (nz.eq.0) ip(jm,2) = lcol
         iw(jm,2) = nz + 1
         ind(lcol,1) = i
         nz = iw(i,1)
         kpl = ip(i,1) + nz
         if (kpl.gt.lrow) go to 50
         if (ind(kpl,2).eq.0) go to 90
c new entry has to be created.
   50    if (lenl+lrow+nz.lt.ia) go to 60
         if (ncp.ge.mcp .or. lenl+lenu+nz.ge.ia) go to 610
c compress row file if necessary.
         call la05es(a, ind(1,2), ip, n, iw, ia, .true.)
   60    kp = ip(i,1)
         ip(i,1) = lrow + 1
         if (nz.eq.0) go to 80
         kpl = kp + nz - 1
         do 70 k=kp,kpl
            lrow = lrow + 1
            a(lrow) = a(k)
            ind(lrow,2) = ind(k,2)
            ind(k,2) = 0
   70    continue
   80    lrow = lrow + 1
         kpl = lrow
c place new element at end of row.
   90    iw(i,1) = nz + 1
         a(kpl) = w(i)
         ind(kpl,2) = jm
  100    w(i) = 0.0e0
  110 continue
      if (iw(im,1).eq.0 .or. iw(jm,2).eq.0 .or. m.gt.last) go to 590
c
c find column singletons, other than the spike. non-singletons are
c     marked with w(j)=1. only iw(.,3) is revised and iw(.,4) is used
c     for workspace.
      ins = m
      m1 = m
      w(jm) = 1.0e0
      do 140 ii=m,last
         i = iw(ii,3)
         j = iw(ii,4)
         if (w(j).eq.0.0e0) go to 130
         kp = ip(i,1)
         kl = kp + iw(i,1) - 1
         do 120 k=kp,kl
            j = ind(k,2)
            w(j) = 1.0e0
  120    continue
         iw(ins,4) = i
         ins = ins + 1
         go to 140
c place singletons in new position.
  130    iw(m1,3) = i
         m1 = m1 + 1
  140 continue
c place non-singletons in new position.
      ij = m + 1
      do 150 ii=m1,last-1
         iw(ii,3) = iw(ij,4)
         ij = ij + 1
  150 continue
c place spike at end.
      iw(last,3) = im
c
c find row singletons, apart from spike row. non-singletons are marked
c     with w(i)=2. again only iw(.,3) is revised and iw(.,4) is used
c     for workspace.
      last1 = last
      jns = last
      w(im) = 2.0e0
      j = jm
      do 180 ij=m1,last
         ii = last + m1 - ij
         i = iw(ii,3)
         if (w(i).ne.2.0e0) go to 170
         k = ip(i,1)
         if (ii.ne.last) j = ind(k,2)
         kp = ip(j,2)
         kl = kp + iw(j,2) - 1
         iw(jns,4) = i
         jns = jns - 1
         do 160 k=kp,kl
            i = ind(k,1)
            w(i) = 2.0e0
  160    continue
         go to 180
  170    iw(last1,3) = i
         last1 = last1 - 1
  180 continue
      do 190 ii=m1,last1
         jns = jns + 1
         i = iw(jns,4)
         w(i) = 3.0e0
         iw(ii,3) = i
  190 continue
c
c deal with singleton spike column. note that bump rows are marked by
c    w(i)=3.0e0
      do 230 ii=m1,last1
         kp = ip(jm,2)
         kl = kp + iw(jm,2) - 1
         is = 0
         do 200 k=kp,kl
            l = ind(k,1)
            if (w(l).ne.3.0e0) go to 200
            if (is.ne.0) go to 240
            i = l
            knp = k
            is = 1
  200    continue
         if (is.eq.0) go to 590
c make a(i,jm) a pivot.
         ind(knp,1) = ind(kp,1)
         ind(kp,1) = i
         kp = ip(i,1)
         do 210 k=kp,ia
            if (ind(k,2).eq.jm) go to 220
  210    continue
  220    am = a(kp)
         a(kp) = a(k)
         a(k) = am
         ind(k,2) = ind(kp,2)
         ind(kp,2) = jm
         jm = ind(k,2)
         iw(ii,4) = i
         w(i) = 2.0e0
  230 continue
      ii = last1
      go to 260
  240 in = m1
      do 250 ij=ii,last1
         iw(ij,4) = iw(in,3)
         in = in + 1
  250 continue
  260 last2 = last1 - 1
      if (m1.eq.last1) go to 570
      do 270 i=m1,last2
         iw(i,3) = iw(i,4)
  270 continue
      m1 = ii
      if (m1.eq.last1) go to 570
c
c clear w
      do 280 i=1,n
         w(i) = 0.0e0
  280 continue
c
c perform elimination
      ir = iw(last1,3)
      do 560 ii=m1,last1
         ipp = iw(ii,3)
         kp = ip(ipp,1)
         kr = ip(ir,1)
         jp = ind(kp,2)
         if (ii.eq.last1) jp = jm
c search non-pivot row for element to be eliminated.
c  and bring it to front of its row
         krl = kr + iw(ir,1) - 1
         do 290 knp=kr,krl
            if (jp.eq.ind(knp,2)) go to 300
  290    continue
         if (ii-last1) 560, 590, 560
c bring element to be eliminated to front of its row.
  300    am = a(knp)
         a(knp) = a(kr)
         a(kr) = am
         ind(knp,2) = ind(kr,2)
         ind(kr,2) = jp
         if (ii.eq.last1) go to 310
         if (abs(a(kp)).lt.u*abs(am)) go to 310
         if (abs(am).lt.u*abs(a(kp))) go to 340
         if (iw(ipp,1).le.iw(ir,1)) go to 340
c perform interchange
  310    iw(last1,3) = ipp
         iw(ii,3) = ir
         ir = ipp
         ipp = iw(ii,3)
         k = kr
         kr = kp
         kp = k
         kj = ip(jp,2)
         do 320 k=kj,ia
            if (ind(k,1).eq.ipp) go to 330
  320    continue
  330    ind(k,1) = ind(kj,1)
         ind(kj,1) = ipp
  340    if (a(kp).eq.0.0e0) go to 590
         if (ii.eq.last1) go to 560
         am = -a(kr)/a(kp)
c compress row file unless it is certain that there is room for new row.
         if (lrow+iw(ir,1)+iw(ipp,1)+lenl.le.ia) go to 350
         if (ncp.ge.mcp .or. lenu+iw(ir,1)+iw(ipp,1)+lenl.gt.ia) go to
     *    610
         call la05es(a, ind(1,2), ip, n, iw, ia, .true.)
         kp = ip(ipp,1)
         kr = ip(ir,1)
  350    krl = kr + iw(ir,1) - 1
         kq = kp + 1
         kpl = kp + iw(ipp,1) - 1
c place pivot row (excluding pivot itself) in w.
         if (kq.gt.kpl) go to 370
         do 360 k=kq,kpl
            j = ind(k,2)
            w(j) = a(k)
  360    continue
  370    ip(ir,1) = lrow + 1
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
  380       lenu = lenu - 1
c remove element from col file.
            k = ip(j,2)
            kl = k + iw(j,2) - 1
            iw(j,2) = kl - k
            do 390 kk=k,kl
               if (ind(kk,1).eq.ir) go to 400
  390       continue
  400       ind(kk,1) = ind(kl,1)
            ind(kl,1) = 0
  410       w(j) = 0.0e0
  420    continue
c
c scan pivot row for fills.
  430    if (kq.gt.kpl) go to 520
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
c if possible place new element at end of present entry.
            if (kl.ne.lcol) go to 440
            if (lcol+lenl.ge.ia) go to 460
            lcol = lcol + 1
            go to 450
  440       if (ind(kl+1,1).ne.0) go to 460
  450       ind(kl+1,1) = ir
            go to 490
c new entry has to be created.
  460       if (lcol+lenl+nz+1.lt.ia) go to 470
c compress column file if there is not room for new entry.
            if (ncp.ge.mcp .or. lenu+lenl+nz+1.ge.ia) go to 610
            call la05es(a, ind, ip(1,2), n, iw(1,2), ia, .false.)
            k = ip(j,2)
            kl = k + nz - 1
c transfer old entry into new.
  470       ip(j,2) = lcol + 1
            do 480 kk=k,kl
               lcol = lcol + 1
               ind(lcol,1) = ind(kk,1)
               ind(kk,1) = 0
  480       continue
c add new element.
            lcol = lcol + 1
            ind(lcol,1) = ir
  490       g = max(g,abs(au))
            iw(j,2) = nz + 1
  500       w(j) = 0.0e0
  510    continue
  520    iw(ir,1) = lrow + 1 - ip(ir,1)
c
c store multiplier
         if (lenl+lcol+1.le.ia) go to 530
c compress col file if necessary.
         if (ncp.ge.mcp) go to 610
         call la05es(a, ind, ip(1,2), n, iw(1,2), ia, .false.)
  530    k = ia - lenl
         lenl = lenl + 1
         a(k) = am
         ind(k,1) = ipp
         ind(k,2) = ir
c create blank in pivotal column.
         kp = ip(jp,2)
         nz = iw(jp,2) - 1
         kl = kp + nz
         do 540 k=kp,kl
            if (ind(k,1).eq.ir) go to 550
  540    continue
  550    ind(k,1) = ind(kl,1)
         iw(jp,2) = nz
         ind(kl,1) = 0
         lenu = lenu - 1
  560 continue
c
c construct column permutation and store it in iw(.,4)
  570 do 580 ii=m,last
         i = iw(ii,3)
         k = ip(i,1)
         j = ind(k,2)
         iw(ii,4) = j
  580 continue
      return
c
c     the following instructions implement the failure exits.
c
  590 if (lp.gt.0) then
         write (xern1, '(i8)') mm
         call xermsg ('slatec', 'la05cs', 'singular matrix after ' //
     *      'replacement of column.  index = ' // xern1, -6, 1)
      endif
      g = -6.0e0
      return
c
  610 if (lp.gt.0) call xermsg ('slatec', 'la05cs',
     *   'lengths of arrays a(*) and ind(*,2) are too small.', -7, 1)
      g = -7.0e0
      return
c
  620 if (lp.gt.0) call xermsg ('slatec', 'la05cs',
     *   'earlier entry gave error return.', -8, 2)
      g = -8.0e0
      return
      end
