*deck wnlit
      subroutine wnlit (w, mdw, m, n, l, ipivot, itype, h, scale, rnorm,
     +   idope, dope, done)
c***begin prologue  wnlit
c***subsidiary
c***purpose  subsidiary to wnnls
c***library   slatec
c***type      single precision (wnlit-s, dwnlit-d)
c***author  hanson, r. j., (snla)
c           haskell, k. h., (snla)
c***description
c
c     this is a companion subprogram to wnnls( ).
c     the documentation for wnnls( ) has complete usage instructions.
c
c     note  the m by (n+1) matrix w( , ) contains the rt. hand side
c           b as the (n+1)st col.
c
c     triangularize l1 by l1 subsystem, where l1=min(m,l), with
c     col interchanges.
c
c***see also  wnnls
c***routines called  h12, isamax, scopy, srotm, srotmg, sscal, sswap,
c                    wnlt1, wnlt2, wnlt3
c***revision history  (yymmdd)
c   790701  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890618  completely restructured and revised.  (wrb & rwc)
c   890620  revised to make wnlt1, wnlt2, and wnlt3 subroutines.  (rwc)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  wnlit
      integer idope(*), ipivot(*), itype(*), l, m, mdw, n
      real             dope(*), h(*), rnorm, scale(*), w(mdw,*)
      logical done
c
      external h12, isamax, scopy, srotm, srotmg, sscal, sswap, wnlt1,
     *   wnlt2, wnlt3
      integer isamax
      logical wnlt2
c
      real             alsq, amax, eanorm, factor, hbar, rn, sparam(5),
     *   t, tau
      integer i, i1, imax, ir, j, j1, jj, jp, krank, l1, lb, lend, me,
     *   mend, niv, nsoln
      logical indep, recalc
c
c***first executable statement  wnlit
      me    = idope(1)
      nsoln = idope(2)
      l1    = idope(3)
c
      alsq   = dope(1)
      eanorm = dope(2)
      tau    = dope(3)
c
      lb     = min(m-1,l)
      recalc = .true.
      rnorm  = 0.e0
      krank  = 0
c
c     we set factor=1.0 so that the heavy weight alamda will be
c     included in the test for column independence.
c
      factor = 1.e0
      lend = l
      do 180 i=1,lb
c
c        set ir to point to the i-th row.
c
         ir = i
         mend = m
         call wnlt1 (i, lend, m, ir, mdw, recalc, imax, hbar, h, scale,
     +                w)
c
c        update column ss and find pivot column.
c
         call wnlt3 (i, imax, m, mdw, ipivot, h, w)
c
c        perform column interchange.
c        test independence of incoming column.
c
  130    if (wnlt2(me, mend, ir, factor, tau, scale, w(1,i))) then
c
c           eliminate i-th column below diagonal using modified givens
c           transformations applied to (a b).
c
c           when operating near the me line, use the largest element
c           above it as the pivot.
c
            do 160 j=m,i+1,-1
               jp = j-1
               if (j.eq.me+1) then
                  imax = me
                  amax = scale(me)*w(me,i)**2
                  do 150 jp=j-1,i,-1
                     t = scale(jp)*w(jp,i)**2
                     if (t.gt.amax) then
                        imax = jp
                        amax = t
                     endif
  150             continue
                  jp = imax
               endif
c
               if (w(j,i).ne.0.e0) then
                  call srotmg (scale(jp), scale(j), w(jp,i), w(j,i),
     +                         sparam)
                  w(j,i) = 0.e0
                  call srotm (n+1-i, w(jp,i+1), mdw, w(j,i+1), mdw,
     +                        sparam)
               endif
  160       continue
         else if (lend.gt.i) then
c
c           column i is dependent.  swap with column lend.
c           perform column interchange,
c           and find column in remaining set with largest ss.
c
            call wnlt3 (i, lend, m, mdw, ipivot, h, w)
            lend = lend - 1
            imax = isamax(lend-i+1, h(i), 1) + i - 1
            hbar = h(imax)
            go to 130
         else
            krank = i - 1
            go to 190
         endif
  180 continue
      krank = l1
c
  190 if (krank.lt.me) then
         factor = alsq
         do 200 i=krank+1,me
            call scopy (l, 0.e0, 0, w(i,1), mdw)
  200    continue
c
c        determine the rank of the remaining equality constraint
c        equations by eliminating within the block of constrained
c        variables.  remove any redundant constraints.
c
         recalc = .true.
         lb = min(l+me-krank, n)
         do 270 i=l+1,lb
            ir = krank + i - l
            lend = n
            mend = me
            call wnlt1 (i, lend, me, ir, mdw, recalc, imax, hbar, h,
     +                   scale, w)
c
c           update col ss and find pivot col
c
            call wnlt3 (i, imax, m, mdw, ipivot, h, w)
c
c           perform column interchange
c           eliminate elements in the i-th col.
c
            do 240 j=me,ir+1,-1
               if (w(j,i).ne.0.e0) then
                 call srotmg (scale(j-1), scale(j), w(j-1,i), w(j,i),
     +                        sparam)
                  w(j,i) = 0.e0
                  call srotm (n+1-i, w(j-1,i+1), mdw,w(j,i+1), mdw,
     +                        sparam)
               endif
  240       continue
c
c           i=column being eliminated.
c           test independence of incoming column.
c           remove any redundant or dependent equality constraints.
c
            if (.not.wnlt2(me, mend, ir, factor,tau,scale,w(1,i))) then
               jj = ir
               do 260 ir=jj,me
                  call scopy (n, 0.e0, 0, w(ir,1), mdw)
                  rnorm = rnorm + (scale(ir)*w(ir,n+1)/alsq)*w(ir,n+1)
                  w(ir,n+1) = 0.e0
                  scale(ir) = 1.e0
c
c                 reclassify the zeroed row as a least squares equation.
c
                  itype(ir) = 1
  260          continue
c
c              reduce me to reflect any discovered dependent equality
c              constraints.
c
               me = jj - 1
               go to 280
            endif
  270    continue
      endif
c
c     try to determine the variables krank+1 through l1 from the
c     least squares equations.  continue the triangularization with
c     pivot element w(me+1,i).
c
  280 if (krank.lt.l1) then
         recalc = .true.
c
c        set factor=alsq to remove effect of heavy weight from
c        test for column independence.
c
         factor = alsq
         do 350 i=krank+1,l1
c
c           set ir to point to the me+1-st row.
c
            ir = me+1
            lend = l
            mend = m
            call wnlt1 (i, l, m, ir, mdw, recalc, imax, hbar, h, scale,
     +                   w)
c
c           update column ss and find pivot column.
c
            call wnlt3 (i, imax, m, mdw, ipivot, h, w)
c
c           perform column interchange.
c           eliminate i-th column below the ir-th element.
c
            do 320 j=m,ir+1,-1
               if (w(j,i).ne.0.e0) then
                 call srotmg (scale(j-1), scale(j), w(j-1,i), w(j,i),
     +                        sparam)
                  w(j,i) = 0.e0
                  call srotm (n+1-i, w(j-1,i+1),  mdw, w(j,i+1), mdw,
     +                        sparam)
               endif
  320       continue
c
c           test if new pivot element is near zero.
c           if so, the column is dependent.
c           then check row norm test to be classified as independent.
c
            t = scale(ir)*w(ir,i)**2
            indep = t .gt. (tau*eanorm)**2
            if (indep) then
               rn = 0.e0
               do 340 i1=ir,m
                  do 330 j1=i+1,n
                     rn = max(rn, scale(i1)*w(i1,j1)**2)
  330             continue
  340          continue
               indep = t .gt. rn*tau**2
            endif
c
c           if independent, swap the ir-th and krank+1-th rows to
c           maintain the triangular form.  update the rank indicator
c           krank and the equality constraint pointer me.
c
            if (.not.indep) go to 360
            call sswap(n+1, w(krank+1,1), mdw, w(ir,1), mdw)
            call sswap(1, scale(krank+1), 1, scale(ir), 1)
c
c           reclassify the least square equation as an equality
c           constraint and rescale it.
c
            itype(ir) = 0
            t = sqrt(scale(krank+1))
            call sscal(n+1, t, w(krank+1,1), mdw)
            scale(krank+1) = alsq
            me = me+1
            krank = krank+1
  350    continue
      endif
c
c     if pseudorank is less than l, apply householder transformation.
c     from right.
c
  360 if (krank.lt.l) then
         do 370 j=krank,1,-1
            call h12 (1, j, krank+1, l, w(j,1), mdw, h(j), w, mdw, 1,
     +                j-1)
  370    continue
      endif
c
      niv = krank + nsoln - l
      if (l.eq.n) done = .true.
c
c     end of initial triangularization.
c
      idope(1) = me
      idope(2) = krank
      idope(3) = niv
      return
      end
