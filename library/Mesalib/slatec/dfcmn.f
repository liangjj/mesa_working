*deck dfcmn
      subroutine dfcmn (ndata, xdata, ydata, sddata, nord, nbkpt,
     +   bkptin, nconst, xconst, yconst, nderiv, mode, coeff, bf, xtemp,
     +   ptemp, bkpt, g, mdg, w, mdw, work, iwork)
c***begin prologue  dfcmn
c***subsidiary
c***purpose  subsidiary to fc
c***library   slatec
c***type      double precision (fcmn-s, dfcmn-d)
c***author  (unknown)
c***description
c
c     this is a companion subprogram to dfc( ).
c     the documentation for dfc( ) has complete usage instructions.
c
c***see also  dfc
c***routines called  daxpy, dbndac, dbndsl, dcopy, dfspvd, dfspvn,
c                    dlsei, dscal, dsort, xermsg
c***revision history  (yymmdd)
c   780801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890618  completely restructured and extensively revised (wrb & rwc)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c   900604  dp version created from sp version.  (rwc)
c***end prologue  dfcmn
      integer iwork(*), mdg, mdw, mode, nbkpt, nconst, ndata, nderiv(*),
     *   nord
      double precision bf(nord,*), bkpt(*), bkptin(*), coeff(*),
     *   g(mdg,*), ptemp(*), sddata(*), w(mdw,*), work(*),
     *   xconst(*), xdata(*), xtemp(*), yconst(*), ydata(*)
c
      external daxpy, dbndac, dbndsl, dcopy, dfspvd, dfspvn, dlsei,
     *   dscal, dsort, xermsg
c
      double precision dummy, prgopt(10), rnorm, rnorme, rnorml, xmax,
     *   xmin, xval, yval
      integer i, idata, ideriv, ileft, intrvl, intw1, ip, ir, irow,
     *   itype, iw1, iw2, l, lw, mt, n, nb, neqcon, nincon, nordm1,
     *   nordp1, np1
      logical band, new, var
      character*8 xern1
c
c***first executable statement  dfcmn
c
c     analyze input.
c
      if (nord.lt.1 .or. nord.gt.20) then
         call xermsg ('slatec', 'dfcmn',
     +      'in dfc, the order of the b-spline must be 1 thru 20.',
     +      2, 1)
         mode = -1
         return
c
      elseif (nbkpt.lt.2*nord) then
         call xermsg ('slatec', 'dfcmn',
     +      'in dfc, the number of knots must be at least twice ' //
     +      'the b-spline order.', 2, 1)
         mode = -1
         return
      endif
c
      if (ndata.lt.0) then
         call xermsg ('slatec', 'dfcmn',
     +      'in dfc, the number of data points must be nonnegative.',
     +      2, 1)
         mode = -1
         return
      endif
c
c     amount of storage allocated for w(*), iw(*).
c
      iw1 = iwork(1)
      iw2 = iwork(2)
      nb = (nbkpt-nord+3)*(nord+1) + 2*max(ndata,nbkpt) + nbkpt +
     +     nord**2
c
c     see if sufficient storage has been allocated.
c
      if (iw1.lt.nb) then
         write (xern1, '(i8)') nb
         call xermsg ('slatec', 'dfcmn',
     *      'in dfc, insufficient storage for w(*).  check nb = ' //
     *      xern1, 2, 1)
         mode = -1
         return
      endif
c
      if (mode.eq.1) then
         band = .true.
         var = .false.
         new = .true.
      elseif (mode.eq.2) then
         band = .false.
         var = .true.
         new = .true.
      elseif (mode.eq.3) then
         band = .true.
         var = .false.
         new = .false.
      elseif (mode.eq.4) then
         band = .false.
         var = .true.
         new = .false.
      else
         call xermsg ('slatec', 'dfcmn',
     +      'in dfc, input value of mode must be 1-4.', 2, 1)
         mode = -1
         return
      endif
      mode = 0
c
c     sort the breakpoints.
c
      call dcopy (nbkpt, bkptin, 1, bkpt, 1)
      call dsort (bkpt, dummy, nbkpt, 1)
c
c     initialize variables.
c
      neqcon = 0
      nincon = 0
      do 100 i = 1,nconst
         l = nderiv(i)
         itype = mod(l,4)
         if (itype.lt.2) then
            nincon = nincon + 1
         else
            neqcon = neqcon + 1
         endif
  100 continue
c
c     compute the number of variables.
c
      n = nbkpt - nord
      np1 = n + 1
      lw = nb + (np1+nconst)*np1 + 2*(neqcon+np1) + (nincon+np1) +
     +     (nincon+2)*(np1+6)
      intw1 = nincon + 2*np1
c
c     save interval containing knots.
c
      xmin = bkpt(nord)
      xmax = bkpt(np1)
c
c     find the smallest referenced independent variable value in any
c     constraint.
c
      do 110 i = 1,nconst
         xmin = min(xmin,xconst(i))
         xmax = max(xmax,xconst(i))
  110 continue
      nordm1 = nord - 1
      nordp1 = nord + 1
c
c     define the option vector prgopt(1-10) for use in dlsei( ).
c
      prgopt(1) = 4
c
c     set the covariance matrix computation flag.
c
      prgopt(2) = 1
      if (var) then
         prgopt(3) = 1
      else
         prgopt(3) = 0
      endif
c
c     increase the rank determination tolerances for both equality
c     constraint equations and least squares equations.
c
      prgopt(4) = 7
      prgopt(5) = 4
      prgopt(6) = 1.d-4
c
      prgopt(7) = 10
      prgopt(8) = 5
      prgopt(9) = 1.d-4
c
      prgopt(10) = 1
c
c     turn off work array length checking in dlsei( ).
c
      iwork(1) = 0
      iwork(2) = 0
c
c     initialize variables and analyze input.
c
      if (new) then
c
c        to process least squares equations sort data and an array of
c        pointers.
c
         call dcopy (ndata, xdata, 1, xtemp, 1)
         do 120 i = 1,ndata
            ptemp(i) = i
  120    continue
c
         if (ndata.gt.0) then
            call dsort (xtemp, ptemp, ndata, 2)
            xmin = min(xmin,xtemp(1))
            xmax = max(xmax,xtemp(ndata))
         endif
c
c        fix breakpoint array if needed.
c
         do 130 i = 1,nord
            bkpt(i) = min(bkpt(i),xmin)
  130    continue
c
         do 140 i = np1,nbkpt
            bkpt(i) = max(bkpt(i),xmax)
  140    continue
c
c        initialize parameters of banded matrix processor, dbndac( ).
c
         mt = 0
         ip = 1
         ir = 1
         ileft = nord
         do 160 idata = 1,ndata
c
c           sorted indices are in ptemp(*).
c
            l = ptemp(idata)
            xval = xdata(l)
c
c           when interval changes, process equations in the last block.
c
            if (xval.ge.bkpt(ileft+1)) then
               call dbndac (g, mdg, nord, ip, ir, mt, ileft-nordm1)
               mt = 0
c
c              move pointer up to have bkpt(ileft).le.xval,
c                 ileft.lt.np1.
c
  150          if (xval.ge.bkpt(ileft+1) .and. ileft.lt.n) then
                  ileft = ileft + 1
                  go to 150
               endif
            endif
c
c           obtain b-spline function value.
c
            call dfspvn (bkpt, nord, 1, xval, ileft, bf)
c
c           move row into place.
c
            irow = ir + mt
            mt = mt + 1
            call dcopy (nord, bf, 1, g(irow,1), mdg)
            g(irow,nordp1) = ydata(l)
c
c           scale data if uncertainty is nonzero.
c
            if (sddata(l).ne.0.d0) call dscal (nordp1, 1.d0/sddata(l),
     +                                  g(irow,1), mdg)
c
c           when staging work area is exhausted, process rows.
c
            if (irow.eq.mdg-1) then
               call dbndac (g, mdg, nord, ip, ir, mt, ileft-nordm1)
               mt = 0
            endif
  160    continue
c
c        process last block of equations.
c
         call dbndac (g, mdg, nord, ip, ir, mt, ileft-nordm1)
c
c        last call to adjust block positioning.
c
         call dcopy (nordp1, 0.d0, 0, g(ir,1), mdg)
         call dbndac (g, mdg, nord, ip, ir, 1, np1)
      endif
c
      band = band .and. nconst.eq.0
      do 170 i = 1,n
         band = band .and. g(i,1).ne.0.d0
  170 continue
c
c     process banded least squares equations.
c
      if (band) then
         call dbndsl (1, g, mdg, nord, ip, ir, coeff, n, rnorm)
         return
      endif
c
c     check further for sufficient storage in working arrays.
c
      if (iw1.lt.lw) then
         write (xern1, '(i8)') lw
         call xermsg ('slatec', 'dfcmn',
     *      'in dfc, insufficient storage for w(*).  check lw = ' //
     *      xern1, 2, 1)
         mode = -1
         return
      endif
c
      if (iw2.lt.intw1) then
         write (xern1, '(i8)') intw1
         call xermsg ('slatec', 'dfcmn',
     *      'in dfc, insufficient storage for iw(*).  check iw1 = ' //
     *      xern1, 2, 1)
         mode = -1
         return
      endif
c
c     write equality constraints.
c     analyze constraint indicators for an equality constraint.
c
      neqcon = 0
      do 220 idata = 1,nconst
         l = nderiv(idata)
         itype = mod(l,4)
         if (itype.gt.1) then
            ideriv = l/4
            neqcon = neqcon + 1
            ileft = nord
            xval = xconst(idata)
c
  180       if (xval.lt.bkpt(ileft+1) .or. ileft.ge.n) go to 190
            ileft = ileft + 1
            go to 180
c
  190       call dfspvd (bkpt, nord, xval, ileft, bf, ideriv+1)
            call dcopy (np1, 0.d0, 0, w(neqcon,1), mdw)
            call dcopy (nord, bf(1,ideriv+1), 1, w(neqcon,ileft-nordm1),
     +                  mdw)
c
            if (itype.eq.2) then
               w(neqcon,np1) = yconst(idata)
            else
               ileft = nord
               yval = yconst(idata)
c
  200          if (yval.lt.bkpt(ileft+1) .or. ileft.ge.n) go to 210
               ileft = ileft + 1
               go to 200
c
  210          call dfspvd (bkpt, nord, yval, ileft, bf, ideriv+1)
               call daxpy (nord, -1.d0, bf(1, ideriv+1), 1,
     +                     w(neqcon, ileft-nordm1), mdw)
            endif
         endif
  220 continue
c
c     transfer least squares data.
c
      do 230 i = 1,np1
         irow = i + neqcon
         call dcopy (n, 0.d0, 0, w(irow,1), mdw)
         call dcopy (min(np1-i, nord), g(i,1), mdg, w(irow,i), mdw)
         w(irow,np1) = g(i,nordp1)
  230 continue
c
c     write inequality constraints.
c     analyze constraint indicators for inequality constraints.
c
      nincon = 0
      do 260 idata = 1,nconst
         l = nderiv(idata)
         itype = mod(l,4)
         if (itype.lt.2) then
            ideriv = l/4
            nincon = nincon + 1
            ileft = nord
            xval = xconst(idata)
c
  240       if (xval.lt.bkpt(ileft+1) .or. ileft.ge.n) go to 250
            ileft = ileft + 1
            go to 240
c
  250       call dfspvd (bkpt, nord, xval, ileft, bf, ideriv+1)
            irow = neqcon + np1 + nincon
            call dcopy (n, 0.d0, 0, w(irow,1), mdw)
            intrvl = ileft - nordm1
            call dcopy (nord, bf(1, ideriv+1), 1, w(irow, intrvl), mdw)
c
            if (itype.eq.1) then
               w(irow,np1) = yconst(idata)
            else
               w(irow,np1) = -yconst(idata)
               call dscal (nord, -1.d0, w(irow, intrvl), mdw)
            endif
         endif
  260 continue
c
c     solve constrained least squares equations.
c
      call dlsei(w, mdw, neqcon, np1, nincon, n, prgopt, coeff, rnorme,
     +          rnorml, mode, work, iwork)
      return
      end
