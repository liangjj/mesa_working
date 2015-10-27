*deck efcmn
      subroutine efcmn (ndata, xdata, ydata, sddata, nord, nbkpt,
     +   bkptin, mdein, mdeout, coeff, bf, xtemp, ptemp, bkpt, g, mdg,
     +   w, mdw, lw)
c***begin prologue  efcmn
c***subsidiary
c***purpose  subsidiary to efc
c***library   slatec
c***type      single precision (efcmn-s, defcmn-d)
c***author  hanson, r. j., (snla)
c***description
c
c     this is a companion subprogram to efc( ).
c     this subprogram does weighted least squares fitting of data by
c     b-spline curves.
c     the documentation for efc( ) has complete usage instructions.
c
c***see also  efc
c***routines called  bndacc, bndsol, bsplvn, scopy, sscal, ssort, xermsg
c***revision history  (yymmdd)
c   800801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890618  completely restructured and extensively revised (wrb & rwc)
c   890831  modified array declarations.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900315  calls to xerror changed to calls to xermsg.  (thj)
c   900328  added type section.  (wrb)
c   900510  convert xerrwv calls to xermsg calls.  (rwc)
c***end prologue  efcmn
      integer lw, mdein, mdeout, mdg, mdw, nbkpt, ndata, nord
      real             bf(nord,*), bkpt(*), bkptin(*), coeff(*),
     *   g(mdg,*), ptemp(*), sddata(*), w(mdw,*), xdata(*), xtemp(*),
     *   ydata(*)
c
      external bndacc, bndsol, bsplvn, scopy, sscal, ssort, xermsg
c
      real             dummy, rnorm, xmax, xmin, xval
      integer i, idata, ileft, intseq, ip, ir, irow, l, mt, n, nb,
     *   nordm1, nordp1, np1
      character*8 xern1, xern2
c
c***first executable statement  efcmn
c
c     initialize variables and analyze input.
c
      n = nbkpt - nord
      np1 = n + 1
c
c     initially set all output coefficients to zero.
c
      call scopy (n, 0.e0, 0, coeff, 1)
      mdeout = -1
      if (nord.lt.1 .or. nord.gt.20) then
         call xermsg ('slatec', 'efcmn',
     +      'in efc, the order of the b-spline must be 1 thru 20.',
     +      3, 1)
         return
      endif
c
      if (nbkpt.lt.2*nord) then
         call xermsg ('slatec', 'efcmn',
     +      'in efc, the number of knots must be at least twice ' //
     +      'the b-spline order.', 4, 1)
         return
      endif
c
      if (ndata.lt.0) then
         call xermsg ('slatec', 'efcmn',
     +      'in efc, the number of data points must be nonnegative.',
     +      5, 1)
         return
      endif
c
      nb = (nbkpt-nord+3)*(nord+1) + (nbkpt+1)*(nord+1) +
     +     2*max(nbkpt,ndata) + nbkpt + nord**2
      if (lw .lt. nb) then
         write (xern1, '(i8)') nb
         write (xern2, '(i8)') lw
         call xermsg ('slatec', 'efcmn',
     *      'in efc, insufficient storage for w(*).  check formula ' //
     *      'that reads lw.ge. ... .  need = ' // xern1 //
     *      ' given = ' // xern2, 6, 1)
         mdeout = -1
         return
      endif
c
      if (mdein.ne.1 .and. mdein.ne.2) then
         call xermsg ('slatec', 'efcmn',
     +      'in efc, input value of mdein must be 1-2.', 7, 1)
         return
      endif
c
c     sort the breakpoints.
c
      call scopy (nbkpt, bkptin, 1, bkpt, 1)
      call ssort (bkpt, dummy, nbkpt, 1)
c
c     save interval containing knots.
c
      xmin = bkpt(nord)
      xmax = bkpt(np1)
      nordm1 = nord - 1
      nordp1 = nord + 1
c
c     process least squares equations.
c
c     sort data and an array of pointers.
c
      call scopy (ndata, xdata, 1, xtemp, 1)
      do 100 i = 1,ndata
         ptemp(i) = i
  100 continue
c
      if (ndata.gt.0) then
         call ssort (xtemp, ptemp, ndata, 2)
         xmin = min(xmin,xtemp(1))
         xmax = max(xmax,xtemp(ndata))
      endif
c
c     fix breakpoint array if needed. this should only involve very
c     minor differences with the input array of breakpoints.
c
      do 110 i = 1,nord
         bkpt(i) = min(bkpt(i),xmin)
  110 continue
c
      do 120 i = np1,nbkpt
         bkpt(i) = max(bkpt(i),xmax)
  120 continue
c
c     initialize parameters of banded matrix processor, bndacc( ).
c
      mt = 0
      ip = 1
      ir = 1
      ileft = nord
      intseq = 1
      do 150 idata = 1,ndata
c
c        sorted indices are in ptemp(*).
c
         l = ptemp(idata)
         xval = xdata(l)
c
c        when interval changes, process equations in the last block.
c
         if (xval.ge.bkpt(ileft+1)) then
            call bndacc (g, mdg, nord, ip, ir, mt, ileft-nordm1)
            mt = 0
c
c           move pointer up to have bkpt(ileft).le.xval, ileft.le.n.
c
            do 130 ileft = ileft,n
               if (xval.lt.bkpt(ileft+1)) go to 140
               if (mdein.eq.2) then
c
c                 data is being sequentially accumulated.
c                 transfer previously accumulated rows from w(*,*) to
c                 g(*,*) and process them.
c
                  call scopy (nordp1, w(intseq,1), mdw, g(ir,1), mdg)
                  call bndacc (g, mdg, nord, ip, ir, 1, intseq)
                  intseq = intseq + 1
               endif
  130       continue
         endif
c
c        obtain b-spline function value.
c
  140    call bsplvn (bkpt, nord, 1, xval, ileft, bf)
c
c        move row into place.
c
         irow = ir + mt
         mt = mt + 1
         call scopy (nord, bf, 1, g(irow,1), mdg)
         g(irow,nordp1) = ydata(l)
c
c        scale data if uncertainty is nonzero.
c
         if (sddata(l).ne.0.e0) call sscal (nordp1, 1.e0/sddata(l),
     +                               g(irow,1), mdg)
c
c        when staging work area is exhausted, process rows.
c
         if (irow.eq.mdg-1) then
            call bndacc (g, mdg, nord, ip, ir, mt, ileft-nordm1)
            mt = 0
         endif
  150 continue
c
c     process last block of equations.
c
      call bndacc (g, mdg, nord, ip, ir, mt, ileft-nordm1)
c
c     finish processing any previously accumulated rows from w(*,*)
c     to g(*,*).
c
      if (mdein.eq.2) then
         do 160 i = intseq,np1
            call scopy (nordp1, w(i,1), mdw, g(ir,1), mdg)
            call bndacc (g, mdg, nord, ip, ir, 1, min(n,i))
  160    continue
      endif
c
c     last call to adjust block positioning.
c
      call scopy (nordp1, 0.e0, 0, g(ir,1), mdg)
      call bndacc (g, mdg, nord, ip, ir, 1, np1)
c
c     transfer accumulated rows from g(*,*) to w(*,*) for
c     possible later sequential accumulation.
c
      do 170 i = 1,np1
         call scopy (nordp1, g(i,1), mdg, w(i,1), mdw)
  170 continue
c
c     solve for coefficients when possible.
c
      do 180 i = 1,n
         if (g(i,1).eq.0.e0) then
            mdeout = 2
            return
         endif
  180 continue
c
c     all the diagonal terms in the accumulated triangular
c     matrix are nonzero.  the solution can be computed but
c     it may be unsuitable for further use due to poor
c     conditioning or the lack of constraints.  no checking
c     for either of these is done here.
c
      call bndsol (1, g, mdg, nord, ip, ir, coeff, n, rnorm)
      mdeout = 1
      return
      end
