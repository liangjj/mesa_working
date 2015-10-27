*deck pjac
      subroutine pjac (neq, y, yh, nyh, ewt, ftem, savf, wm, iwm, f,
     +   jac, rpar, ipar)
c***begin prologue  pjac
c***subsidiary
c***purpose  subsidiary to debdf
c***library   slatec
c***type      single precision (pjac-s, dpjac-d)
c***author  watts, h. a., (snla)
c***description
c
c   pjac sets up the iteration matrix (involving the jacobian) for the
c   integration package debdf.
c
c***see also  debdf
c***routines called  sgbfa, sgefa, vnwrms
c***common blocks    debdf1
c***revision history  (yymmdd)
c   800901  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c   920422  changed dimension statement.  (wrb)
c***end prologue  pjac
c
clll. optimize
      integer neq, nyh, iwm, i, i1, i2, ier, ii, iownd, iowns, j, j1,
     1   jj, jstart, kflag, l, lenp, maxord, mba, mband, meb1, meband,
     2   meth, miter, ml, ml3, mu, n, nfe, nje, nq, nqu, nst
      external f, jac
      real y, yh, ewt, ftem, savf, wm,
     1   rownd, rowns, el0, h, hmin, hmxi, hu, tn, uround,
     2   con, di, fac, hl0, r, r0, srur, yi, yj, yjj, vnwrms
      dimension         y(*), yh(nyh,*), ewt(*), ftem(*), savf(*),
     1   wm(*), iwm(*), rpar(*), ipar(*)
      common /debdf1/ rownd, rowns(210),
     1   el0, h, hmin, hmxi, hu, tn, uround, iownd(14), iowns(6),
     2   ier, jstart, kflag, l, meth, miter, maxord, n, nq, nst, nfe,
     3   nje, nqu
c-----------------------------------------------------------------------
c pjac is called by stod  to compute and process the matrix
c p = i - h*el(1)*j , where j is an approximation to the jacobian.
c here j is computed by the user-supplied routine jac if
c miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
c if miter = 3, a diagonal approximation to j is used.
c j is stored in wm and replaced by p.  if miter .ne. 3, p is then
c subjected to lu decomposition in preparation for later solution
c of linear systems with p as coefficient matrix. this is done
c by sgefa if miter = 1 or 2, and by sgbfa if miter = 4 or 5.
c
c in addition to variables described previously, communication
c with pjac uses the following..
c y    = array containing predicted values on entry.
c ftem = work array of length n (acor in stod ).
c savf = array containing f evaluated at predicted y.
c wm   = real work space for matrices.  on output it contains the
c        inverse diagonal matrix if miter = 3 and the lu decomposition
c        of p if miter is 1, 2 , 4, or 5.
c        storage of matrix elements starts at wm(3).
c        wm also contains the following matrix-related data..
c        wm(1) = sqrt(uround), used in numerical jacobian increments.
c        wm(2) = h*el0, saved for later use if miter = 3.
c iwm  = integer work space containing pivot information, starting at
c        iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains the
c        band parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c el0  = el(1) (input).
c ier  = output error flag,  = 0 if no trouble, .ne. 0 if
c        p matrix found to be singular.
c this routine also uses the common variables el0, h, tn, uround,
c miter, n, nfe, and nje.
c-----------------------------------------------------------------------
c***first executable statement  pjac
      nje = nje + 1
      hl0 = h*el0
      go to (100, 200, 300, 400, 500), miter
c if miter = 1, call jac and multiply by scalar. -----------------------
 100  lenp = n*n
      do 110 i = 1,lenp
 110    wm(i+2) = 0.0e0
      call jac (tn, y, wm(3), n, rpar, ipar)
      con = -hl0
      do 120 i = 1,lenp
 120    wm(i+2) = wm(i+2)*con
      go to 240
c if miter = 2, make n calls to f to approximate j. --------------------
 200  fac = vnwrms (n, savf, ewt)
      r0 = 1000.0e0*abs(h)*uround*n*fac
      if (r0 .eq. 0.0e0) r0 = 1.0e0
      srur = wm(1)
      j1 = 2
      do 230 j = 1,n
        yj = y(j)
        r = max(srur*abs(yj),r0*ewt(j))
        y(j) = y(j) + r
        fac = -hl0/r
        call f (tn, y, ftem, rpar, ipar)
        do 220 i = 1,n
 220      wm(i+j1) = (ftem(i) - savf(i))*fac
        y(j) = yj
        j1 = j1 + n
 230    continue
      nfe = nfe + n
c add identity matrix. -------------------------------------------------
 240  j = 3
      do 250 i = 1,n
        wm(j) = wm(j) + 1.0e0
 250    j = j + (n + 1)
c do lu decomposition on p. --------------------------------------------
      call sgefa (wm(3), n, n, iwm(21), ier)
      return
c if miter = 3, construct a diagonal approximation to j and p. ---------
 300  wm(2) = hl0
      ier = 0
      r = el0*0.1e0
      do 310 i = 1,n
 310    y(i) = y(i) + r*(h*savf(i) - yh(i,2))
      call f (tn, y, wm(3), rpar, ipar)
      nfe = nfe + 1
      do 320 i = 1,n
        r0 = h*savf(i) - yh(i,2)
        di = 0.1e0*r0 - h*(wm(i+2) - savf(i))
        wm(i+2) = 1.0e0
        if (abs(r0) .lt. uround*ewt(i)) go to 320
        if (abs(di) .eq. 0.0e0) go to 330
        wm(i+2) = 0.1e0*r0/di
 320    continue
      return
 330  ier = -1
      return
c if miter = 4, call jac and multiply by scalar. -----------------------
 400  ml = iwm(1)
      mu = iwm(2)
      ml3 =  3
      mband = ml + mu + 1
      meband = mband + ml
      lenp = meband*n
      do 410 i = 1,lenp
 410    wm(i+2) = 0.0e0
      call jac (tn, y, wm(ml3), meband, rpar, ipar)
      con = -hl0
      do 420 i = 1,lenp
 420    wm(i+2) = wm(i+2)*con
      go to 570
c if miter = 5, make mband calls to f to approximate j. ----------------
 500  ml = iwm(1)
      mu = iwm(2)
      mband = ml + mu + 1
      mba = min(mband,n)
      meband = mband + ml
      meb1 = meband - 1
      srur = wm(1)
      fac = vnwrms (n, savf, ewt)
      r0 = 1000.0e0*abs(h)*uround*n*fac
      if (r0 .eq. 0.0e0) r0 = 1.0e0
      do 560 j = 1,mba
        do 530 i = j,n,mband
          yi = y(i)
          r = max(srur*abs(yi),r0*ewt(i))
 530      y(i) = y(i) + r
        call f (tn, y, ftem, rpar, ipar)
        do 550 jj = j,n,mband
          y(jj) = yh(jj,1)
          yjj = y(jj)
          r = max(srur*abs(yjj),r0*ewt(jj))
          fac = -hl0/r
          i1 = max(jj-mu,1)
          i2 = min(jj+ml,n)
          ii = jj*meb1 - ml + 2
          do 540 i = i1,i2
 540        wm(ii+i) = (ftem(i) - savf(i))*fac
 550      continue
 560    continue
      nfe = nfe + mba
c add identity matrix. -------------------------------------------------
 570  ii = mband + 2
      do 580 i = 1,n
        wm(ii) = wm(ii) + 1.0e0
 580    ii = ii + meband
c do lu decomposition of p. --------------------------------------------
      call sgbfa (wm(3), meband, n, ml, mu, iwm(21), ier)
      return
c----------------------- end of subroutine pjac -----------------------
      end
