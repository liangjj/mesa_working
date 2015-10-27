*deck dpjac
      subroutine dpjac (neq, y, yh, nyh, ewt, ftem, savf, wm, iwm, df,
     +   djac, rpar, ipar)
c***begin prologue  dpjac
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (pjac-s, dpjac-d)
c***author  watts, h. a., (snla)
c***description
c
c   dpjac sets up the iteration matrix (involving the jacobian) for the
c   integration package ddebdf.
c
c***see also  ddebdf
c***routines called  dgbfa, dgefa, dvnrms
c***common blocks    ddebd1
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c   920422  changed dimension statement.  (wrb)
c***end prologue  dpjac
c
      integer i, i1, i2, ier, ii, iownd, iowns, ipar, iwm, j, j1,
     1      jj, jstart, kflag, l, lenp, maxord, mba, mband,
     2      meb1, meband, meth, miter, ml, ml3, mu, n, neq,
     3      nfe, nje, nq, nqu, nst, nyh
      double precision con, di, dvnrms, el0, ewt,
     1      fac, ftem, h, hl0, hmin, hmxi, hu, r, r0, rownd, rowns,
     2      rpar, savf, srur, tn, uround, wm, y, yh, yi, yj, yjj
      external df, djac
      dimension y(*),yh(nyh,*),ewt(*),ftem(*),savf(*),wm(*),iwm(*),
     1          rpar(*),ipar(*)
      common /ddebd1/ rownd,rowns(210),el0,h,hmin,hmxi,hu,tn,uround,
     1                iownd(14),iowns(6),ier,jstart,kflag,l,meth,miter,
     2                maxord,n,nq,nst,nfe,nje,nqu
c     ------------------------------------------------------------------
c      dpjac is called by dstod  to compute and process the matrix
c      p = i - h*el(1)*j , where j is an approximation to the jacobian.
c      here j is computed by the user-supplied routine djac if
c      miter = 1 or 4, or by finite differencing if miter = 2, 3, or 5.
c      if miter = 3, a diagonal approximation to j is used.
c      j is stored in wm and replaced by p.  if miter .ne. 3, p is then
c      subjected to lu decomposition in preparation for later solution
c      of linear systems with p as coefficient matrix. this is done
c      by dgefa if miter = 1 or 2, and by dgbfa if miter = 4 or 5.
c
c      in addition to variables described previously, communication
c      with dpjac uses the following..
c      y    = array containing predicted values on entry.
c      ftem = work array of length n (acor in dstod ).
c      savf = array containing df evaluated at predicted y.
c      wm   = double precision work space for matrices.  on output it
c      contains the
c             inverse diagonal matrix if miter = 3 and the lu
c             decomposition of p if miter is 1, 2 , 4, or 5.
c             storage of matrix elements starts at wm(3).
c             wm also contains the following matrix-related data..
c             wm(1) = sqrt(uround), used in numerical jacobian
c             increments.  wm(2) = h*el0, saved for later use if miter =
c             3.
c      iwm  = integer work space containing pivot information, starting
c             at iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains
c             the band parameters ml = iwm(1) and mu = iwm(2) if miter
c             is 4 or 5.
c      el0  = el(1) (input).
c      ier  = output error flag,  = 0 if no trouble, .ne. 0 if
c             p matrix found to be singular.
c      this routine also uses the common variables el0, h, tn, uround,
c      miter, n, nfe, and nje.
c-----------------------------------------------------------------------
c     begin block permitting ...exits to 240
c        begin block permitting ...exits to 220
c           begin block permitting ...exits to 130
c              begin block permitting ...exits to 70
c***first executable statement  dpjac
                  nje = nje + 1
                  hl0 = h*el0
                  go to (10,40,90,140,170), miter
c                 if miter = 1, call djac and multiply by scalar.
c                 -----------------------
   10             continue
                  lenp = n*n
                  do 20 i = 1, lenp
                     wm(i+2) = 0.0d0
   20             continue
                  call djac(tn,y,wm(3),n,rpar,ipar)
                  con = -hl0
                  do 30 i = 1, lenp
                     wm(i+2) = wm(i+2)*con
   30             continue
c              ...exit
                  go to 70
c                 if miter = 2, make n calls to df to approximate j.
c                 --------------------
   40             continue
                  fac = dvnrms(n,savf,ewt)
                  r0 = 1000.0d0*abs(h)*uround*n*fac
                  if (r0 .eq. 0.0d0) r0 = 1.0d0
                  srur = wm(1)
                  j1 = 2
                  do 60 j = 1, n
                     yj = y(j)
                     r = max(srur*abs(yj),r0*ewt(j))
                     y(j) = y(j) + r
                     fac = -hl0/r
                     call df(tn,y,ftem,rpar,ipar)
                     do 50 i = 1, n
                        wm(i+j1) = (ftem(i) - savf(i))*fac
   50                continue
                     y(j) = yj
                     j1 = j1 + n
   60             continue
                  nfe = nfe + n
   70          continue
c              add identity matrix.
c              -------------------------------------------------
               j = 3
               do 80 i = 1, n
                  wm(j) = wm(j) + 1.0d0
                  j = j + (n + 1)
   80          continue
c              do lu decomposition on p.
c              --------------------------------------------
               call dgefa(wm(3),n,n,iwm(21),ier)
c     .........exit
               go to 240
c              if miter = 3, construct a diagonal approximation to j and
c              p. ---------
   90          continue
               wm(2) = hl0
               ier = 0
               r = el0*0.1d0
               do 100 i = 1, n
                  y(i) = y(i) + r*(h*savf(i) - yh(i,2))
  100          continue
               call df(tn,y,wm(3),rpar,ipar)
               nfe = nfe + 1
               do 120 i = 1, n
                  r0 = h*savf(i) - yh(i,2)
                  di = 0.1d0*r0 - h*(wm(i+2) - savf(i))
                  wm(i+2) = 1.0d0
                  if (abs(r0) .lt. uround*ewt(i)) go to 110
c           .........exit
                     if (abs(di) .eq. 0.0d0) go to 130
                     wm(i+2) = 0.1d0*r0/di
  110             continue
  120          continue
c     .........exit
               go to 240
  130       continue
            ier = -1
c     ......exit
            go to 240
c           if miter = 4, call djac and multiply by scalar.
c           -----------------------
  140       continue
            ml = iwm(1)
            mu = iwm(2)
            ml3 = 3
            mband = ml + mu + 1
            meband = mband + ml
            lenp = meband*n
            do 150 i = 1, lenp
               wm(i+2) = 0.0d0
  150       continue
            call djac(tn,y,wm(ml3),meband,rpar,ipar)
            con = -hl0
            do 160 i = 1, lenp
               wm(i+2) = wm(i+2)*con
  160       continue
c        ...exit
            go to 220
c           if miter = 5, make mband calls to df to approximate j.
c           ----------------
  170       continue
            ml = iwm(1)
            mu = iwm(2)
            mband = ml + mu + 1
            mba = min(mband,n)
            meband = mband + ml
            meb1 = meband - 1
            srur = wm(1)
            fac = dvnrms(n,savf,ewt)
            r0 = 1000.0d0*abs(h)*uround*n*fac
            if (r0 .eq. 0.0d0) r0 = 1.0d0
            do 210 j = 1, mba
               do 180 i = j, n, mband
                  yi = y(i)
                  r = max(srur*abs(yi),r0*ewt(i))
                  y(i) = y(i) + r
  180          continue
               call df(tn,y,ftem,rpar,ipar)
               do 200 jj = j, n, mband
                  y(jj) = yh(jj,1)
                  yjj = y(jj)
                  r = max(srur*abs(yjj),r0*ewt(jj))
                  fac = -hl0/r
                  i1 = max(jj-mu,1)
                  i2 = min(jj+ml,n)
                  ii = jj*meb1 - ml + 2
                  do 190 i = i1, i2
                     wm(ii+i) = (ftem(i) - savf(i))*fac
  190             continue
  200          continue
  210       continue
            nfe = nfe + mba
  220    continue
c        add identity matrix.
c        -------------------------------------------------
         ii = mband + 2
         do 230 i = 1, n
            wm(ii) = wm(ii) + 1.0d0
            ii = ii + meband
  230    continue
c        do lu decomposition of p.
c        --------------------------------------------
         call dgbfa(wm(3),meband,n,ml,mu,iwm(21),ier)
  240 continue
      return
c     ----------------------- end of subroutine dpjac
c     -----------------------
      end
