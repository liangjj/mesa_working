*deck dbvpor
      subroutine dbvpor (y, nrowy, ncomp, xpts, nxpts, a, nrowa, alpha,
     +   nic, b, nrowb, beta, nfc, iflag, z, mxnon, p, ntp, ip, w, niv,
     +   yhp, u, v, coef, s, stowa, g, work, iwork, nfcc)
c***begin prologue  dbvpor
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (bvpor-s, dbvpor-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c     input to dbvpor    (items not defined in dbvsup comments)
c **********************************************************************
c
c     nopg = 0 -- orthonormalization points not pre-assigned
c          = 1 -- orthonormalization points pre-assigned
c
c     mxnon = maximum number of orthogonalizations allowed.
c
c     ndisk = 0 -- in-core storage
c           = 1 -- disk storage.  value of ntape in data statement
c                  is set to 13.  if another value is desired,
c                  the data statement must be changed.
c
c     integ = type of integrator and associated test to be used
c             to determine when to orthonormalize.
c
c             1 -- use gram-schmidt test and dderkf
c             2 -- use gram-schmidt test and ddeabm
c
c     tol = tolerance for allowable error in orthogonalization test.
c
c     nps = 0 normalize particular solution to unit length at each
c             point of orthonormalization.
c         = 1 do not normalize particular solution.
c
c     ntp = must be .ge. nfc*(nfc+1)/2.
c
c     nfcc = 2*nfc for special treatment of a complex*16 valued problem
c
c     icoco = 0 skip final computations (superposition coefficients
c               and, hence, boundary problem solution)
c           = 1 calculate superposition coefficients and obtain
c               solution to the boundary value problem
c
c **********************************************************************
c     output from dbvpor
c **********************************************************************
c
c     y(nrowy,nxpts) = solution at specified output points.
c
c     mxnon = number of orthonormalizations performed by dbvpor.
c
c     z(mxnon+1) = locations of orthonormalizations performed by dbvpor.
c
c     niv = number of independent vectors returned from dmgsbv. normally
c           this parameter will be meaningful only when dmgsbv returns
c           with mflag = 2.
c
c **********************************************************************
c
c     the following variables are in the argument list because of
c     variable dimensioning.  in general, they contain no information of
c     use to the user.  the amount of storage set aside by the user must
c     be greater than or equal to that indicated by the dimension
c     statements.  for the disk storage mode, non = 0 and kpts = 1,
c     while for the in-core storage mode, non = mxnon and kpts = nxpts.
c
c     p(ntp,non+1)
c     ip(nfcc,non+1)
c     yhp(ncomp,nfc+1)  plus an additional column of the length  neqivp
c     u(ncomp,nfc,kpts)
c     v(ncomp,kpts)
c     w(nfcc,non+1)
c     coef(nfcc)
c     s(nfc+1)
c     stowa(ncomp*(nfc+1)+neqivp+1)
c     g(ncomp)
c     work(kkkws)
c     iwork(llliws)
c
c **********************************************************************
c     subroutines used by dbvpor
c         dlssud -- solves an underdetermined system of linear
c                   equations.  this routine is used to get a full
c                   set of initial conditions for integration.
c                   called by dbvpor.
c
c         dvecs -- obtains starting vectors for special treatment
c                   of complex*16 valued problems, called by dbvpor.
c
c         drkfab -- routine which conducts integration using dderkf or
c                   ddeabm.
c
c         dstway -- storage for backup capability, called by
c                   dbvpor and dreort.
c
c         dstor1 -- storage at output points, called by dbvpor,
c                   drkfab, dreort and dstway.
c
c         ddot -- single precision vector inner product routine,
c                   called by dbvpor, dcoef, dlssud, dmgsbv,
c                   dbksol, dreort and dprvec.
c         ** note **
c         a considerable improvement in speed can be achieved if a
c         machine language version is used for ddot.
c
c         dcoef -- computes the superposition constants from the
c                   boundary conditions at xfinal.
c
c         dbksol -- solves an upper triangular set of linear equations.
c
c **********************************************************************
c
c***see also  dbvsup
c***routines called  dbksol, dcoef, ddot, dlssud, drkfab, dstor1,
c                    dstway, dvecs
c***common blocks    dml15t, dml18j, dml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dbvpor
c
      double precision ddot
      integer i, i1, i2, ic, icoco, iflag, igofx, indpvt, info, inhomo,
     1     integ, ira, isflg, istkop, ivp, j,
     2     k, knswot, kod, kop, kpts, kwc, kwd, kws, kwt, l, lotjp, m,
     3     mnswot, mxnon, mxnond, n, ncomp, ncomp2, ncompd, ndisk, ndw,
     4     neq, neqivp, nfc, nfcc, nfccd, nfcd, nfcp1, nfcp2, nic,
     5     nicd, niv, nn, non, nopg, nps, nrowa, nrowb, nrowy, nswot,
     6     ntape, ntp, ntpd, numort, nxpts, nxptsd,
     7     ip(nfcc,*), iwork(*)
      double precision a(nrowa,*), ae, alpha(*), b(nrowb,*),
     1     beta(*), c, coef(*), g(*), p(ntp,*), pwcnd, px,
     2     re, s(*), stowa(*), tnd, tol, u(ncomp,nfc,*),
     3     v(ncomp,*), w(nfcc,*), work(*), x, xbeg, xend, xop,
     4     xot, xpts(*), xsav, y(nrowy,*), yhp(ncomp,*),
     5     z(*)
c
c     ******************************************************************
c
      common /dml8sz/ c,xsav,igofx,inhomo,ivp,ncompd,nfcd
      common /dml15t/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /dml18j/ ae,re,tol,nxptsd,nicd,nopg,mxnond,ndisk,ntape,
     1                neq,indpvt,integ,nps,ntpd,neqivp,numort,nfccd,
     2                icoco
c
c      *****************************************************************
c
c***first executable statement  dbvpor
      nfcp1 = nfc + 1
      numort = 0
      c = 1.0d0
c
c     ******************************************************************
c         calculate initial conditions which satisfy
c                       a*yh(xinitial)=0  and  a*yp(xinitial)=alpha.
c         when nfc .ne. nfcc dlssud defines values yhp in a matrix of
c         size (nfcc+1)*ncomp and ,hence, overflows the storage
c         allocation into the u array. however, this is okay since
c         plenty of space is available in u and it has not yet been
c         used.
c
      ndw = nrowa*ncomp
      kws = ndw + nic + 1
      kwd = kws + nic
      kwt = kwd + nic
      kwc = kwt + nic
      iflag = 0
      call dlssud(a,yhp(1,nfcc+1),alpha,nic,ncomp,nrowa,yhp,ncomp,iflag,
     1            1,ira,0,work(1),work(ndw+1),iwork,work(kws),work(kwd),
     2            work(kwt),isflg,work(kwc))
      if (iflag .eq. 1) go to 10
         iflag = -4
      go to 200
   10 continue
         if (nfc .ne. nfcc)
     1      call dvecs(ncomp,nfc,yhp,work,iwork,inhomo,iflag)
         if (iflag .eq. 1) go to 20
            iflag = -5
         go to 190
   20    continue
c
c           ************************************************************
c               determine the number of differential equations to be
c               integrated, initialize variables for auxiliary initial
c               value problem and store initial conditions.
c
            neq = ncomp*nfc
            if (inhomo .eq. 1) neq = neq + ncomp
            ivp = 0
            if (neqivp .eq. 0) go to 40
               ivp = neq
               neq = neq + neqivp
               nfcp2 = nfcp1
               if (inhomo .eq. 1) nfcp2 = nfcp1 + 1
               do 30 k = 1, neqivp
                  yhp(k,nfcp2) = alpha(nic+k)
   30          continue
   40       continue
            call dstor1(u,yhp,v,yhp(1,nfcp1),0,ndisk,ntape)
c
c           ************************************************************
c               set up data for the orthonormalization testing procedure
c               and save initial conditions in case a restart is
c               necessary.
c
            nswot = 1
            knswot = 0
            lotjp = 1
            tnd = log10(10.0d0*tol)
            pwcnd = log10(sqrt(tol))
            x = xbeg
            px = x
            xot = xend
            xop = x
            kop = 1
            call dstway(u,v,yhp,0,stowa)
c
c           ************************************************************
c           ******** forward integration of all initial value equations
c           **********
c           ************************************************************
c
            call drkfab(ncomp,xpts,nxpts,nfc,iflag,z,mxnon,p,ntp,ip,yhp,
     1                  niv,u,v,w,s,stowa,g,work,iwork,nfcc)
            if (iflag .ne. 0 .or. icoco .eq. 0) go to 180
c
c              *********************************************************
c              **************** backward sweep to obtain solution
c              *******************
c              *********************************************************
c
c                  calculate superposition coefficients at xfinal.
c
c                for the disk storage version, it is not necessary to
c                read  u  and  v at the last output point, since the
c                local copy of each still exists.
c
               kod = 1
               if (ndisk .eq. 0) kod = nxpts
               i1 = 1 + nfcc*nfcc
               i2 = i1 + nfcc
               call dcoef(u(1,1,kod),v(1,kod),ncomp,nrowb,nfc,nic,b,
     1                     beta,coef,inhomo,re,ae,work,work(i1),
     2                     work(i2),iwork,iflag,nfcc)
c
c              *********************************************************
c                  calculate solution at output points by recurring
c                  backwards.  as we recur backwards from xfinal to
c                  xinitial we must calculate new superposition
c                  coefficients each time we cross a point of
c                  orthonormalization.
c
               k = numort
               ncomp2 = ncomp/2
               ic = 1
               if (nfc .ne. nfcc) ic = 2
               do 170 j = 1, nxpts
                  kpts = nxpts - j + 1
                  kod = kpts
                  if (ndisk .eq. 1) kod = 1
   50             continue
c                 ...exit
                     if (k .eq. 0) go to 120
c                 ...exit
                     if (xend .gt. xbeg .and. xpts(kpts) .ge. z(k))
     1                  go to 120
c                 ...exit
                     if (xend .lt. xbeg .and. xpts(kpts) .le. z(k))
     1                  go to 120
                     non = k
                     if (ndisk .eq. 0) go to 60
                        non = 1
                        backspace ntape
                        read (ntape)
     1                       (ip(i,1), i = 1, nfcc),(p(i,1), i = 1, ntp)
                        backspace ntape
   60                continue
                     if (inhomo .ne. 1) go to 90
                        if (ndisk .eq. 0) go to 70
                           backspace ntape
                           read (ntape) (w(i,1), i = 1, nfcc)
                           backspace ntape
   70                   continue
                        do 80 n = 1, nfcc
                           coef(n) = coef(n) - w(n,non)
   80                   continue
   90                continue
                     call dbksol(nfcc,p(1,non),coef)
                     do 100 m = 1, nfcc
                        work(m) = coef(m)
  100                continue
                     do 110 m = 1, nfcc
                        l = ip(m,non)
                        coef(l) = work(m)
  110                continue
                     k = k - 1
                  go to 50
  120             continue
                  if (ndisk .eq. 0) go to 130
                     backspace ntape
                     read (ntape)
     1                    (v(i,1), i = 1, ncomp),
     2                    ((u(i,m,1), i = 1, ncomp), m = 1, nfc)
                     backspace ntape
  130             continue
                  do 140 n = 1, ncomp
                     y(n,kpts) = v(n,kod)
     1                           + ddot(nfc,u(n,1,kod),ncomp,coef,ic)
  140             continue
                  if (nfc .eq. nfcc) go to 160
                     do 150 n = 1, ncomp2
                        nn = ncomp2 + n
                        y(n,kpts) = y(n,kpts)
     1                              - ddot(nfc,u(nn,1,kod),ncomp,
     2                                     coef(2),2)
                        y(nn,kpts) = y(nn,kpts)
     1                               + ddot(nfc,u(n,1,kod),ncomp,
     2                                      coef(2),2)
  150                continue
  160             continue
  170          continue
  180       continue
  190    continue
  200 continue
c
c     ******************************************************************
c
      mxnon = numort
      return
      end
