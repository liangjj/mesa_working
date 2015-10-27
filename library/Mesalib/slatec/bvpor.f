*deck bvpor
      subroutine bvpor (y, nrowy, ncomp, xpts, nxpts, a, nrowa, alpha,
     +   nic, b, nrowb, beta, nfc, iflag, z, mxnon, p, ntp, ip, w, niv,
     +   yhp, u, v, coef, s, stowa, g, work, iwork, nfcc)
c***begin prologue  bvpor
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (bvpor-s, dbvpor-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c     input to bvpor    (items not defined in bvsup comments)
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
c             1 -- use gram-schmidt test and derkf
c             2 -- use gram-schmidt test and deabm
c
c     tol = tolerance for allowable error in orthogonalization test.
c
c     nps = 0 normalize particular solution to unit length at each
c             point of orthonormalization.
c         = 1 do not normalize particular solution.
c
c     ntp = must be .ge. nfc*(nfc+1)/2.
c
c
c     nfcc = 2*nfc for special treatment of a complex valued problem
c
c     icoco = 0 skip final computations (superposition coefficients
c               and ,hence, boundary problem solution)
c           = 1 calculate superposition coefficients and obtain
c               solution to the boundary value problem
c
c **********************************************************************
c     output from bvpor
c **********************************************************************
c
c     y(nrowy,nxpts) = solution at specified output points.
c
c     mxnon = number of orthonormalizations performed by bvpor.
c
c     z(mxnon+1) = locations of orthonormalizations performed by bvpor.
c
c     niv = number of independent vectors returned from mgsbv. normally
c        this parameter will be meaningful only when mgsbv returns with
c           mflag = 2.
c
c **********************************************************************
c
c     the following variables are in the argument list because of
c     variable dimensioning. in general, they contain no information of
c     use to the user.  the amount of storage set aside by the user must
c     be greater than or equal to that indicated by the dimension
c     statements.   for the disk storage mode, non = 0 and kpts = 1,
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
c     subroutines used by bvpor
c         lssuds -- solves an underdetermined system of linear
c                   equations.  this routine is used to get a full
c                   set of initial conditions for integration.
c                   called by bvpor
c
c         svecs -- obtains starting vectors for special treatment
c                  of complex valued problems , called by bvpor
c
c         rkfab -- routine which conducts integration using derkf or
c                   deabm
c
c         stway -- storage for backup capability, called by
c                   bvpor and reort
c
c         stor1 -- storage at output points, called by bvpor,
c                  rkfab, reort and stway.
c
c         sdot -- single precision vector inner product routine,
c                   called by bvpor, scoef, lssuds, mgsbv,
c                   bksol, reort and prvec.
c         ** note **
c         a considerable improvement in speed can be achieved if a
c         machine language version is used for sdot.
c
c         scoef -- computes the superposition constants from the
c                  boundary conditions at xfinal.
c
c         bksol -- solves an upper triangular set of linear equations.
c
c **********************************************************************
c
c***see also  bvsup
c***routines called  bksol, lssuds, rkfab, scoef, sdot, stor1, stway,
c                    svecs
c***common blocks    ml15to, ml18jr, ml8sz
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  bvpor
c
      dimension y(nrowy,*),a(nrowa,*),alpha(*),b(nrowb,*),
     1          beta(*),p(ntp,*),ip(nfcc,*),
     2          u(ncomp,nfc,*),v(ncomp,*),w(nfcc,*),
     3          coef(*),z(*),yhp(ncomp,*),xpts(*),s(*),
     4          work(*),iwork(*),stowa(*),g(*)
c
c **********************************************************************
c
      common /ml8sz/ c,xsav,igofx,inhomo,ivp,ncompd,nfcd
      common /ml15to/ px,pwcnd,tnd,x,xbeg,xend,xot,xop,info(15),istkop,
     1                knswot,kop,lotjp,mnswot,nswot
      common /ml18jr/ ae,re,tol,nxptsd,nicd,nopg,mxnond,ndisk,ntape,
     1                neq,indpvt,integ,nps,ntpd,neqivp,numort,nfccd,
     2                icoco
c
c **********************************************************************
c
c***first executable statement  bvpor
      nfcp1 = nfc + 1
      numort = 0
      c = 1.0
c
c **********************************************************************
c     calculate initial conditions which satisfy
c                   a*yh(xinitial)=0  and  a*yp(xinitial)=alpha.
c     when nfc .ne. nfcc lssuds defines values yhp in a matrix of size
c     (nfcc+1)*ncomp and ,hence, overflows the storage allocation into
c     the u array. however, this is okay since plenty of space is
c     available in u and it has not yet been used.
c
      ndw = nrowa * ncomp
      kws = ndw + nic + 1
      kwd = kws + nic
      kwt = kwd + nic
      kwc = kwt + nic
      iflag = 0
      call lssuds(a,yhp(1,nfcc+1),alpha,nic,ncomp,nrowa,yhp,ncomp,
     1            iflag,1,ira,0,work(1),work(ndw+1),iwork,work(kws),
     2            work(kwd),work(kwt),isflg,work(kwc))
      if (iflag .eq. 1) go to 3
      iflag=-4
      go to 250
    3 if (nfc .ne. nfcc) call svecs(ncomp,nfc,yhp,work,iwork,
     1                   inhomo,iflag)
      if (iflag .eq. 1)  go to 5
      iflag=-5
      go to 250
c
c **********************************************************************
c     determine the number of differential equations to be integrated,
c     initialize variables for auxiliary initial value problem and
c     store initial conditions.
c
    5 neq = ncomp * nfc
      if (inhomo .eq. 1)  neq = neq + ncomp
      ivp = 0
      if (neqivp .eq. 0)  go to 10
      ivp = neq
      neq = neq + neqivp
      nfcp2 = nfcp1
      if (inhomo .eq. 1)  nfcp2 = nfcp1 + 1
      do 7 k = 1,neqivp
    7 yhp(k,nfcp2) = alpha(nic+k)
   10 call stor1(u,yhp,v,yhp(1,nfcp1),0,ndisk,ntape)
c
c **********************************************************************
c     set up data for the orthonormalization testing procedure and
c     save initial conditions in case a restart is necessary.
c
      nswot=1
      knswot=0
      lotjp=1
      tnd=log10(10.*tol)
      pwcnd=log10(sqrt(tol))
      x=xbeg
      px=x
      xot=xend
      xop=x
      kop=1
      call stway(u,v,yhp,0,stowa)
c
c **********************************************************************
c ******** forward integration of all initial value equations **********
c **********************************************************************
c
      call rkfab(ncomp,xpts,nxpts,nfc,iflag,z,mxnon,p,ntp,ip,
     1            yhp,niv,u,v,w,s,stowa,g,work,iwork,nfcc)
      if (iflag .ne. 0  .or.  icoco .eq. 0)  go to 250
c
c **********************************************************************
c **************** backward sweep to obtain solution *******************
c **********************************************************************
c
c     calculate superposition coefficients at xfinal.
c
c   for the disk storage version, it is not necessary to read  u  and  v
c   at the last output point, since the local copy of each still exists.
c
      kod = 1
      if (ndisk .eq. 0)  kod = nxpts
      i1=1+nfcc*nfcc
      i2=i1+nfcc
      call scoef(u(1,1,kod),v(1,kod),ncomp,nrowb,nfc,nic,b,beta,coef,
     1           inhomo,re,ae,work,work(i1),work(i2),iwork,iflag,nfcc)
c
c **********************************************************************
c     calculate solution at output points by recurring backwards.
c     as we recur backwards from xfinal to xinitial we must calculate
c     new superposition coefficients each time we cross a point of
c     orthonormalization.
c
      k = numort
      ncomp2=ncomp/2
      ic=1
      if (nfc .ne. nfcc) ic=2
      do 200 j = 1,nxpts
      kpts = nxpts - j + 1
      kod = kpts
      if (ndisk .eq. 1)  kod = 1
  135 if (k .eq. 0)  go to 170
      if (xend.gt.xbeg .and. xpts(kpts).ge.z(k))  go to 170
      if (xend.lt.xbeg .and. xpts(kpts).le.z(k))  go to 170
      non = k
      if (ndisk .eq. 0)  go to 136
      non = 1
      backspace ntape
      read (ntape) (ip(i,1), i = 1,nfcc),(p(i,1), i = 1,ntp)
      backspace ntape
  136 if (inhomo .ne. 1)  go to 150
      if (ndisk .eq. 0)  go to 138
      backspace ntape
      read (ntape) (w(i,1), i = 1,nfcc)
      backspace ntape
  138 do 140 n = 1,nfcc
  140 coef(n) = coef(n) - w(n,non)
  150 call bksol(nfcc,p(1,non),coef)
      do 155 m = 1,nfcc
  155 work(m) = coef(m)
      do 160 m = 1,nfcc
      l = ip(m,non)
  160 coef(l) = work(m)
      k = k - 1
      go to 135
  170 if (ndisk .eq. 0)  go to 175
      backspace ntape
      read (ntape) (v(i,1), i = 1,ncomp),
     1             ((u(i,m,1), i = 1,ncomp), m = 1,nfc)
      backspace ntape
  175 do 180 n = 1,ncomp
  180 y(n,kpts) = v(n,kod) + sdot(nfc,u(n,1,kod),ncomp,coef,ic)
      if (nfc .eq. nfcc) go to 200
      do 190 n=1,ncomp2
      nn=ncomp2+n
      y(n,kpts)=y(n,kpts) - sdot(nfc,u(nn,1,kod),ncomp,coef(2),2)
  190 y(nn,kpts)=y(nn,kpts) + sdot(nfc,u(n,1,kod),ncomp,coef(2),2)
  200 continue
c
c **********************************************************************
c
  250 mxnon = numort
      return
      end
