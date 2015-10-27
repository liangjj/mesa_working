*deck dcoef
      subroutine dcoef (yh, yp, ncomp, nrowb, nfc, nic, b, beta, coef,
     +   inhomo, re, ae, by, cvec, work, iwork, iflag, nfcc)
c***begin prologue  dcoef
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (scoef-s, dcoef-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c input to dcoef
c **********************************************************************
c
c     yh = matrix of homogeneous solutions.
c     yp = vector containing particular solution.
c     ncomp = number of components per solution vector.
c     nrowb = first dimension of b in calling program.
c     nfc = number of base solution vectors.
c     nfcc = 2*nfc for the special treatment of complex*16 valued
c            equations. otherwise, nfcc=nfc.
c     nic = number of specified initial conditions.
c     b = boundary condition matrix at x = xfinal.
c     beta = vector of nonhomogeneous boundary conditions at x = xfinal.
c              1 - nonzero particular solution
c     inhomo = 2 - zero particular solution
c              3 - eigenvalue problem
c     re = relative error tolerance.
c     ae = absolute error tolerance.
c     by = storage space for the matrix  b*yh
c     cvec = storage space for the vector  beta-b*yp
c     work = double precision array of internal storage. dimension must
c     be ge
c            nfcc*(nfcc+4)
c     iwork = integer array of internal storage. dimension must be ge
c             3+nfcc
c
c **********************************************************************
c output from dcoef
c **********************************************************************
c
c     coef = array containing superposition constants.
c     iflag = indicator of success from dsuds in solving the
c             boundary equations.
c           = 0 boundary equations are solved.
c           = 1 boundary equations appear to have many solutions.
c           = 2 boundary equations appear to be inconsistent.
c           = 3 for this value of an eigenparameter, the boundary
c               equations have only the zero solution.
c
c **********************************************************************
c
c     subroutine dcoef solves for the superposition constants from the
c     linear equations defined by the boundary conditions at x = xfinal.
c
c                          b*yp + b*yh*coef = beta
c
c **********************************************************************
c
c***see also  dbvsup
c***routines called  ddot, dsuds, xgetf, xsetf
c***common blocks    dml5mc
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   890921  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  dcoef
c
      double precision ddot
      integer i, iflag, inhomo, iwork(*), j, k, kflag, ki, l, lpar,
     1     mlso, ncomp, ncomp2, nf, nfc, nfcc, nfccm1, nic,
     2     nrowb
      double precision ae, b(nrowb,*), bbn, beta(*), bn, brn,
     1     by(nfcc,*), bykl, bys, coef(*), cons, cvec(*), eps,
     2     fouru, gam, re, sqovfl, sru, twou, un, uro, work(*),
     3     yh(ncomp,*), yp(*), ypn
c
      common /dml5mc/ uro,sru,eps,sqovfl,twou,fouru,lpar
c***first executable statement  dcoef
c
c     set up matrix  b*yh  and vector  beta - b*yp
c
      ncomp2 = ncomp/2
      do 80 k = 1, nfcc
         do 10 j = 1, nfc
            l = j
            if (nfc .ne. nfcc) l = 2*j - 1
            by(k,l) = ddot(ncomp,b(k,1),nrowb,yh(1,j),1)
   10    continue
         if (nfc .eq. nfcc) go to 30
            do 20 j = 1, nfc
               l = 2*j
               bykl = ddot(ncomp2,b(k,1),nrowb,yh(ncomp2+1,j),1)
               by(k,l) = ddot(ncomp2,b(k,ncomp2+1),nrowb,yh(1,j),1)
     1                   - bykl
   20       continue
   30    continue
         go to (40,50,60), inhomo
c        case 1
   40    continue
            cvec(k) = beta(k) - ddot(ncomp,b(k,1),nrowb,yp,1)
         go to 70
c        case 2
   50    continue
            cvec(k) = beta(k)
         go to 70
c        case 3
   60    continue
            cvec(k) = 0.0d0
   70    continue
   80 continue
      cons = abs(cvec(1))
      bys = abs(by(1,1))
c
c     ******************************************************************
c         solve linear system
c
      iflag = 0
      mlso = 0
      if (inhomo .eq. 3) mlso = 1
      kflag = 0.5d0 * log10(eps)
      call xgetf(nf)
      call xsetf(0)
   90 continue
         call dsuds(by,coef,cvec,nfcc,nfcc,nfcc,kflag,mlso,work,iwork)
         if (kflag .ne. 3) go to 100
         kflag = 1
         iflag = 1
      go to 90
  100 continue
      if (kflag .eq. 4) iflag = 2
      call xsetf(nf)
      if (nfcc .eq. 1) go to 180
         if (inhomo .ne. 3) go to 170
            if (iwork(1) .lt. nfcc) go to 140
               iflag = 3
               do 110 k = 1, nfcc
                  coef(k) = 0.0d0
  110          continue
               coef(nfcc) = 1.0d0
               nfccm1 = nfcc - 1
               do 130 k = 1, nfccm1
                  j = nfcc - k
                  l = nfcc - j + 1
                  gam = ddot(l,by(j,j),nfcc,coef(j),1)/(work(j)*by(j,j))
                  do 120 i = j, nfcc
                     coef(i) = coef(i) + gam*by(j,i)
  120             continue
  130          continue
            go to 160
  140       continue
               do 150 k = 1, nfcc
                  ki = 4*nfcc + k
                  coef(k) = work(ki)
  150          continue
  160       continue
  170    continue
      go to 220
  180 continue
c
c        ***************************************************************
c            testing for existence and uniqueness of boundary-value
c            problem solution in a scalar case
c
         bn = 0.0d0
         un = 0.0d0
         ypn = 0.0d0
         do 190 k = 1, ncomp
            un = max(un,abs(yh(k,1)))
            ypn = max(ypn,abs(yp(k)))
            bn = max(bn,abs(b(1,k)))
  190    continue
         bbn = max(bn,abs(beta(1)))
         if (bys .gt. 10.0d0*(re*un + ae)*bn) go to 200
            brn = bbn/bn*bys
            if (cons .ge. 0.1d0*brn .and. cons .le. 10.0d0*brn)
     1         iflag = 1
            if (cons .gt. 10.0d0*brn) iflag = 2
            if (cons .le. re*abs(beta(1)) + ae + (re*ypn + ae)*bn)
     1         iflag = 1
            if (inhomo .eq. 3) coef(1) = 1.0d0
         go to 210
  200    continue
         if (inhomo .ne. 3) go to 210
            iflag = 3
            coef(1) = 1.0d0
  210    continue
  220 continue
      return
      end
