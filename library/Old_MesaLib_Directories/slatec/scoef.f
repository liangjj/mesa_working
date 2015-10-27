*deck scoef
      subroutine scoef (yh, yp, ncomp, nrowb, nfc, nic, b, beta, coef,
     +   inhomo, re, ae, by, cvec, work, iwork, iflag, nfcc)
c***begin prologue  scoef
c***subsidiary
c***purpose  subsidiary to bvsup
c***library   slatec
c***type      single precision (scoef-s, dcoef-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c input to scoef
c **********************************************************************
c
c     yh = matrix of homogeneous solutions.
c     yp = vector containing particular solution.
c     ncomp = number of components per solution vector.
c     nrowb = first dimension of b in calling program.
c     nfc = number of base solution vectors.
c     nfcc = 2*nfc for the special treatment of complex valued
c            equations. otherwise, nfcc=nfc.
c     nic = number of specified initial conditions.
c     b = boundary condition matrix at x = xfinal.
c     beta = vector of nonhomogeneous boundary conditions at x = xfinal.
c              1 - nonzero particular solution
c     inhomo = 2 - zero particular solution
c              3 - eigenvalue problem
c     re = relative error tolerance
c     ae = absolute error tolerance
c     by = storage space for the matrix  b*yh
c     cvec = storage space for the vector  beta-b*yp
c     work = real array of internal storage. dimension must be .ge.
c            nfcc*(nfcc+4)
c     iwork = integer array of internal storage. dimension must be .ge.
c             3+nfcc
c
c **********************************************************************
c output from scoef
c **********************************************************************
c
c     coef = array containing superposition constants.
c     iflag = indicator of success from suds in solving the
c             boundary equations
c           = 0 boundary equations are solved
c           = 1 boundary equations appear to have many solutions
c           = 2 boundary equations appear to be inconsistent
c           = 3 for this value of an eigenparameter, the boundary
c               equations have only the zero solution.
c
c **********************************************************************
c
c     subroutine scoef solves for the superposition constants from the
c     linear equations defined by the boundary conditions at x = xfinal.
c
c                          b*yp + b*yh*coef = beta
c
c **********************************************************************
c
c***see also  bvsup
c***routines called  sdot, suds, xgetf, xsetf
c***common blocks    ml5mco
c***revision history  (yymmdd)
c   750601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c***end prologue  scoef
c
      dimension yh(ncomp,*),yp(*),b(nrowb,*),beta(*),
     1          coef(*),by(nfcc,*),cvec(*),work(*),iwork(*)
c
      common /ml5mco/ uro,sru,eps,sqovfl,twou,fouru,lpar
c
c     set up matrix  b*yh  and vector  beta - b*yp
c
c***first executable statement  scoef
      ncomp2=ncomp/2
      do 7 k = 1,nfcc
      do 1 j = 1,nfc
      l=j
      if (nfc .ne. nfcc) l=2*j-1
    1 by(k,l) = sdot(ncomp,b(k,1),nrowb,yh(1,j),1)
      if (nfc .eq. nfcc) go to 3
      do 2 j=1,nfc
      l=2*j
      bykl=sdot(ncomp2,b(k,1),nrowb,yh(ncomp2+1,j),1)
      by(k,l)=sdot(ncomp2,b(k,ncomp2+1),nrowb,yh(1,j),1) - bykl
    2 continue
    3 go to (4,5,6), inhomo
c     case 1
    4 cvec(k) = beta(k) - sdot(ncomp,b(k,1),nrowb,yp,1)
      go to 7
c     case 2
    5 cvec(k) = beta(k)
      go to 7
c     case 3
    6 cvec(k) = 0.
    7 continue
      cons=abs(cvec(1))
      bys=abs(by(1,1))
c
c **********************************************************************
c     solve linear system
c
      iflag=0
      mlso=0
      if (inhomo .eq. 3) mlso=1
      kflag = 0.5 * log10(eps)
      call xgetf(nf)
      call xsetf(0)
   10 call suds(by,coef,cvec,nfcc,nfcc,nfcc,kflag,mlso,work,iwork)
      if (kflag .ne. 3) go to 13
      kflag=1
      iflag=1
      go to 10
   13 if (kflag .eq. 4) iflag=2
      call xsetf(nf)
      if (nfcc .eq. 1) go to 25
      if (inhomo .ne. 3) return
      if (iwork(1) .lt. nfcc) go to 17
      iflag=3
      do 14 k=1,nfcc
   14 coef(k)=0.
      coef(nfcc)=1.
      nfccm1=nfcc-1
      do 15 k=1,nfccm1
      j=nfcc-k
      l=nfcc-j+1
      gam=sdot(l,by(j,j),nfcc,coef(j),1)/(work(j)*by(j,j))
      do 15 i=j,nfcc
   15 coef(i)=coef(i)+gam*by(j,i)
      return
   17 do 20 k=1,nfcc
      ki=4*nfcc+k
   20 coef(k)=work(ki)
      return
c
c **********************************************************************
c     testing for existence and uniqueness of boundary-value problem
c     solution in a scalar case
c
   25 bn = 0.
      un = 0.
      ypn=0.
      do 30 k = 1,ncomp
      un = max(un,abs(yh(k,1)))
      ypn=max(ypn,abs(yp(k)))
   30 bn = max(bn,abs(b(1,k)))
      bbn = max(bn,abs(beta(1)))
      if (bys .gt. 10.*(re*un + ae)*bn)  go to 35
      brn = bbn / bn * bys
      if (cons .ge. 0.1*brn  .and.  cons .le. 10.*brn) iflag=1
      if (cons .gt. 10.*brn) iflag=2
      if (cons  .le.  re*abs(beta(1))+ae + (re*ypn+ae)*bn) iflag=1
      if (inhomo .eq. 3) coef(1)=1.
      return
   35 if (inhomo .ne. 3) return
      iflag=3
      coef(1)=1.
      return
      end
