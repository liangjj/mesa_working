*deck polyvl
      subroutine polyvl (nder, xx, yfit, yp, n, x, c, work, ierr)
c***begin prologue  polyvl
c***purpose  calculate the value of a polynomial and its first nder
c            derivatives where the polynomial was produced by a previous
c            call to polint.
c***library   slatec
c***category  e3
c***type      single precision (polyvl-s, dpolvl-d)
c***keywords  polynomial evaluation
c***author  huddleston, r. e., (snll)
c***description
c
c     written by robert e. huddleston, sandia laboratories, livermore
c
c     abstract -
c        subroutine polyvl calculates the value of the polynomial and
c     its first nder derivatives where the polynomial was produced by
c     a previous call to polint.
c        the variable n and the arrays x and c must not be altered
c     between the call to polint and the call to polyvl.
c
c     ******  dimensioning information *******
c
c     yp   must be dimensioned by at least nder
c     x    must be dimensioned by at least n (see the abstract )
c     c    must be dimensioned by at least n (see the abstract )
c     work must be dimensioned by at least 2*n if nder is .gt. 0.
c
c     *** note ***
c       if nder=0, neither yp nor work need to be dimensioned variables.
c       if nder=1, yp does not need to be a dimensioned variable.
c
c
c     *****  input parameters
c
c     nder - the number of derivatives to be evaluated
c
c     xx   - the argument at which the polynomial and its derivatives
c            are to be evaluated.
c
c     n    - *****
c            *       n, x, and c must not be altered between the call
c     x    - *       to polint and the call to polyvl.
c     c    - *****
c
c
c     *****  output parameters
c
c     yfit - the value of the polynomial at xx
c
c     yp   - the derivatives of the polynomial at xx.  the derivative of
c            order j at xx is stored in  yp(j) , j = 1,...,nder.
c
c     ierr - output error flag with the following possible values.
c          = 1  indicates normal execution
c
c     ***** storage parameters
c
c     work  = this is an array to provide internal working storage for
c             polyvl.  it must be dimensioned by at least 2*n if nder is
c             .gt. 0.  if nder=0, work does not need to be a dimensioned
c             variable.
c
c***references  l. f. shampine, s. m. davenport and r. e. huddleston,
c                 curve fitting by polynomials in one variable, report
c                 sla-74-0270, sandia laboratories, june 1974.
c***routines called  (none)
c***revision history  (yymmdd)
c   740601  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  polyvl
      dimension  yp(*),x(*),c(*),work(*)
c***first executable statement  polyvl
      ierr=1
         if (nder.gt.0) go to 10020
c
c     *****   coding for the case nder = 0
c
      pione=1.0
      pone=c(1)
      yfit=pone
      if (n.eq.1) return
      do 10010 k=2,n
      pitwo=(xx-x(k-1))*pione
      pione=pitwo
      ptwo=pone+pitwo*c(k)
      pone=ptwo
10010 continue
      yfit=ptwo
      return
c
c     *****   end of nder = 0 case
c
10020 continue
         if (n.gt.1) go to 10040
      yfit=c(1)
c
c     *****  coding for the case  n=1 and nder .gt. 0
c
      do 10030 k=1,nder
      yp(k)=0.0
10030 continue
      return
c
c     *****  end of the case  n = 1 and  nder .gt. 0
c
10040 continue
         if (nder.lt.n) go to 10050
c
c     *****  set flags for number of derivatives and for derivatives
c            in excess of the degree (n-1) of the polynomial.
c
      izero=1
      ndr=n-1
         go to 10060
10050 continue
      izero=0
      ndr=nder
10060 continue
      m=ndr+1
      mm=m
c
c     *****  start of the case nder .gt. 0  and n .gt. 1
c     *****  the polynomial and its derivatives will be evaluated at xx
c
      do 10070 k=1,ndr
      yp(k)=c(k+1)
10070 continue
c
c     *****  the following section of code is easier to read if one
c            breaks work into two arrays w and v. the code would then
c            read
c                w(1) = 1.
c                pone = c(1)
c               *do   k = 2,n
c               *   v(k-1) =  xx - x(k-1)
c               *   w(k)   =  v(k-1)*w(k-1)
c               *   ptwo   =  pone + w(k)*c(k)
c               *   pone   =  pwo
c
c               yfit = ptwo
c
      work(1)=1.0
      pone=c(1)
      do 10080 k=2,n
      km1=k-1
      npkm1=n+k-1
      work(npkm1)=xx-x(km1)
      work(k)=work(npkm1)*work(km1)
      ptwo=pone+work(k)*c(k)
      pone=ptwo
10080 continue
      yfit=ptwo
c
c     ** at this point the polynomial has been evaluated and information
c        for the derivative evaluations have been stored in the array
c        work
         if (n.eq.2) go to 10110
      if (m.eq.n) mm=ndr
c
c     ***** evaluate the derivatives at xx
c
c                  ******  do k=2,mm   (for most cases, mm = nder + 1)
c                  *  ******  do i=2,n-k+1
c                  *  *       w(i) = v(k-2+i)*w(i-1) + w(i)
c                  *  *       yp(k-1) = yp(k-1) + w(i)*c(k-1+i)
c                  ******  continue
c
      do 10090 k=2,mm
      nmkp1=n-k+1
      km1=k-1
      km2pn=k-2+n
      do 10090 i=2,nmkp1
      km2pni=km2pn+i
      im1=i-1
      km1pi=km1+i
      work(i)=work(km2pni)*work(im1)+work(i)
      yp(km1)=yp(km1)+work(i)*c(km1pi)
10090 continue
         if (ndr.eq.1) go to 10110
      fac=1.0
      do 10100 k=2,ndr
      xk=k
      fac=xk*fac
      yp(k)=fac*yp(k)
10100 continue
c
c     ***** end of derivative evaluations
c
10110 continue
      if (izero.eq.0) return
c
c     *****  set excess derivatives to zero.
c
      do 10120 k=n,nder
      yp(k)=0.0
10120 continue
      return
      end
