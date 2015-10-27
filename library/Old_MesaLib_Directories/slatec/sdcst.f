*deck sdcst
      subroutine sdcst (maxord, mint, iswflg, el, tq)
c***begin prologue  sdcst
c***subsidiary
c***purpose  sdcst sets coefficients used by the core integrator sdstp.
c***library   slatec (sdrive)
c***type      single precision (sdcst-s, ddcst-d, cdcst-c)
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***description
c
c  sdcst is called by sdntl.  the array el determines the basic method.
c  the array tq is involved in adjusting the step size in relation
c  to truncation error.  el and tq depend upon mint, and are calculated
c  for orders 1 to maxord(.le. 12).  for each order nq, the coefficients
c  el are calculated from the generating polynomial:
c    l(t) = el(1,nq) + el(2,nq)*t + ... + el(nq+1,nq)*t**nq.
c  for the implicit adams methods, l(t) is given by
c    dl/dt = (1+t)*(2+t)* ... *(nq-1+t)/k,   l(-1) = 0,
c    where      k = factorial(nq-1).
c  for the gear methods,
c    l(t) = (1+t)*(2+t)* ... *(nq+t)/k,
c    where      k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c  for each order nq, there are three components of tq.
c***routines called  (none)
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  sdcst
      real el(13,12), factrl(12), gamma(14), sum, tq(3,12)
      integer i, iswflg, j, maxord, mint, mxrd
c***first executable statement  sdcst
      factrl(1) = 1.e0
      do 10 i = 2,maxord
 10     factrl(i) = i*factrl(i-1)
c                                             compute adams coefficients
      if (mint .eq. 1) then
        gamma(1) = 1.e0
        do 40 i = 1,maxord+1
          sum = 0.e0
          do 30 j = 1,i
 30         sum = sum - gamma(j)/(i-j+2)
 40       gamma(i+1) = sum
        el(1,1) = 1.e0
        el(2,1) = 1.e0
        el(2,2) = 1.e0
        el(3,2) = 1.e0
        do 60 j = 3,maxord
          el(2,j) = factrl(j-1)
          do 50 i = 3,j
 50         el(i,j) = (j-1)*el(i,j-1) + el(i-1,j-1)
 60       el(j+1,j) = 1.e0
        do 80 j = 2,maxord
          el(1,j) = el(1,j-1) + gamma(j)
          el(2,j) = 1.e0
          do 80 i = 3,j+1
 80         el(i,j) = el(i,j)/((i-1)*factrl(j-1))
        do 100 j = 1,maxord
          tq(1,j) = -1.e0/(factrl(j)*gamma(j))
          tq(2,j) = -1.e0/gamma(j+1)
 100      tq(3,j) = -1.e0/gamma(j+2)
c                                              compute gear coefficients
      else if (mint .eq. 2) then
        el(1,1) = 1.e0
        el(2,1) = 1.e0
        do 130 j = 2,maxord
          el(1,j) = factrl(j)
          do 120 i = 2,j
 120        el(i,j) = j*el(i,j-1) + el(i-1,j-1)
 130      el(j+1,j) = 1.e0
        sum = 1.e0
        do 150 j = 2,maxord
          sum = sum + 1.e0/j
          do 150 i = 1,j+1
 150        el(i,j) = el(i,j)/(factrl(j)*sum)
        do 170 j = 1,maxord
          if (j .gt. 1) tq(1,j) = 1.e0/factrl(j-1)
          tq(2,j) = (j+1)/el(1,j)
 170      tq(3,j) = (j+2)/el(1,j)
      end if
c                          compute constants used in the stiffness test.
c                          these are the ratio of tq(2,nq) for the gear
c                          methods to those for the adams methods.
      if (iswflg .eq. 3) then
        mxrd = min(maxord, 5)
        if (mint .eq. 2) then
          gamma(1) = 1.e0
          do 190 i = 1,mxrd
            sum = 0.e0
            do 180 j = 1,i
 180          sum = sum - gamma(j)/(i-j+2)
 190        gamma(i+1) = sum
        end if
        sum = 1.e0
        do 200 i = 2,mxrd
          sum = sum + 1.e0/i
 200      el(1+i,1) = -(i+1)*sum*gamma(i+1)
      end if
      return
      end
