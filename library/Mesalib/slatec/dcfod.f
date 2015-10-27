*deck dcfod
      subroutine dcfod (meth, elco, tesco)
c***begin prologue  dcfod
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (cfod-s, dcfod-d)
c***author  (unknown)
c***description
c
c   dcfod defines coefficients needed in the integrator package ddebdf
c
c***see also  ddebdf
c***routines called  (none)
c***revision history  (yymmdd)
c   820301  date written
c   890911  removed unnecessary intrinsics.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c***end prologue  dcfod
c
c
      integer i, ib, meth, nq, nqm1, nqp1
      double precision agamq, elco, fnq, fnqm1, pc, pint, ragq,
     1      rq1fac, rqfac, tesco, tsign, xpin
      dimension elco(13,12),tesco(3,12)
c     ------------------------------------------------------------------
c      dcfod  is called by the integrator routine to set coefficients
c      needed there.  the coefficients for the current method, as
c      given by the value of meth, are set for all orders and saved.
c      the maximum order assumed here is 12 if meth = 1 and 5 if meth =
c      2.  (a smaller value of the maximum order is also allowed.)
c      dcfod  is called once at the beginning of the problem,
c      and is not called again unless and until meth is changed.
c
c      the elco array contains the basic method coefficients.
c      the coefficients el(i), 1 .le. i .le. nq+1, for the method of
c      order nq are stored in elco(i,nq).  they are given by a
c      generating polynomial, i.e.,
c          l(x) = el(1) + el(2)*x + ... + el(nq+1)*x**nq.
c      for the implicit adams methods, l(x) is given by
c          dl/dx = (x+1)*(x+2)*...*(x+nq-1)/factorial(nq-1),    l(-1) =
c      0.  for the bdf methods, l(x) is given by
c          l(x) = (x+1)*(x+2)* ... *(x+nq)/k,
c      where         k = factorial(nq)*(1 + 1/2 + ... + 1/nq).
c
c      the tesco array contains test constants used for the
c      local error test and the selection of step size and/or order.
c      at order nq, tesco(k,nq) is used for the selection of step
c      size at order nq - 1 if k = 1, at order nq if k = 2, and at order
c      nq + 1 if k = 3.
c     ------------------------------------------------------------------
      dimension pc(12)
c
c***first executable statement  dcfod
      go to (10,60), meth
c
   10 continue
         elco(1,1) = 1.0d0
         elco(2,1) = 1.0d0
         tesco(1,1) = 0.0d0
         tesco(2,1) = 2.0d0
         tesco(1,2) = 1.0d0
         tesco(3,12) = 0.0d0
         pc(1) = 1.0d0
         rqfac = 1.0d0
         do 50 nq = 2, 12
c           ------------------------------------------------------------
c            the pc array will contain the coefficients of the
c                polynomial p(x) = (x+1)*(x+2)*...*(x+nq-1).
c            initially, p(x) = 1.
c           ------------------------------------------------------------
            rq1fac = rqfac
            rqfac = rqfac/nq
            nqm1 = nq - 1
            fnqm1 = nqm1
            nqp1 = nq + 1
c           form coefficients of p(x)*(x+nq-1).
c           ----------------------------------
            pc(nq) = 0.0d0
            do 20 ib = 1, nqm1
               i = nqp1 - ib
               pc(i) = pc(i-1) + fnqm1*pc(i)
   20       continue
            pc(1) = fnqm1*pc(1)
c           compute integral, -1 to 0, of p(x) and x*p(x).
c           -----------------------
            pint = pc(1)
            xpin = pc(1)/2.0d0
            tsign = 1.0d0
            do 30 i = 2, nq
               tsign = -tsign
               pint = pint + tsign*pc(i)/i
               xpin = xpin + tsign*pc(i)/(i+1)
   30       continue
c           store coefficients in elco and tesco.
c           --------------------------------
            elco(1,nq) = pint*rq1fac
            elco(2,nq) = 1.0d0
            do 40 i = 2, nq
               elco(i+1,nq) = rq1fac*pc(i)/i
   40       continue
            agamq = rqfac*xpin
            ragq = 1.0d0/agamq
            tesco(2,nq) = ragq
            if (nq .lt. 12) tesco(1,nqp1) = ragq*rqfac/nqp1
            tesco(3,nqm1) = ragq
   50    continue
      go to 100
c
   60 continue
         pc(1) = 1.0d0
         rq1fac = 1.0d0
         do 90 nq = 1, 5
c           ------------------------------------------------------------
c            the pc array will contain the coefficients of the
c                polynomial p(x) = (x+1)*(x+2)*...*(x+nq).
c            initially, p(x) = 1.
c           ------------------------------------------------------------
            fnq = nq
            nqp1 = nq + 1
c           form coefficients of p(x)*(x+nq).
c           ------------------------------------
            pc(nqp1) = 0.0d0
            do 70 ib = 1, nq
               i = nq + 2 - ib
               pc(i) = pc(i-1) + fnq*pc(i)
   70       continue
            pc(1) = fnq*pc(1)
c           store coefficients in elco and tesco.
c           --------------------------------
            do 80 i = 1, nqp1
               elco(i,nq) = pc(i)/pc(2)
   80       continue
            elco(2,nq) = 1.0d0
            tesco(1,nq) = rq1fac
            tesco(2,nq) = nqp1/elco(1,nq)
            tesco(3,nq) = (nq+2)/elco(1,nq)
            rq1fac = rq1fac/fnq
   90    continue
  100 continue
      return
c     ----------------------- end of subroutine dcfod
c     -----------------------
      end
