*deck cdntp
      subroutine cdntp (h, k, n, nq, t, tout, yh, y)
c***begin prologue  cdntp
c***subsidiary
c***purpose  subroutine cdntp interpolates the k-th derivative of y at
c            tout, using the data in the yh array.  if k has a value
c            greater than nq, the nq-th derivative is calculated.
c***library   slatec (sdrive)
c***type      complex (sdntp-s, ddntp-d, cdntp-c)
c***author  kahaner, d. k., (nist)
c             national institute of standards and technology
c             gaithersburg, md  20899
c           sutherland, c. d., (lanl)
c             mail stop d466
c             los alamos national laboratory
c             los alamos, nm  87545
c***routines called  (none)
c***revision history  (yymmdd)
c   790601  date written
c   900329  initial submission to slatec.
c***end prologue  cdntp
      integer i, j, jj, k, kk, kused, n, nq
      complex y(*), yh(n,*)
      real factor, h, r, t, tout
c***first executable statement  cdntp
      if (k .eq. 0) then
        do 10 i = 1,n
 10       y(i) = yh(i,nq+1)
        r = ((tout - t)/h)
        do 20 jj = 1,nq
          j = nq + 1 - jj
          do 20 i = 1,n
 20         y(i) = yh(i,j) + r*y(i)
      else
        kused = min(k, nq)
        factor = 1.e0
        do 40 kk = 1,kused
 40       factor = factor*(nq+1-kk)
        do 50 i = 1,n
 50       y(i) = factor*yh(i,nq+1)
        r = ((tout - t)/h)
        do 80 jj = kused+1,nq
          j = kused + 1 + nq - jj
          factor = 1.e0
          do 60 kk = 1,kused
 60         factor = factor*(j-kk)
          do 70 i = 1,n
 70         y(i) = factor*yh(i,j) + r*y(i)
 80       continue
        do 100 i = 1,n
 100      y(i) = y(i)*h**(-kused)
      end if
      return
      end
