*deck sdpsc
      subroutine sdpsc (ksgn, n, nq, yh)
c***begin prologue  sdpsc
c***subsidiary
c***purpose  subroutine sdpsc computes the predicted yh values by
c            effectively multiplying the yh array by the pascal triangle
c            matrix when ksgn is +1, and performs the inverse function
c            when ksgn is -1.
c***library   slatec (sdrive)
c***type      single precision (sdpsc-s, ddpsc-d, cdpsc-c)
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
c***end prologue  sdpsc
      integer i, j, j1, j2, ksgn, n, nq
      real yh(n,*)
c***first executable statement  sdpsc
      if (ksgn .gt. 0) then
        do 10 j1 = 1,nq
          do 10 j2 = j1,nq
            j = nq - j2 + j1
            do 10 i = 1,n
 10           yh(i,j) = yh(i,j) + yh(i,j+1)
      else
        do 30 j1 = 1,nq
          do 30 j2 = j1,nq
            j = nq - j2 + j1
            do 30 i = 1,n
 30           yh(i,j) = yh(i,j) - yh(i,j+1)
      end if
      return
      end
