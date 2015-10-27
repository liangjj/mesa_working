*deck dslvs
      subroutine dslvs (wm, iwm, x, tem)
c***begin prologue  dslvs
c***subsidiary
c***purpose  subsidiary to ddebdf
c***library   slatec
c***type      double precision (slvs-s, dslvs-d)
c***author  watts, h. a., (snla)
c***description
c
c   dslvs solves the linear system in the iteration scheme for the
c   integrator package ddebdf.
c
c***see also  ddebdf
c***routines called  dgbsl, dgesl
c***common blocks    ddebd1
c***revision history  (yymmdd)
c   820301  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c   920422  changed dimension statement.  (wrb)
c***end prologue  dslvs
c
      integer i, ier, iownd, iowns, iwm, jstart, kflag, l, maxord,
     1      meband, meth, miter, ml, mu, n, nfe, nje, nq, nqu, nst
      double precision di, el0, h, hl0, hmin, hmxi, hu, phl0,
     1      r, rownd, rowns, tem, tn, uround, wm, x
      dimension wm(*), iwm(*), x(*), tem(*)
      common /ddebd1/ rownd,rowns(210),el0,h,hmin,hmxi,hu,tn,uround,
     1                iownd(14),iowns(6),ier,jstart,kflag,l,meth,miter,
     2                maxord,n,nq,nst,nfe,nje,nqu
c     ------------------------------------------------------------------
c      this routine manages the solution of the linear system arising
c      from a chord iteration.  it is called by dstod  if miter .ne. 0.
c      if miter is 1 or 2, it calls dgesl to accomplish this.
c      if miter = 3 it updates the coefficient h*el0 in the diagonal
c      matrix, and then computes the solution.
c      if miter is 4 or 5, it calls dgbsl.
c      communication with dslvs uses the following variables..
c      wm  = double precision work space containing the inverse diagonal
c      matrix if miter
c            is 3 and the lu decomposition of the matrix otherwise.
c            storage of matrix elements starts at wm(3).
c            wm also contains the following matrix-related data..
c            wm(1) = sqrt(uround) (not used here),
c            wm(2) = hl0, the previous value of h*el0, used if miter =
c            3.
c      iwm = integer work space containing pivot information, starting
c            at iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains
c            the band parameters ml = iwm(1) and mu = iwm(2) if miter is
c            4 or 5.
c      x   = the right-hand side vector on input, and the solution
c            vector on output, of length n.
c      tem = vector of work space of length n, not used in this version.
c      ier = output flag (in common).  ier = 0 if no trouble occurred.
c            ier = -1 if a singular matrix arose with miter = 3.
c      this routine also uses the common variables el0, h, miter, and n.
c-----------------------------------------------------------------------
c     begin block permitting ...exits to 80
c        begin block permitting ...exits to 60
c***first executable statement  dslvs
            ier = 0
            go to (10,10,20,70,70), miter
   10       continue
            call dgesl(wm(3),n,n,iwm(21),x,0)
c     ......exit
            go to 80
c
   20       continue
            phl0 = wm(2)
            hl0 = h*el0
            wm(2) = hl0
            if (hl0 .eq. phl0) go to 40
               r = hl0/phl0
               do 30 i = 1, n
                  di = 1.0d0 - r*(1.0d0 - 1.0d0/wm(i+2))
c        .........exit
                  if (abs(di) .eq. 0.0d0) go to 60
                  wm(i+2) = 1.0d0/di
   30          continue
   40       continue
            do 50 i = 1, n
               x(i) = wm(i+2)*x(i)
   50       continue
c     ......exit
            go to 80
   60    continue
         ier = -1
c     ...exit
         go to 80
c
   70    continue
         ml = iwm(1)
         mu = iwm(2)
         meband = 2*ml + mu + 1
         call dgbsl(wm(3),meband,n,ml,mu,iwm(21),x,0)
   80 continue
      return
c     ----------------------- end of subroutine dslvs
c     -----------------------
      end
