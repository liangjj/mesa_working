*deck slvs
      subroutine slvs (wm, iwm, x, tem)
c***begin prologue  slvs
c***subsidiary
c***purpose  subsidiary to debdf
c***library   slatec
c***type      single precision (slvs-s, dslvs-d)
c***author  watts, h. a., (snla)
c***description
c
c   slvs solves the linear system in the iteration scheme for the
c   integrator package debdf.
c
c***see also  debdf
c***routines called  sgbsl, sgesl
c***common blocks    debdf1
c***revision history  (yymmdd)
c   800901  date written
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910722  updated author section.  (als)
c   920422  changed dimension statement.  (wrb)
c***end prologue  slvs
c
clll. optimize
      integer iwm, i, ier, iownd, iowns, jstart, kflag, l, maxord,
     1   meband, meth, miter, ml, mu, n, nfe, nje, nq, nqu, nst
      real wm, x, tem,
     1   rownd, rowns, el0, h, hmin, hmxi, hu, tn, uround,
     2   di, hl0, phl0, r
      dimension wm(*), iwm(*), x(*), tem(*)
      common /debdf1/ rownd, rowns(210),
     1   el0, h, hmin, hmxi, hu, tn, uround, iownd(14), iowns(6),
     2   ier, jstart, kflag, l, meth, miter, maxord, n, nq, nst, nfe,
     3   nje, nqu
c-----------------------------------------------------------------------
c this routine manages the solution of the linear system arising from
c a chord iteration.  it is called by stod  if miter .ne. 0.
c if miter is 1 or 2, it calls sgesl to accomplish this.
c if miter = 3 it updates the coefficient h*el0 in the diagonal
c matrix, and then computes the solution.
c if miter is 4 or 5, it calls sgbsl.
c communication with slvs uses the following variables..
c wm  = real work space containing the inverse diagonal matrix if miter
c       is 3 and the lu decomposition of the matrix otherwise.
c       storage of matrix elements starts at wm(3).
c       wm also contains the following matrix-related data..
c       wm(1) = sqrt(uround) (not used here),
c       wm(2) = hl0, the previous value of h*el0, used if miter = 3.
c iwm = integer work space containing pivot information, starting at
c       iwm(21), if miter is 1, 2, 4, or 5.  iwm also contains the
c       band parameters ml = iwm(1) and mu = iwm(2) if miter is 4 or 5.
c x   = the right-hand side vector on input, and the solution vector
c       on output, of length n.
c tem = vector of work space of length n, not used in this version.
c ier = output flag (in common).  ier = 0 if no trouble occurred.
c       ier = -1 if a singular matrix arose with miter = 3.
c this routine also uses the common variables el0, h, miter, and n.
c-----------------------------------------------------------------------
c***first executable statement  slvs
      ier = 0
      go to (100, 100, 300, 400, 400), miter
 100  call sgesl (wm(3), n, n, iwm(21), x, 0)
      return
c
 300  phl0 = wm(2)
      hl0 = h*el0
      wm(2) = hl0
      if (hl0 .eq. phl0) go to 330
      r = hl0/phl0
      do 320 i = 1,n
        di = 1.0e0 - r*(1.0e0 - 1.0e0/wm(i+2))
        if (abs(di) .eq. 0.0e0) go to 390
 320    wm(i+2) = 1.0e0/di
 330  do 340 i = 1,n
 340    x(i) = wm(i+2)*x(i)
      return
 390  ier = -1
      return
c
 400  ml = iwm(1)
      mu = iwm(2)
      meband = 2*ml + mu + 1
      call sgbsl (wm(3), meband, n, ml, mu, iwm(21), x, 0)
      return
c----------------------- end of subroutine slvs -----------------------
      end
