*deck dbvder
      subroutine dbvder (x, y, yp, g, ipar)
c***begin prologue  dbvder
c***subsidiary
c***purpose  subsidiary to dbvsup
c***library   slatec
c***type      double precision (bvder-s, dbvder-d)
c***author  watts, h. a., (snla)
c***description
c
c **********************************************************************
c     nfc = number of base solution vectors
c
c     ncomp = number of components per solution vector
c
c              1 -- nonzero particular solution
c     inhomo =
c              2 or 3 -- zero particular solution
c
c             0 -- inhomogeneous vector term g(x) identically zero
c     igofx =
c             1 -- inhomogeneous vector term g(x) not identically zero
c
c     g = inhomogeneous vector term g(x)
c
c     xsav = previous value of x
c
c     c = normalization factor for the particular solution
c
c           0   ( if  neqivp = 0 )
c     ivp =
c           number of differential equations integrated due to
c           the original boundary value problem   ( if  neqivp .gt. 0 )
c
c     nofst - for problems with auxiliary initial value equations,
c             nofst communicates to the routine dfmat how to access
c             the dependent variables corresponding to this initial
c             value problem.  for example, during any call to dfmat,
c             the first dependent variable for the initial value
c             problem is in position  y(nofst + 1).
c             see example in sand77-1328.
c **********************************************************************
c
c***see also  dbvsup
c***routines called  (none)
c***common blocks    dml8sz, dmlivp
c***revision history  (yymmdd)
c   750601  date written
c   890831  modified array declarations.  (wrb)
c   890921  realigned order of variables in certain common blocks.
c           (wrb)
c   891214  prologue converted to version 4.0 format.  (bab)
c   900328  added type section.  (wrb)
c   910701  corrected routines called section.  (wrb)
c   910722  updated author section.  (als)
c   920618  minor restructuring of code.  (rwc, wrb)
c***end prologue  dbvder
      integer igofx, inhomo, ipar, ivp, j, k, l, na, ncomp, nfc, nofst
      double precision c, g(*), x, xsav, y(*), yp(*)
c
c **********************************************************************
c
      common /dml8sz/ c,xsav,igofx,inhomo,ivp,ncomp,nfc
c
c **********************************************************************
c     the common block below is used to communicate with the user
c     supplied subroutine dfmat.  the user should not alter this
c     common block.
c
      common /dmlivp/ nofst
c **********************************************************************
c
c***first executable statement  dbvder
      if (ivp .gt. 0) call duivp(x,y(ivp+1),yp(ivp+1))
      nofst = ivp
      na = 1
      do 10 k=1,nfc
         call dfmat(x,y(na),yp(na))
         nofst = nofst - ncomp
         na = na + ncomp
   10 continue
c
      if (inhomo .ne. 1) return
      call dfmat(x,y(na),yp(na))
c
      if (igofx .eq. 0) return
      if (x .ne. xsav) then
         if (ivp .eq. 0) call dgvec(x,g)
         if (ivp .gt. 0) call duvec(x,y(ivp+1),g)
         xsav = x
      endif
c
c     if the user has chosen not to normalize the particular
c     solution, then c is defined in dbvpor to be 1.0
c
c     the following loop is just
c     call daxpy (ncomp, 1.0d0/c, g, 1, yp(na), 1)
c
      do 20 j=1,ncomp
         l = na + j - 1
         yp(l) = yp(l) + g(j)/c
   20 continue
      return
      end
