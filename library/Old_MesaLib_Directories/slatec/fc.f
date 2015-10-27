*deck fc
      subroutine fc (ndata, xdata, ydata, sddata, nord, nbkpt, bkpt,
     +   nconst, xconst, yconst, nderiv, mode, coeff, w, iw)
c***begin prologue  fc
c***purpose  fit a piecewise polynomial curve to discrete data.
c            the piecewise polynomials are represented as b-splines.
c            the fitting is done in a weighted least squares sense.
c            equality and inequality constraints can be imposed on the
c            fitted curve.
c***library   slatec
c***category  k1a1a1, k1a2a, l8a3
c***type      single precision (fc-s, dfc-d)
c***keywords  b-spline, constrained least squares, curve fitting,
c             weighted least squares
c***author  hanson, r. j., (snla)
c***description
c
c      this subprogram fits a piecewise polynomial curve
c      to discrete data.  the piecewise polynomials are
c      represented as b-splines.
c      the fitting is done in a weighted least squares sense.
c      equality and inequality constraints can be imposed on the
c      fitted curve.
c
c      for a description of the b-splines and usage instructions to
c      evaluate them, see
c
c      c. w. de boor, package for calculating with b-splines.
c                     siam j. numer. anal., p. 441, (june, 1977).
c
c      for further documentation and discussion of constrained
c      curve fitting using b-splines, see
c
c      r. j. hanson, constrained least squares curve fitting
c                   to discrete data using b-splines, a user's
c                   guide. sandia labs. tech. rept. sand-78-1291,
c                   december, (1978).
c
c  input..
c      ndata,xdata(*),
c      ydata(*),
c      sddata(*)
c                         the ndata discrete (x,y) pairs and the y value
c                         standard deviation or uncertainty, sd, are in
c                         the respective arrays xdata(*), ydata(*), and
c                         sddata(*).  no sorting of xdata(*) is
c                         required.  any non-negative value of ndata is
c                         allowed.  a negative value of ndata is an
c                         error.  a zero value for any entry of
c                         sddata(*) will weight that data point as 1.
c                         otherwise the weight of that data point is
c                         the reciprocal of this entry.
c
c      nord,nbkpt,
c      bkpt(*)
c                         the nbkpt knots of the b-spline of order nord
c                         are in the array bkpt(*).  normally the
c                         problem data interval will be included between
c                         the limits bkpt(nord) and bkpt(nbkpt-nord+1).
c                         the additional end knots bkpt(i),i=1,...,
c                         nord-1 and i=nbkpt-nord+2,...,nbkpt, are
c                         required to compute the functions used to fit
c                         the data.  no sorting of bkpt(*) is required.
c                         internal to  fc( ) the extreme end knots may
c                         be reduced and increased respectively to
c                         accommodate any data values that are exterior
c                         to the given knot values.  the contents of
c                         bkpt(*) is not changed.
c
c                         nord must be in the range 1 .le. nord .le. 20.
c                         the value of nbkpt must satisfy the condition
c                         nbkpt .ge. 2*nord.
c                         other values are considered errors.
c
c                         (the order of the spline is one more than the
c                         degree of the piecewise polynomial defined on
c                         each interval.  this is consistent with the
c                         b-spline package convention.  for example,
c                         nord=4 when we are using piecewise cubics.)
c
c      nconst,xconst(*),
c      yconst(*),nderiv(*)
c                         the number of conditions that constrain the
c                         b-spline is nconst.  a constraint is specified
c                         by an (x,y) pair in the arrays xconst(*) and
c                         yconst(*), and by the type of constraint and
c                         derivative value encoded in the array
c                         nderiv(*).  no sorting of xconst(*) is
c                         required.  the value of nderiv(*) is
c                         determined as follows.  suppose the i-th
c                         constraint applies to the j-th derivative
c                         of the b-spline.  (any non-negative value of
c                         j < nord is permitted.  in particular the
c                         value j=0 refers to the b-spline itself.)
c                         for this i-th constraint, set
c                          xconst(i)=x,
c                          yconst(i)=y, and
c                          nderiv(i)=itype+4*j, where
c
c                          itype = 0,      if (j-th deriv. at x) .le. y.
c                                = 1,      if (j-th deriv. at x) .ge. y.
c                                = 2,      if (j-th deriv. at x) .eq. y.
c                                = 3,      if (j-th deriv. at x) .eq.
c                                             (j-th deriv. at y).
c                          (a value of nderiv(i)=-1 will cause this
c                          constraint to be ignored.  this subprogram
c                          feature is often useful when temporarily
c                          suppressing a constraint while still
c                          retaining the source code of the calling
c                          program.)
c
c        mode
c                         an input flag that directs the least squares
c                         solution method used by fc( ).
c
c                         the variance function, referred to below,
c                         defines the square of the probable error of
c                         the fitted curve at any point, xval.
c                         this feature of  fc( ) allows one to use the
c                         square root of this variance function to
c                         determine a probable error band around the
c                         fitted curve.
c
c                         =1  a new problem.  no variance function.
c
c                         =2  a new problem.  want variance function.
c
c                         =3  an old problem.  no variance function.
c
c                         =4  an old problem.  want variance function.
c
c                         any value of mode other than 1-4 is an error.
c
c                         the user with a new problem can skip directly
c                         to the description of the input parameters
c                         iw(1), iw(2).
c
c                         if the user correctly specifies the new or old
c                         problem status, the subprogram fc( ) will
c                         perform more efficiently.
c                         by an old problem it is meant that subprogram
c                         fc( ) was last called with this same set of
c                         knots, data points and weights.
c
c                         another often useful deployment of this old
c                         problem designation can occur when one has
c                         previously obtained a q-r orthogonal
c                         decomposition of the matrix resulting from
c                         b-spline fitting of data (without constraints)
c                         at the breakpoints bkpt(i), i=1,...,nbkpt.
c                         for example, this matrix could be the result
c                         of sequential accumulation of the least
c                         squares equations for a very large data set.
c                         the user writes this code in a manner
c                         convenient for the application.  for the
c                         discussion here let
c
c                                      n=nbkpt-nord, and k=n+3
c
c                         let us assume that an equivalent least squares
c                         system
c
c                                      rc=d
c
c                         has been obtained.  here r is an n+1 by n
c                         matrix and d is a vector with n+1 components.
c                         the last row of r is zero.  the matrix r is
c                         upper triangular and banded.  at most nord of
c                         the diagonals are nonzero.
c                         the contents of r and d can be copied to the
c                         working array w(*) as follows.
c
c                         the i-th diagonal of r, which has n-i+1
c                         elements, is copied to w(*) starting at
c
c                                      w((i-1)*k+1),
c
c                         for i=1,...,nord.
c                         the vector d is copied to w(*) starting at
c
c                                      w(nord*k+1)
c
c                         the input value used for ndata is arbitrary
c                         when an old problem is designated.  because
c                         of the feature of fc( ) that checks the
c                         working storage array lengths, a value not
c                         exceeding nbkpt should be used.  for example,
c                         use ndata=0.
c
c                         (the constraints or variance function request
c                         can change in each call to fc( ).)  a new
c                         problem is anything other than an old problem.
c
c      iw(1),iw(2)
c                         the amounts of working storage actually
c                         allocated for the working arrays w(*) and
c                         iw(*).  these quantities are compared with the
c                         actual amounts of storage needed in fc( ).
c                         insufficient storage allocated for either
c                         w(*) or iw(*) is an error.  this feature was
c                         included in fc( ) because misreading the
c                         storage formulas for w(*) and iw(*) might very
c                         well lead to subtle and hard-to-find
c                         programming bugs.
c
c                         the length of w(*) must be at least
c
c                           nb=(nbkpt-nord+3)*(nord+1)+
c                               2*max(ndata,nbkpt)+nbkpt+nord**2
c
c                         whenever possible the code uses banded matrix
c                         processors bndacc( ) and bndsol( ).  these
c                         are utilized if there are no constraints,
c                         no variance function is required, and there
c                         is sufficient data to uniquely determine the
c                         b-spline coefficients.  if the band processors
c                         cannot be used to determine the solution,
c                         then the constrained least squares code lsei
c                         is used.  in this case the subprogram requires
c                         an additional block of storage in w(*).  for
c                         the discussion here define the integers neqcon
c                         and nincon respectively as the number of
c                         equality (itype=2,3) and inequality
c                         (itype=0,1) constraints imposed on the fitted
c                         curve.  define
c
c                           l=nbkpt-nord+1
c
c                         and note that
c
c                           nconst=neqcon+nincon.
c
c                         when the subprogram fc( ) uses lsei( ) the
c                         length of the working array w(*) must be at
c                         least
c
c                           lw=nb+(l+nconst)*l+
c                              2*(neqcon+l)+(nincon+l)+(nincon+2)*(l+6)
c
c                         the length of the array iw(*) must be at least
c
c                           iw1=nincon+2*l
c
c                         in any case.
c
c  output..
c      mode
c                         an output flag that indicates the status
c                         of the constrained curve fit.
c
c                         =-1  a usage error of fc( ) occurred.  the
c                              offending condition is noted with the
c                              slatec library error processor, xermsg.
c                              in case the working arrays w(*) or iw(*)
c                              are not long enough, the minimal
c                              acceptable length is printed.
c
c                         = 0  successful constrained curve fit.
c
c                         = 1  the requested equality constraints
c                              are contradictory.
c
c                         = 2  the requested inequality constraints
c                              are contradictory.
c
c                         = 3  both equality and inequality constraints
c                              are contradictory.
c
c      coeff(*)
c                         if the output value of mode=0 or 1, this array
c                         contains the unknowns obtained from the least
c                         squares fitting process.  these n=nbkpt-nord
c                         parameters are the b-spline coefficients.
c                         for mode=1, the equality constraints are
c                         contradictory.  to make the fitting process
c                         more robust, the equality constraints are
c                         satisfied in a least squares sense.  in this
c                         case the array coeff(*) contains b-spline
c                         coefficients for this extended concept of a
c                         solution.  if mode=-1,2 or 3 on output, the
c                         array coeff(*) is undefined.
c
c  working arrays..
c      w(*),iw(*)
c                         these arrays are respectively typed real and
c                         integer.
c                         their required lengths are specified as input
c                         parameters in iw(1), iw(2) noted above.  the
c                         contents of w(*) must not be modified by the
c                         user if the variance function is desired.
c
c  evaluating the
c  variance function..
c                         to evaluate the variance function (assuming
c                         that the uncertainties of the y values were
c                         provided to  fc( ) and an input value of
c                         mode=2 or 4 was used), use the function
c                         subprogram  cv( )
c
c                           var=cv(xval,ndata,nconst,nord,nbkpt,
c                                  bkpt,w)
c
c                         here xval is the point where the variance is
c                         desired.  the other arguments have the same
c                         meaning as in the usage of fc( ).
c
c                         for those users employing the old problem
c                         designation, let mdata be the number of data
c                         points in the problem.  (this may be different
c                         from ndata if the old problem designation
c                         feature was used.)  the value, var, should be
c                         multiplied by the quantity
c
c                         real(max(ndata-n,1))/max(mdata-n,1)
c
c                         the output of this subprogram is not defined
c                         if an input value of mode=1 or 3 was used in
c                         fc( ) or if an output value of mode=-1, 2, or
c                         3 was obtained.  the variance function, except
c                         for the scaling factor noted above, is given
c                         by
c
c                           var=(transpose of b(xval))*c*b(xval)
c
c                         the vector b(xval) is the b-spline basis
c                         function values at x=xval.
c                         the covariance matrix, c, of the solution
c                         coefficients accounts only for the least
c                         squares equations and the explicitly stated
c                         equality constraints.  this fact must be
c                         considered when interpreting the variance
c                         function from a data fitting problem that has
c                         inequality constraints on the fitted curve.
c
c  evaluating the
c  fitted curve..
c                         to evaluate derivative number ider at xval,
c                         use the function subprogram bvalu( ).
c
c                           f = bvalu(bkpt,coeff,nbkpt-nord,nord,ider,
c                                      xval,inbv,workb)
c
c                         the output of this subprogram will not be
c                         defined unless an output value of mode=0 or 1
c                         was obtained from  fc( ), xval is in the data
c                         interval, and ider is nonnegative and .lt.
c                         nord.
c
c                         the first time bvalu( ) is called, inbv=1
c                         must be specified.  this value of inbv is the
c                         overwritten by bvalu( ).  the array workb(*)
c                         must be of length at least 3*nord, and must
c                         not be the same as the w(*) array used in
c                         the call to fc( ).
c
c                         bvalu( ) expects the breakpoint array bkpt(*)
c                         to be sorted.
c
c***references  r. j. hanson, constrained least squares curve fitting
c                 to discrete data using b-splines, a users guide,
c                 report sand78-1291, sandia laboratories, december
c                 1978.
c***routines called  fcmn
c***revision history  (yymmdd)
c   780801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  convert references to xerrwv to references to xermsg.  (rwc)
c   900607  editorial changes to prologue to make prologues for efc,
c           defc, fc, and dfc look as much the same as possible.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  fc
      real             bkpt(*), coeff(*), sddata(*), w(*), xconst(*),
     *   xdata(*), yconst(*), ydata(*)
      integer iw(*), mode, nbkpt, nconst, ndata, nderiv(*), nord
c
      external fcmn
c
      integer i1, i2, i3, i4, i5, i6, i7, mdg, mdw
c
c***first executable statement  fc
      mdg = nbkpt - nord + 3
      mdw = nbkpt - nord + 1 + nconst
c                         usage in fcmn( ) of w(*)..
c     i1,...,i2-1      g(*,*)
c
c     i2,...,i3-1      xtemp(*)
c
c     i3,...,i4-1      ptemp(*)
c
c     i4,...,i5-1      bkpt(*) (local to fcmn( ))
c
c     i5,...,i6-1      bf(*,*)
c
c     i6,...,i7-1      w(*,*)
c
c     i7,...           work(*) for lsei( )
c
      i1 = 1
      i2 = i1 + mdg*(nord+1)
      i3 = i2 + max(ndata,nbkpt)
      i4 = i3 + max(ndata,nbkpt)
      i5 = i4 + nbkpt
      i6 = i5 + nord*nord
      i7 = i6 + mdw*(nbkpt-nord+1)
      call  fcmn(ndata, xdata, ydata, sddata, nord, nbkpt, bkpt, nconst,
     1   xconst, yconst, nderiv, mode, coeff, w(i5), w(i2), w(i3),
     2   w(i4), w(i1), mdg, w(i6), mdw, w(i7), iw)
      return
      end
