*deck efc
      subroutine efc (ndata, xdata, ydata, sddata, nord, nbkpt, bkpt,
     +   mdein, mdeout, coeff, lw, w)
c***begin prologue  efc
c***purpose  fit a piecewise polynomial curve to discrete data.
c            the piecewise polynomials are represented as b-splines.
c            the fitting is done in a weighted least squares sense.
c***library   slatec
c***category  k1a1a1, k1a2a, l8a3
c***type      single precision (efc-s, defc-d)
c***keywords  b-spline, curve fitting, weighted least squares
c***author  hanson, r. j., (snla)
c***description
c
c      this subprogram fits a piecewise polynomial curve
c      to discrete data.  the piecewise polynomials are
c      represented as b-splines.
c      the fitting is done in a weighted least squares sense.
c
c      the data can be processed in groups of modest size.
c      the size of the group is chosen by the user.  this feature
c      may be necessary for purposes of using constrained curve fitting
c      with subprogram fc( ) on a very large data set.
c
c      for a description of the b-splines and usage instructions to
c      evaluate them, see
c
c      c. w. de boor, package for calculating with b-splines.
c                     siam j. numer. anal., p. 441, (june, 1977).
c
c      for further discussion of (constrained) curve fitting using
c      b-splines, see
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
c                         internal to  efc( ) the extreme end knots may
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
c        mdein
c                         an integer flag, with one of two possible
c                         values (1 or 2), that directs the subprogram
c                         action with regard to new data points provided
c                         by the user.
c
c                         =1  the first time that efc( ) has been
c                         entered.  there are ndata points to process.
c
c                         =2  this is another entry to efc( ).  the sub-
c                         program efc( ) has been entered with mdein=1
c                         exactly once before for this problem.  there
c                         are ndata new additional points to merge and
c                         process with any previous points.
c                         (when using efc( ) with mdein=2 it is import-
c                         ant that the set of knots remain fixed at the
c                         same values for all entries to efc( ).)
c       lw
c                         the amount of working storage actually
c                         allocated for the working array w(*).
c                         this quantity is compared with the
c                         actual amount of storage needed in efc( ).
c                         insufficient storage allocated for w(*) is
c                         an error.  this feature was included in efc( )
c                         because misreading the storage formula
c                         for w(*) might very well lead to subtle
c                         and hard-to-find programming bugs.
c
c                         the length of the array w(*) must satisfy
c
c                         lw .ge. (nbkpt-nord+3)*(nord+1)+
c                                 (nbkpt+1)*(nord+1)+
c                               2*max(ndata,nbkpt)+nbkpt+nord**2
c
c  output..
c      mdeout
c                         an output flag that indicates the status
c                         of the curve fit.
c
c                         =-1  a usage error of efc( ) occurred.  the
c                         offending condition is noted with the slatec
c                         library error processor, xermsg( ).  in case
c                         the working array w(*) is not long enough, the
c                         minimal acceptable length is printed.
c
c                         =1  the b-spline coefficients for the fitted
c                         curve have been returned in array coeff(*).
c
c                         =2  not enough data has been processed to
c                         determine the b-spline coefficients.
c                         the user has one of two options.  continue
c                         to process more data until a unique set
c                         of coefficients is obtained, or use the
c                         subprogram fc( ) to obtain a specific
c                         set of coefficients.  the user should read
c                         the usage instructions for fc( ) for further
c                         details if this second option is chosen.
c      coeff(*)
c                         if the output value of mdeout=1, this array
c                         contains the unknowns obtained from the least
c                         squares fitting process.  these n=nbkpt-nord
c                         parameters are the b-spline coefficients.
c                         for mdeout=2, not enough data was processed to
c                         uniquely determine the b-spline coefficients.
c                         in this case, and also when mdeout=-1, all
c                         values of coeff(*) are set to zero.
c
c                         if the user is not satisfied with the fitted
c                         curve returned by efc( ), the constrained
c                         least squares curve fitting subprogram fc( )
c                         may be required.  the work done within efc( )
c                         to accumulate the data can be utilized by
c                         the user, if so desired.  this involves
c                         saving the first (nbkpt-nord+3)*(nord+1)
c                         entries of w(*) and providing this data
c                         to fc( ) with the "old problem" designation.
c                         the user should read the usage instructions
c                         for subprogram fc( ) for further details.
c
c  working array..
c      w(*)
c                         this array is typed real.
c                         its length is  specified as an input parameter
c                         in lw as noted above.  the contents of w(*)
c                         must not be modified by the user between calls
c                         to efc( ) with values of mdein=1,2,2,... .
c                         the first (nbkpt-nord+3)*(nord+1) entries of
c                         w(*) are acceptable as direct input to fc( )
c                         for an "old problem" only when mdeout=1 or 2.
c
c  evaluating the
c  fitted curve..
c                         to evaluate derivative number ider at xval,
c                         use the function subprogram bvalu( ).
c
c                         f = bvalu(bkpt,coeff,nbkpt-nord,nord,ider,
c                                      xval,inbv,workb)
c
c                         the output of this subprogram will not be
c                         defined unless an output value of mdeout=1
c                         was obtained from efc( ), xval is in the data
c                         interval, and ider is nonnegative and .lt.
c                         nord.
c
c                         the first time bvalu( ) is called, inbv=1
c                         must be specified.  this value of inbv is the
c                         overwritten by bvalu( ).  the array workb(*)
c                         must be of length at least 3*nord, and must
c                         not be the same as the w(*) array used in the
c                         call to efc( ).
c
c                         bvalu( ) expects the breakpoint array bkpt(*)
c                         to be sorted.
c
c***references  r. j. hanson, constrained least squares curve fitting
c                 to discrete data using b-splines, a users guide,
c                 report sand78-1291, sandia laboratories, december
c                 1978.
c***routines called  efcmn
c***revision history  (yymmdd)
c   800801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   900510  change prologue comments to refer to xermsg.  (rwc)
c   900607  editorial changes to prologue to make prologues for efc,
c           defc, fc, and dfc look as much the same as possible.  (rwc)
c   920501  reformatted the references section.  (wrb)
c***end prologue  efc
c
c      subroutine           function/remarks
c
c      bsplvn( )          compute function values of b-splines.  from
c                         the b-spline package of de boor noted above.
c
c      bndacc( ),         banded least squares matrix processors.
c      bndsol( )          from lawson-hanson, solving least
c                         squares problems.
c
c      ssort( )           data sorting subroutine, from the
c                         sandia math. library, sand77-1441.
c
c      xermsg( )          error handling routine
c                         for the slatec math. library.
c                         see sand78-1189, by r. e. jones.
c
c      scopy( ),sscal( )  subroutines from the blas package.
c
c                         written by r. hanson, sandia natl. labs.,
c                         alb., n. m., august-september, 1980.
c
      real bkpt(*),coeff(*),sddata(*),w(*),xdata(*),ydata(*)
      integer lw, mdein, mdeout, nbkpt, ndata, nord
c
      external efcmn
c
      integer lbf, lbkpt, lg, lptemp, lww, lxtemp, mdg, mdw
c
c***first executable statement  efc
c     lww=1               usage in efcmn( ) of w(*)..
c     lww,...,lg-1        w(*,*)
c
c     lg,...,lxtemp-1     g(*,*)
c
c     lxtemp,...,lptemp-1 xtemp(*)
c
c     lptemp,...,lbkpt-1  ptemp(*)
c
c     lbkpt,...,lbf       bkpt(*) (local to efcmn( ))
c
c     lbf,...,lbf+nord**2 bf(*,*)
c
      mdg = nbkpt+1
      mdw = nbkpt-nord+3
      lww = 1
      lg = lww + mdw*(nord+1)
      lxtemp = lg + mdg*(nord+1)
      lptemp = lxtemp + max(ndata,nbkpt)
      lbkpt  = lptemp + max(ndata,nbkpt)
      lbf = lbkpt + nbkpt
      call efcmn(ndata,xdata,ydata,sddata,
     1         nord,nbkpt,bkpt,
     2         mdein,mdeout,
     3         coeff,
     4         w(lbf),w(lxtemp),w(lptemp),w(lbkpt),
     5         w(lg),mdg,w(lww),mdw,lw)
      return
      end
