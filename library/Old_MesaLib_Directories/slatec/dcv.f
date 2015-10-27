*deck dcv
      double precision function dcv (xval, ndata, nconst, nord, nbkpt,
     +   bkpt, w)
c***begin prologue  dcv
c***purpose  evaluate the variance function of the curve obtained
c            by the constrained b-spline fitting subprogram dfc.
c***library   slatec
c***category  l7a3
c***type      double precision (cv-s, dcv-d)
c***keywords  analysis of covariance, b-spline,
c             constrained least squares, curve fitting
c***author  hanson, r. j., (snla)
c***description
c
c     dcv( ) is a companion function subprogram for dfc( ).  the
c     documentation for dfc( ) has complete usage instructions.
c
c     dcv( ) is used to evaluate the variance function of the curve
c     obtained by the constrained b-spline fitting subprogram, dfc( ).
c     the variance function defines the square of the probable error
c     of the fitted curve at any point, xval.  one can use the square
c     root of this variance function to determine a probable error band
c     around the fitted curve.
c
c     dcv( ) is used after a call to dfc( ).  mode, an input variable to
c     dfc( ), is used to indicate if the variance function is desired.
c     in order to use dcv( ), mode must equal 2 or 4 on input to dfc( ).
c     mode is also used as an output flag from dfc( ).  check to make
c     sure that mode = 0 after calling dfc( ), indicating a successful
c     constrained curve fit.  the array sddata, as input to dfc( ), must
c     also be defined with the standard deviation or uncertainty of the
c     y values to use dcv( ).
c
c     to evaluate the variance function after calling dfc( ) as stated
c     above, use dcv( ) as shown here
c
c          var=dcv(xval,ndata,nconst,nord,nbkpt,bkpt,w)
c
c     the variance function is given by
c
c      var=(transpose of b(xval))*c*b(xval)/dble(max(ndata-n,1))
c
c     where n = nbkpt - nord.
c
c     the vector b(xval) is the b-spline basis function values at
c     x=xval.  the covariance matrix, c, of the solution coefficients
c     accounts only for the least squares equations and the explicitly
c     stated equality constraints.  this fact must be considered when
c     interpreting the variance function from a data fitting problem
c     that has inequality constraints on the fitted curve.
c
c     all the variables in the calling sequence for dcv( ) are used in
c     dfc( ) except the variable xval.  do not change the values of
c     these variables between the call to dfc( ) and the use of dcv( ).
c
c     the following is a brief description of the variables
c
c     xval    the point where the variance is desired, a double
c             precision variable.
c
c     ndata   the number of discrete (x,y) pairs for which dfc( )
c             calculated a piece-wise polynomial curve.
c
c     nconst  the number of conditions that constrained the b-spline in
c             dfc( ).
c
c     nord    the order of the b-spline used in dfc( ).
c             the value of nord must satisfy 1 < nord < 20 .
c
c             (the order of the spline is one more than the degree of
c             the piece-wise polynomial defined on each interval.  this
c             is consistent with the b-spline package convention.  for
c             example, nord=4 when we are using piece-wise cubics.)
c
c     nbkpt   the number of knots in the array bkpt(*).
c             the value of nbkpt must satisfy nbkpt .ge. 2*nord.
c
c     bkpt(*) the double precision array of knots.  normally the problem
c             data interval will be included between the limits
c             bkpt(nord) and bkpt(nbkpt-nord+1).  the additional end
c             knots bkpt(i),i=1,...,nord-1 and i=nbkpt-nord+2,...,nbkpt,
c             are required by dfc( ) to compute the functions used to
c             fit the data.
c
c     w(*)    double precision work array as used in dfc( ).  see dfc( )
c             for the required length of w(*).  the contents of w(*)
c             must not be modified by the user if the variance function
c             is desired.
c
c***references  r. j. hanson, constrained least squares curve fitting
c                 to discrete data using b-splines, a users guide,
c                 report sand78-1291, sandia laboratories, december
c                 1978.
c***routines called  ddot, dfspvn
c***revision history  (yymmdd)
c   780801  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890831  modified array declarations.  (wrb)
c   890911  removed unnecessary intrinsics.  (wrb)
c   891006  cosmetic changes to prologue.  (wrb)
c   891006  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  dcv
      integer i, ileft, ip, is, last, mdg, mdw, n, nbkpt, nconst,
     *     ndata, nord
      double precision bkpt, ddot, v, w, xval, zero
      dimension bkpt(*),w(*),v(40)
c***first executable statement  dcv
      zero = 0.0d0
      mdg = nbkpt - nord + 3
      mdw = nbkpt - nord + 1 + nconst
      is = mdg*(nord + 1) + 2*max(ndata,nbkpt) + nbkpt + nord**2
      last = nbkpt - nord + 1
      ileft = nord
   10 if (xval .lt. bkpt(ileft+1) .or. ileft .ge. last - 1) go to 20
         ileft = ileft + 1
      go to 10
   20 continue
      call dfspvn(bkpt,nord,1,xval,ileft,v(nord+1))
      ileft = ileft - nord + 1
      ip = mdw*(ileft - 1) + ileft + is
      n = nbkpt - nord
      do 30 i = 1, nord
         v(i) = ddot(nord,w(ip),1,v(nord+1),1)
         ip = ip + mdw
   30 continue
      dcv = max(ddot(nord,v,1,v(nord+1),1),zero)
c
c     scale the variance so it is an unbiased estimate.
      dcv = dcv/max(ndata-n,1)
      return
      end
