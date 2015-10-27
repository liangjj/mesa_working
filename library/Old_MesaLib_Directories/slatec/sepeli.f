*deck sepeli
      subroutine sepeli (intl, iorder, a, b, m, mbdcnd, bda, alpha, bdb,
     +   beta, c, d, n, nbdcnd, bdc, gama, bdd, xnu, cofx, cofy, grhs,
     +   usol, idmn, w, pertrb, ierror)
c***begin prologue  sepeli
c***purpose  discretize and solve a second and, optionally, a fourth
c            order finite difference approximation on a uniform grid to
c            the general separable elliptic partial differential
c            equation on a rectangle with any combination of periodic or
c            mixed boundary conditions.
c***library   slatec (fishpack)
c***category  i2b1a2
c***type      single precision (sepeli-s)
c***keywords  elliptic, fishpack, helmholtz, pde, separable
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c dimension of           bda(n+1), bdb(n+1), bdc(m+1), bdd(m+1),
c arguments              usol(idmn,n+1), grhs(idmn,n+1),
c                        w (see argument list)
c
c latest revision        march 1977
c
c purpose                sepeli solves for either the second-order
c                        finite difference approximation or a
c                        fourth-order approximation to a separable
c                        elliptic equation.
c
c                                    2    2
c                             af(x)*d u/dx + bf(x)*du/dx  + cf(x)*u +
c                                    2    2
c                             df(y)*d u/dy  + ef(y)*du/dy + ff(y)*u
c
c                             = g(x,y)
c
c                        on a rectangle (x greater than or equal to a
c                        and less than or equal to b; y greater than
c                        or equal to c and less than or equal to d).
c                        any combination of periodic or mixed boundary
c                        conditions is allowed.
c
c purpose                the possible boundary conditions are:
c                        in the x-direction:
c                         (0) periodic, u(x+b-a,y)=u(x,y) for all y,x
c                         (1) u(a,y), u(b,y) are specified for all y
c                         (2) u(a,y), du(b,y)/dx+beta*u(b,y) are
c                             specified for all y
c                         (3) du(a,y)/dx+alpha*u(a,y),du(b,y)/dx+
c                             beta*u(b,y) are specified for all y
c                         (4) du(a,y)/dx+alpha*u(a,y),u(b,y) are
c                             specified for all y
c
c                        in the y-direction:
c                         (0) periodic, u(x,y+d-c)=u(x,y) for all x,y
c                         (1) u(x,c),u(x,d) are specified for all x
c                         (2) u(x,c),du(x,d)/dy+xnu*u(x,d) are specified
c                             for all x
c                         (3) du(x,c)/dy+gama*u(x,c),du(x,d)/dy+
c                             xnu*u(x,d) are specified for all x
c                         (4) du(x,c)/dy+gama*u(x,c),u(x,d) are
c                             specified for all x
c
c arguments
c
c on input               intl
c                          = 0 on initial entry to sepeli or if any of
c                              the arguments c, d, n, nbdcnd, cofy are
c                              changed from a previous call
c                          = 1 if c, d, n, nbdcnd, cofy are unchanged
c                              from the previous call.
c
c                        iorder
c                          = 2 if a second-order approximation is sought
c                          = 4 if a fourth-order approximation is sought
c
c                        a,b
c                          the range of the x-independent variable;
c                          i.e., x is greater than or equal to a and
c                          less than or equal to b.  a must be less than
c                          b.
c
c                        m
c                          the number of panels into which the interval
c                          [a,b] is subdivided.  hence, there will be
c                          m+1 grid points in the x-direction given by
c                          xi=a+(i-1)*dlx for i=1,2,...,m+1 where
c                          dlx=(b-a)/m is the panel width.  m must be
c                          less than idmn and greater than 5.
c
c                        mbdcnd
c                          indicates the type of boundary condition at
c                          x=a and x=b
c                          = 0 if the solution is periodic in x; i.e.,
c                              u(x+b-a,y)=u(x,y) for all y,x
c                          = 1 if the solution is specified at x=a and
c                              x=b; i.e., u(a,y) and u(b,y) are
c                              specified for all y
c                          = 2 if the solution is specified at x=a and
c                              the boundary condition is mixed at x=b;
c                              i.e., u(a,y) and du(b,y)/dx+beta*u(b,y)
c                              are specified for all y
c                          = 3 if the boundary conditions at x=a and x=b
c                              are mixed; i.e., du(a,y)/dx+alpha*u(a,y)
c                              and du(b,y)/dx+beta*u(b,y) are specified
c                              for all y
c                          = 4 if the boundary condition at x=a is mixed
c                              and the solution is specified at x=b;
c                              i.e., du(a,y)/dx+alpha*u(a,y) and u(b,y)
c                              are specified for all y
c
c                        bda
c                          a one-dimensional array of length n+1 that
c                          specifies the values of du(a,y)/dx+
c                          alpha*u(a,y) at x=a, when mbdcnd=3 or 4.
c                               bda(j) = du(a,yj)/dx+alpha*u(a,yj);
c                               j=1,2,...,n+1
c                          when mbdcnd has any other value, bda is a
c                          dummy parameter.
c
c on input               alpha
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at x=a (see
c                          argument bda).  if mbdcnd = 3,4 then alpha is
c                          a dummy parameter.
c
c                        bdb
c                          a one-dimensional array of length n+1 that
c                          specifies the values of du(b,y)/dx+
c                          beta*u(b,y) at x=b.  when mbdcnd=2 or 3
c                               bdb(j) = du(b,yj)/dx+beta*u(b,yj);
c                               j=1,2,...,n+1
c                          when mbdcnd has any other value, bdb is a
c                          dummy parameter.
c
c                        beta
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at x=b (see
c                          argument bdb).  if mbdcnd=2,3 then beta is a
c                          dummy parameter.
c
c                        c,d
c                          the range of the y-independent variable;
c                          i.e., y is greater than or equal to c and
c                          less than or equal to d.  c must be less than
c                          d.
c
c                        n
c                          the number of panels into which the interval
c                          [c,d] is subdivided.  hence, there will be
c                          n+1 grid points in the y-direction given by
c                          yj=c+(j-1)*dly for j=1,2,...,n+1 where
c                          dly=(d-c)/n is the panel width.  in addition,
c                          n must be greater than 4.
c
c                        nbdcnd
c                          indicates the types of boundary conditions at
c                          y=c and y=d
c                          = 0 if the solution is periodic in y; i.e.,
c                              u(x,y+d-c)=u(x,y) for all x,y
c                          = 1 if the solution is specified at y=c and
c                              y = d, i.e., u(x,c) and u(x,d) are
c                              specified for all x
c                          = 2 if the solution is specified at y=c and
c                              the boundary condition is mixed at y=d;
c                              i.e., u(x,c) and du(x,d)/dy+xnu*u(x,d)
c                              are specified for all x
c                          = 3 if the boundary conditions are mixed at
c                              y=c and y=d; i.e., du(x,d)/dy+gama*u(x,c)
c                              and du(x,d)/dy+xnu*u(x,d) are specified
c                              for all x
c                          = 4 if the boundary condition is mixed at y=c
c                              and the solution is specified at y=d;
c                              i.e. du(x,c)/dy+gama*u(x,c) and u(x,d)
c                              are specified for all x
c
c                        bdc
c                          a one-dimensional array of length m+1 that
c                          specifies the value of du(x,c)/dy+gama*u(x,c)
c                          at y=c.  when nbdcnd=3 or 4
c                             bdc(i) = du(xi,c)/dy + gama*u(xi,c);
c                             i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdc is a
c                          dummy parameter.
c
c                        gama
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at y=c (see
c                          argument bdc).  if nbdcnd=3,4 then gama is a
c                          dummy parameter.
c
c                        bdd
c                          a one-dimensional array of length m+1 that
c                          specifies the value of du(x,d)/dy +
c                          xnu*u(x,d) at y=c.  when nbdcnd=2 or 3
c                            bdd(i) = du(xi,d)/dy + xnu*u(xi,d);
c                            i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdd is a
c                          dummy parameter.
c
c                        xnu
c                          the scalar multiplying the solution in case
c                          of a mixed boundary condition at y=d (see
c                          argument bdd).  if nbdcnd=2 or 3 then xnu is
c                          a dummy parameter.
c
c                        cofx
c                          a user-supplied subprogram with
c                          parameters x, afun, bfun, cfun which
c                          returns the values of the x-dependent
c                          coefficients af(x), bf(x), cf(x) in
c                          the elliptic equation at x.
c
c                        cofy
c                          a user-supplied subprogram with
c                          parameters y, dfun, efun, ffun which
c                          returns the values of the y-dependent
c                          coefficients df(y), ef(y), ff(y) in
c                          the elliptic equation at y.
c
c                        note:  cofx and cofy must be declared external
c                        in the calling routine.  the values returned in
c                        afun and dfun must satisfy afun*dfun greater
c                        than 0 for a less than x less than b,
c                        c less than y less than d (see ierror=10).
c                        the coefficients provided may lead to a matrix
c                        equation which is not diagonally dominant in
c                        which case solution may fail (see ierror=4).
c
c                        grhs
c                          a two-dimensional array that specifies the
c                          values of the right-hand side of the elliptic
c                          equation; i.e., grhs(i,j)=g(xi,yi), for
c                          i=2,...,m; j=2,...,n.  at the boundaries,
c                          grhs is defined by
c
c                          mbdcnd   grhs(1,j)   grhs(m+1,j)
c                          ------   ---------   -----------
c                            0      g(a,yj)     g(b,yj)
c                            1         *           *
c                            2         *        g(b,yj)  j=1,2,...,n+1
c                            3      g(a,yj)     g(b,yj)
c                            4      g(a,yj)        *
c
c                          nbdcnd   grhs(i,1)   grhs(i,n+1)
c                          ------   ---------   -----------
c                            0      g(xi,c)     g(xi,d)
c                            1         *           *
c                            2         *        g(xi,d)  i=1,2,...,m+1
c                            3      g(xi,c)     g(xi,d)
c                            4      g(xi,c)        *
c
c                          where * means these quantities are not used.
c                          grhs should be dimensioned idmn by at least
c                          n+1 in the calling routine.
c
c                        usol
c                          a two-dimensional array that specifies the
c                          values of the solution along the boundaries.
c                          at the boundaries, usol is defined by
c
c                          mbdcnd   usol(1,j)   usol(m+1,j)
c                          ------   ---------   -----------
c                            0         *           *
c                            1      u(a,yj)     u(b,yj)
c                            2      u(a,yj)        *     j=1,2,...,n+1
c                            3         *           *
c                            4         *        u(b,yj)
c
c                          nbdcnd   usol(i,1)   usol(i,n+1)
c                          ------   ---------   -----------
c                            0         *           *
c                            1      u(xi,c)     u(xi,d)
c                            2      u(xi,c)        *     i=1,2,...,m+1
c                            3         *           *
c                            4         *        u(xi,d)
c
c                          where * means the quantities are not used in
c                          the solution.
c
c                          if iorder=2, the user may equivalence grhs
c                          and usol to save space.  note that in this
c                          case the tables specifying the boundaries of
c                          the grhs and usol arrays determine the
c                          boundaries uniquely except at the corners.
c                          if the tables call for both g(x,y) and
c                          u(x,y) at a corner then the solution must be
c                          chosen.  for example, if mbdcnd=2 and
c                          nbdcnd=4, then u(a,c), u(a,d), u(b,d) must be
c                          chosen at the corners in addition to g(b,c).
c
c                          if iorder=4, then the two arrays, usol and
c                          grhs, must be distinct.
c
c                          usol should be dimensioned idmn by at least
c                          n+1 in the calling routine.
c
c                        idmn
c                          the row (or first) dimension of the arrays
c                          grhs and usol as it appears in the program
c                          calling sepeli.  this parameter is used to
c                          specify the variable dimension of grhs and
c                          usol.  idmn must be at least 7 and greater
c                          than or equal to m+1.
c
c                        w
c                          a one-dimensional array that must be provided
c                          by the user for work space.  let
c                          k=int(log2(n+1))+1 and set  l=2**(k+1).
c                          then (k-2)*l+k+10*n+12*m+27 will suffice
c                          as a length of w.  the actual length of w in
c                          the calling routine must be set in w(1) (see
c                          ierror=11).
c
c on output              usol
c                          contains the approximate solution to the
c                          elliptic equation.  usol(i,j) is the
c                          approximation to u(xi,yj) for i=1,2...,m+1
c                          and j=1,2,...,n+1.  the approximation has
c                          error o(dlx**2+dly**2) if called with
c                          iorder=2 and o(dlx**4+dly**4) if called with
c                          iorder=4.
c
c                        w
c                          contains intermediate values that must not be
c                          destroyed if sepeli is called again with
c                          intl=1.  in addition w(1) contains the exact
c                          minimal length (in floating point) required
c                          for the work space (see ierror=11).
c
c                        pertrb
c                          if a combination of periodic or derivative
c                          boundary conditions (i.e., alpha=beta=0 if
c                          mbdcnd=3; gama=xnu=0 if nbdcnd=3) is
c                          specified and if the coefficients of u(x,y)
c                          in the separable elliptic equation are zero
c                          (i.e., cf(x)=0 for x greater than or equal to
c                          a and less than or equal to b; ff(y)=0 for
c                          y greater than or equal to c and less than
c                          or equal to d) then a solution may not exist.
c                          pertrb is a constant calculated and
c                          subtracted from the right-hand side of the
c                          matrix equations generated by sepeli which
c                          insures that a solution exists.  sepeli then
c                          computes this solution which is a weighted
c                          minimal least squares solution to the
c                          original problem.
c
c                        ierror
c                          an error flag that indicates invalid input
c                          parameters or failure to find a solution
c                          = 0 no error
c                          = 1 if a greater than b or c greater than d
c                          = 2 if mbdcnd less than 0 or mbdcnd greater
c                              than 4
c                          = 3 if nbdcnd less than 0 or nbdcnd greater
c                              than 4
c                          = 4 if attempt to find a solution fails.
c                              (the linear system generated is not
c                              diagonally dominant.)
c                          = 5 if idmn is too small (see discussion of
c                              idmn)
c                          = 6 if m is too small or too large (see
c                              discussion of m)
c                          = 7 if n is too small (see discussion of n)
c                          = 8 if iorder is not 2 or 4
c                          = 9 if intl is not 0 or 1
c                          = 10 if afun*dfun less than or equal to 0 for
c                               some interior mesh point (xi,yj)
c                          = 11 if the work space length input in w(1)
c                               is less than the exact minimal work
c                               space length required output in w(1).
c
c                          note (concerning ierror=4):  for the
c                          coefficients input through cofx, cofy, the
c                          discretization may lead to a block
c                          tridiagonal linear system which is not
c                          diagonally dominant (for example, this
c                          happens if cfun=0 and bfun/(2.*dlx) greater
c                          than afun/dlx**2).  in this case solution may
c                          fail.  this cannot happen in the limit as
c                          dlx, dly approach zero.  hence, the condition
c                          may be remedied by taking larger values for m
c                          or n.
c
c entry points           sepeli, spelip, chkprm, chksng, orthog, minsol,
c                        trisp, defer, dx, dy, blktri, blktr1, indxb,
c                        indxa, indxc, prod, prodp, cprod, cprodp,
c                        ppadd, psgf, bsrh, ppsgf, ppspf, compb,
c                        trun1, stor1, tqlrat
c
c special conditions     none
c
c common blocks          splp, cblkt
c
c i/o                    none
c
c precision              single
c
c specialist             john c. adams, ncar, boulder, colorado  80307
c
c language               fortran
c
c history                developed at ncar during 1975-76.
c
c algorithm              sepeli automatically discretizes the separable
c                        elliptic equation which is then solved by a
c                        generalized cyclic reduction algorithm in the
c                        subroutine, blktri.  the fourth-order solution
c                        is obtained using 'deferred corrections' which
c                        is described and referenced in sections,
c                        references and method.
c
c space required         14654 (octal) = 6572 (decimal)
c
c accuracy and timing    the following computational results were
c                        obtained by solving the sample problem at the
c                        end of this write-up on the control data 7600.
c                        the op count is proportional to m*n*log2(n).
c                        in contrast to the other routines in this
c                        chapter, accuracy is tested by computing and
c                        tabulating second- and fourth-order
c                        discretization errors.  below is a table
c                        containing computational results.  the times
c                        given do not include initialization (i.e.,
c                        times are for intl=1).  note that the
c                        fourth-order accuracy is not realized until the
c                        mesh is sufficiently refined.
c
c              second-order    fourth-order   second-order  fourth-order
c    m    n   execution time  execution time    error         error
c               (m sec)         (m sec)
c     6    6         6              14          6.8e-1        1.2e0
c    14   14        23              58          1.4e-1        1.8e-1
c    30   30       100             247          3.2e-2        9.7e-3
c    62   62       445           1,091          7.5e-3        3.0e-4
c   126  126     2,002           4,772          1.8e-3        3.5e-6
c
c portability            there are no machine-dependent constants.
c
c required resident      sqrt, abs, log
c routines
c
c references             keller, h.b., 'numerical methods for two-point
c                          boundary-value problems', blaisdel (1968),
c                          waltham, mass.
c
c                        swarztrauber, p., and r. sweet (1975):
c                          'efficient fortran subprograms for the
c                          solution of elliptic partial differential
c                          equations'.  ncar technical note
c                          ncar-tn/ia-109, pp. 135-137.
c
c***references  h. b. keller, numerical methods for two-point
c                 boundary-value problems, blaisdel, waltham, mass.,
c                 1968.
c               p. n. swarztrauber and r. sweet, efficient fortran
c                 subprograms for the solution of elliptic equations,
c                 ncar tn/ia-109, july 1975, 138 pp.
c***routines called  chkprm, spelip
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sepeli
c
      dimension       grhs(idmn,*)           ,usol(idmn,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     1                w(*)
      external        cofx       ,cofy
c***first executable statement  sepeli
      call chkprm (intl,iorder,a,b,m,mbdcnd,c,d,n,nbdcnd,cofx,cofy,
     1             idmn,ierror)
      if (ierror .ne. 0) return
c
c     compute minimum work space and check work space length input
c
      l = n+1
      if (nbdcnd .eq. 0) l = n
      logb2n = int(log(l+0.5)/log(2.0))+1
      ll = 2**(logb2n+1)
      k = m+1
      l = n+1
      length = (logb2n-2)*ll+logb2n+max(2*l,6*k)+5
      if (nbdcnd .eq. 0) length = length+2*l
      ierror = 11
      linput = int(w(1)+0.5)
      loutpt = length+6*(k+l)+1
      w(1) = loutpt
      if (loutpt .gt. linput) return
      ierror = 0
c
c     set work space indices
c
      i1 = length+2
      i2 = i1+l
      i3 = i2+l
      i4 = i3+l
      i5 = i4+l
      i6 = i5+l
      i7 = i6+l
      i8 = i7+k
      i9 = i8+k
      i10 = i9+k
      i11 = i10+k
      i12 = i11+k
      i13 = 2
      call spelip (intl,iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,d,n,
     1             nbdcnd,bdc,gama,bdd,xnu,cofx,cofy,w(i1),w(i2),w(i3),
     2             w(i4),w(i5),w(i6),w(i7),w(i8),w(i9),w(i10),w(i11),
     3             w(i12),grhs,usol,idmn,w(i13),pertrb,ierror)
      return
      end
