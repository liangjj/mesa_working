*deck sepx4
      subroutine sepx4 (iorder, a, b, m, mbdcnd, bda, alpha, bdb, beta,
     +   c, d, n, nbdcnd, bdc, bdd, cofx, grhs, usol, idmn, w, pertrb,
     +   ierror)
c***begin prologue  sepx4
c***purpose  solve for either the second or fourth order finite
c            difference approximation to the solution of a separable
c            elliptic partial differential equation on a rectangle.
c            any combination of periodic or mixed boundary conditions is
c            allowed.
c***library   slatec (fishpack)
c***category  i2b1a2
c***type      single precision (sepx4-s)
c***keywords  elliptic, fishpack, helmholtz, pde, separable
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c purpose                sepx4 solves for either the second-order
c                        finite difference approximation or a
c                        fourth-order approximation  to the
c                        solution of a separable elliptic equation
c                             af(x)*uxx+bf(x)*ux+cf(x)*u+uyy = g(x,y)
c
c                        on a rectangle (x greater than or equal to a
c                        and less than or equal to b; y greater than
c                        or equal to c and less than or equal to d).
c                        any combination of periodic or mixed boundary
c                        conditions is allowed.
c                        if boundary conditions in the x direction
c                        are periodic (see mbdcnd=0 below) then the
c                        coefficients must satisfy
c                        af(x)=c1,bf(x)=0,cf(x)=c2 for all x.
c                        here c1,c2 are constants, c1.gt.0.
c
c                        the possible boundary conditions are
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
c                         (2) u(x,c),du(x,d)/dy are specified for all x
c                         (3) du(x,c)/dy,du(x,d)/dy are specified for
c                            all x
c                        (4) du(x,c)/dy,u(x,d) are specified for all x
c
c usage                  call sepx4(iorder,a,b,m,mbdcnd,bda,alpha,bdb,
c                                  beta,c,d,n,nbdcnd,bdc,bdd,cofx,
c                                  grhs,usol,idmn,w,pertrb,ierror)
c
c arguments
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
c                              i.e., du(x,c)/dy and u(x,d)
c                              are specified for all x
c                          = 3 if the boundary conditions are mixed at
c                              y= c and y=d i.e., du(x,d)/dy
c                              and du(x,d)/dy are specified
c                              for all x
c                          = 4 if the boundary condition is mixed at y=c
c                              and the solution is specified at y=d;
c                              i.e. du(x,c)/dy+gama*u(x,c) and u(x,d)
c                              are specified for all x
c
c                        bdc
c                          a one-dimensional array of length m+1 that
c                          specifies the value du(x,c)/dy
c                          at y=c.  when nbdcnd=3 or 4
c                            bdc(i) = du(xi,c)/dy
c                             i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdc is a
c                          dummy parameter.
c
c
c                        bdd
c                          a one-dimensional array of length m+1 that
c                          specifies the value of du(x,d)/dy
c                          at y=d.  when nbdcnd=2 or 3
c                            bdd(i)=du(xi,d)/dy
c                             i=1,2,...,m+1.
c                          when nbdcnd has any other value, bdd is a
c                          dummy parameter.
c
c
c                        cofx
c                          a user-supplied subprogram with
c                          parameters x, afun, bfun, cfun which
c                          returns the values of the x-dependent
c                          coefficients af(x), bf(x), cf(x) in
c                          the elliptic equation at x.
c                          if boundary conditions in the x direction
c                          are periodic then the coefficients
c                          must satisfy af(x)=c1,bf(x)=0,cf(x)=c2 for
c                          all x.  here c1.gt.0 and c2 are constants.
c
c                          note that cofx must be declared external
c                          in the calling routine.
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
c                          calling sepx4.  this parameter is used to
c                          specify the variable dimension of grhs and
c                          usol.  idmn must be at least 7 and greater
c                          than or equal to m+1.
c
c                        w
c                          a one-dimensional array that must be provided
c                          by the user for work space.
c                          10*n+(16+int(log2(n)))*(m+1)+23 will suffice
c                          as a length for w.  the actual length of
c                          w in the calling routine must be set in w(1)
c                          (see ierror=11).
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
c                          w(1) contains the exact minimal length (in
c                          floating point) required for the work space
c                          (see ierror=11).
c
c                        pertrb
c                          if a combination of periodic or derivative
c                          boundary conditions (i.e., alpha=beta=0 if
c                          mbdcnd=3) is specified and if cf(x)=0 for all
c                          x, then a solution to the discretized matrix
c                          equation may not exist (reflecting the non-
c                          uniqueness of solutions to the pde).  pertrb
c                          is a constant calculated and subtracted from
c                          the right hand side of the matrix equation
c                          insuring the existence of a solution.
c                          sepx4 computes this solution which is a
c                          weighted minimal least squares solution to
c                          the original problem.  if singularity is
c                          not detected pertrb=0.0 is returned by
c                          sepx4.
c
c                        ierror
c                          an error flag that indicates invalid input
c                          parameters or failure to find a solution
c                          = 0  no error
c                          = 1  if a greater than b or c greater than d
c                          = 2  if mbdcnd less than 0 or mbdcnd greater
c                               than 4
c                          = 3  if nbdcnd less than 0 or nbdcnd greater
c                               than 4
c                          = 4  if attempt to find a solution fails.
c                               (the linear system generated is not
c                               diagonally dominant.)
c                          = 5  if idmn is too small (see discussion of
c                               idmn)
c                          = 6  if m is too small or too large (see
c                               discussion of m)
c                          = 7  if n is too small (see discussion of n)
c                          = 8  if iorder is not 2 or 4
c                          = 10 if afun is less than or equal to zero
c                               for some interior mesh point xi
c                          = 11 if the work space length input in w(1)
c                               is less than the exact minimal work
c                               space length required output in w(1).
c                          = 12 if mbdcnd=0 and af(x)=cf(x)=constant
c                               or bf(x)=0 for all x is not true.
c
c *long description:
c
c dimension of           bda(n+1), bdb(n+1), bdc(m+1), bdd(m+1),
c arguments              usol(idmn,n+1), grhs(idmn,n+1),
c                        w (see argument list)
c
c latest revision        october 1980
c
c special conditions     none
c
c common blocks          spl4
c
c i/o                    none
c
c precision              single
c
c required library       none
c files
c
c specialist             john c. adams, ncar, boulder, colorado  80307
c
c language               fortran
c
c
c entry points           sepx4,speli4,chkpr4,chksn4,ortho4,minso4,tris4,
c                        defe4,dx4,dy4
c
c history                sepx4 was developed by modifying the ulib
c                        routine sepeli during october 1978.
c                        it should be used instead of sepeli whenever
c                        possible.  the increase in speed is at least
c                        a factor of three.
c
c algorithm              sepx4 automatically discretizes the separable
c                        elliptic equation which is then solved by a
c                        generalized cyclic reduction algorithm in the
c                        subroutine pois.  the fourth order solution
c                        is obtained using the technique of
c                        deferred corrections referenced below.
c
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
c***routines called  chkpr4, speli4
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c   920122  minor corrections and modifications to prologue.  (wrb)
c   920501  reformatted the references section.  (wrb)
c***end prologue  sepx4
c
      dimension       grhs(idmn,*)           ,usol(idmn,*)
      dimension       bda(*)     ,bdb(*)     ,bdc(*)     ,bdd(*)     ,
     1                w(*)
      external cofx
c***first executable statement  sepx4
      call chkpr4(iorder,a,b,m,mbdcnd,c,d,n,nbdcnd,cofx,idmn,ierror)
      if (ierror .ne. 0) return
c
c     compute minimum work space and check work space length input
c
      l = n+1
      if (nbdcnd .eq. 0) l = n
      k = m+1
      l = n+1
c     estimate log base 2 of n
      log2n=int(log(real(n+1))/log(2.0)+0.5)
      length=4*(n+1)+(10+log2n)*(m+1)
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
      call speli4(iorder,a,b,m,mbdcnd,bda,alpha,bdb,beta,c,d,n,
     1nbdcnd,bdc,bdd,cofx,w(i1),w(i2),w(i3),
     2             w(i4),w(i5),w(i6),w(i7),w(i8),w(i9),w(i10),w(i11),
     3             w(i12),grhs,usol,idmn,w(i13),pertrb,ierror)
      return
      end
