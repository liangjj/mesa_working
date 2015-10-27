*deck hw3crt
      subroutine hw3crt (xs, xf, l, lbdcnd, bdxs, bdxf, ys, yf, m,
     +   mbdcnd, bdys, bdyf, zs, zf, n, nbdcnd, bdzs, bdzf, elmbda,
     +   ldimf, mdimf, f, pertrb, ierror, w)
c***begin prologue  hw3crt
c***purpose  solve the standard seven-point finite difference
c            approximation to the helmholtz equation in cartesian
c            coordinates.
c***library   slatec (fishpack)
c***category  i2b1a1a
c***type      single precision (hw3crt-s)
c***keywords  cartesian, elliptic, fishpack, helmholtz, pde
c***author  adams, j., (ncar)
c           swarztrauber, p. n., (ncar)
c           sweet, r., (ncar)
c***description
c
c     subroutine hw3crt solves the standard seven-point finite
c     difference approximation to the helmholtz equation in cartesian
c     coordinates:
c
c         (d/dx)(du/dx) + (d/dy)(du/dy) + (d/dz)(du/dz)
c
c                    + lambda*u = f(x,y,z) .
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c    * * * * * * * *    parameter description     * * * * * * * * * *
c
c
c            * * * * * *   on input    * * * * * *
c
c     xs,xf
c        the range of x, i.e. xs .le. x .le. xf .
c        xs must be less than xf.
c
c     l
c        the number of panels into which the interval (xs,xf) is
c        subdivided.  hence, there will be l+1 grid points in the
c        x-direction given by x(i) = xs+(i-1)dx for i=1,2,...,l+1,
c        where dx = (xf-xs)/l is the panel width.  l must be at
c        least 5 .
c
c     lbdcnd
c        indicates the type of boundary conditions at x = xs and x = xf.
c
c        = 0  if the solution is periodic in x, i.e.
c             u(l+i,j,k) = u(i,j,k).
c        = 1  if the solution is specified at x = xs and x = xf.
c        = 2  if the solution is specified at x = xs and the derivative
c             of the solution with respect to x is specified at x = xf.
c        = 3  if the derivative of the solution with respect to x is
c             specified at x = xs and x = xf.
c        = 4  if the derivative of the solution with respect to x is
c             specified at x = xs and the solution is specified at x=xf.
c
c     bdxs
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to x at x = xs.
c        when lbdcnd = 3 or 4,
c
c             bdxs(j,k) = (d/dx)u(xs,y(j),z(k)), j=1,2,...,m+1,
c                                                k=1,2,...,n+1.
c
c        when lbdcnd has any other value, bdxs is a dummy variable.
c        bdxs must be dimensioned at least (m+1)*(n+1).
c
c     bdxf
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to x at x = xf.
c        when lbdcnd = 2 or 3,
c
c             bdxf(j,k) = (d/dx)u(xf,y(j),z(k)), j=1,2,...,m+1,
c                                                k=1,2,...,n+1.
c
c        when lbdcnd has any other value, bdxf is a dummy variable.
c        bdxf must be dimensioned at least (m+1)*(n+1).
c
c     ys,yf
c        the range of y, i.e. ys .le. y .le. yf.
c        ys must be less than yf.
c
c     m
c        the number of panels into which the interval (ys,yf) is
c        subdivided.  hence, there will be m+1 grid points in the
c        y-direction given by y(j) = ys+(j-1)dy for j=1,2,...,m+1,
c        where dy = (yf-ys)/m is the panel width.  m must be at
c        least 5 .
c
c     mbdcnd
c        indicates the type of boundary conditions at y = ys and y = yf.
c
c        = 0  if the solution is periodic in y, i.e.
c             u(i,m+j,k) = u(i,j,k).
c        = 1  if the solution is specified at y = ys and y = yf.
c        = 2  if the solution is specified at y = ys and the derivative
c             of the solution with respect to y is specified at y = yf.
c        = 3  if the derivative of the solution with respect to y is
c             specified at y = ys and y = yf.
c        = 4  if the derivative of the solution with respect to y is
c             specified at y = ys and the solution is specified at y=yf.
c
c     bdys
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to y at y = ys.
c        when mbdcnd = 3 or 4,
c
c             bdys(i,k) = (d/dy)u(x(i),ys,z(k)), i=1,2,...,l+1,
c                                                k=1,2,...,n+1.
c
c        when mbdcnd has any other value, bdys is a dummy variable.
c        bdys must be dimensioned at least (l+1)*(n+1).
c
c     bdyf
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to y at y = yf.
c        when mbdcnd = 2 or 3,
c
c             bdyf(i,k) = (d/dy)u(x(i),yf,z(k)), i=1,2,...,l+1,
c                                                k=1,2,...,n+1.
c
c        when mbdcnd has any other value, bdyf is a dummy variable.
c        bdyf must be dimensioned at least (l+1)*(n+1).
c
c     zs,zf
c        the range of z, i.e. zs .le. z .le. zf.
c        zs must be less than zf.
c
c     n
c        the number of panels into which the interval (zs,zf) is
c        subdivided.  hence, there will be n+1 grid points in the
c        z-direction given by z(k) = zs+(k-1)dz for k=1,2,...,n+1,
c        where dz = (zf-zs)/n is the panel width.  n must be at least 5.
c
c     nbdcnd
c        indicates the type of boundary conditions at z = zs and z = zf.
c
c        = 0  if the solution is periodic in z, i.e.
c             u(i,j,n+k) = u(i,j,k).
c        = 1  if the solution is specified at z = zs and z = zf.
c        = 2  if the solution is specified at z = zs and the derivative
c             of the solution with respect to z is specified at z = zf.
c        = 3  if the derivative of the solution with respect to z is
c             specified at z = zs and z = zf.
c        = 4  if the derivative of the solution with respect to z is
c             specified at z = zs and the solution is specified at z=zf.
c
c     bdzs
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to z at z = zs.
c        when nbdcnd = 3 or 4,
c
c             bdzs(i,j) = (d/dz)u(x(i),y(j),zs), i=1,2,...,l+1,
c                                                j=1,2,...,m+1.
c
c        when nbdcnd has any other value, bdzs is a dummy variable.
c        bdzs must be dimensioned at least (l+1)*(m+1).
c
c     bdzf
c        a two-dimensional array that specifies the values of the
c        derivative of the solution with respect to z at z = zf.
c        when nbdcnd = 2 or 3,
c
c             bdzf(i,j) = (d/dz)u(x(i),y(j),zf), i=1,2,...,l+1,
c                                                j=1,2,...,m+1.
c
c        when nbdcnd has any other value, bdzf is a dummy variable.
c        bdzf must be dimensioned at least (l+1)*(m+1).
c
c     elmbda
c        the constant lambda in the helmholtz equation. if
c        lambda .gt. 0, a solution may not exist.  however, hw3crt will
c        attempt to find a solution.
c
c     f
c        a three-dimensional array that specifies the values of the
c        right side of the helmholtz equation and boundary values (if
c        any).  for i=2,3,...,l, j=2,3,...,m, and k=2,3,...,n
c
c                   f(i,j,k) = f(x(i),y(j),z(k)).
c
c        on the boundaries f is defined by
c
c        lbdcnd      f(1,j,k)         f(l+1,j,k)
c        ------   ---------------   ---------------
c
c          0      f(xs,y(j),z(k))   f(xs,y(j),z(k))
c          1      u(xs,y(j),z(k))   u(xf,y(j),z(k))
c          2      u(xs,y(j),z(k))   f(xf,y(j),z(k))   j=1,2,...,m+1
c          3      f(xs,y(j),z(k))   f(xf,y(j),z(k))   k=1,2,...,n+1
c          4      f(xs,y(j),z(k))   u(xf,y(j),z(k))
c
c        mbdcnd      f(i,1,k)         f(i,m+1,k)
c        ------   ---------------   ---------------
c
c          0      f(x(i),ys,z(k))   f(x(i),ys,z(k))
c          1      u(x(i),ys,z(k))   u(x(i),yf,z(k))
c          2      u(x(i),ys,z(k))   f(x(i),yf,z(k))   i=1,2,...,l+1
c          3      f(x(i),ys,z(k))   f(x(i),yf,z(k))   k=1,2,...,n+1
c          4      f(x(i),ys,z(k))   u(x(i),yf,z(k))
c
c        nbdcnd      f(i,j,1)         f(i,j,n+1)
c        ------   ---------------   ---------------
c
c          0      f(x(i),y(j),zs)   f(x(i),y(j),zs)
c          1      u(x(i),y(j),zs)   u(x(i),y(j),zf)
c          2      u(x(i),y(j),zs)   f(x(i),y(j),zf)   i=1,2,...,l+1
c          3      f(x(i),y(j),zs)   f(x(i),y(j),zf)   j=1,2,...,m+1
c          4      f(x(i),y(j),zs)   u(x(i),y(j),zf)
c
c        f must be dimensioned at least (l+1)*(m+1)*(n+1).
c
c        note:
c
c        if the table calls for both the solution u and the right side f
c        on a boundary, then the solution must be specified.
c
c     ldimf
c        the row (or first) dimension of the arrays f,bdys,bdyf,bdzs,
c        and bdzf as it appears in the program calling hw3crt. this
c        parameter is used to specify the variable dimension of these
c        arrays.  ldimf must be at least l+1.
c
c     mdimf
c        the column (or second) dimension of the array f and the row (or
c        first) dimension of the arrays bdxs and bdxf as it appears in
c        the program calling hw3crt.  this parameter is used to specify
c        the variable dimension of these arrays.
c        mdimf must be at least m+1.
c
c     w
c        a one-dimensional array that must be provided by the user for
c        work space.  the length of w must be at least 30 + l + m + 5*n
c        + max(l,m,n) + 7*(int((l+1)/2) + int((m+1)/2))
c
c
c            * * * * * *   on output   * * * * * *
c
c     f
c        contains the solution u(i,j,k) of the finite difference
c        approximation for the grid point (x(i),y(j),z(k)) for
c        i=1,2,...,l+1, j=1,2,...,m+1, and k=1,2,...,n+1.
c
c     pertrb
c        if a combination of periodic or derivative boundary conditions
c        is specified for a poisson equation (lambda = 0), a solution
c        may not exist.  pertrb is a constant, calculated and subtracted
c        from f, which ensures that a solution exists.  pwscrt then
c        computes this solution, which is a least squares solution to
c        the original approximation.  this solution is not unique and is
c        unnormalized.  the value of pertrb should be small compared to
c        the right side f.  otherwise, a solution is obtained to an
c        essentially different problem.  this comparison should always
c        be made to insure that a meaningful solution has been obtained.
c
c     ierror
c        an error flag that indicates invalid input parameters.  except
c        for numbers 0 and 12, a solution is not attempted.
c
c        =  0  no error
c        =  1  xs .ge. xf
c        =  2  l .lt. 5
c        =  3  lbdcnd .lt. 0 .or. lbdcnd .gt. 4
c        =  4  ys .ge. yf
c        =  5  m .lt. 5
c        =  6  mbdcnd .lt. 0 .or. mbdcnd .gt. 4
c        =  7  zs .ge. zf
c        =  8  n .lt. 5
c        =  9  nbdcnd .lt. 0 .or. nbdcnd .gt. 4
c        = 10  ldimf .lt. l+1
c        = 11  mdimf .lt. m+1
c        = 12  lambda .gt. 0
c
c        since this is the only means of indicating a possibly incorrect
c        call to hw3crt, the user should test ierror after the call.
c
c *long description:
c
c    * * * * * * *   program specifications    * * * * * * * * * * * *
c
c     dimension of   bdxs(mdimf,n+1),bdxf(mdimf,n+1),bdys(ldimf,n+1),
c     arguments      bdyf(ldimf,n+1),bdzs(ldimf,m+1),bdzf(ldimf,m+1),
c                    f(ldimf,mdimf,n+1),w(see argument list)
c
c     latest         december 1, 1978
c     revision
c
c     subprograms    hw3crt,pois3d,pos3d1,tridq,rffti,rfftf,rfftf1,
c     required       rfftb,rfftb1,costi,cost,sinti,sint,cosqi,cosqf,
c                    cosqf1,cosqb,cosqb1,sinqi,sinqf,sinqb,cffti,
c                    cffti1,cfftb,cfftb1,passb2,passb3,passb4,passb,
c                    cfftf,cfftf1,passf1,passf2,passf3,passf4,passf,
c                    pimach
c
c     special        none
c     conditions
c
c     common         none
c     blocks
c
c     i/o            none
c
c     precision      single
c
c     specialist     roland sweet
c
c     language       fortran
c
c     history        written by roland sweet at ncar in july 1977
c
c     algorithm      this subroutine defines the finite difference
c                    equations, incorporates boundary data, and
c                    adjusts the right side of singular systems and
c                    then calls pois3d to solve the system.
c
c     space          7862(decimal) = 17300(octal) locations on the
c     required       ncar control data 7600
c
c     timing and        the execution time t on the ncar control data
c     accuracy       7600 for subroutine hw3crt is roughly proportional
c                    to l*m*n*(log2(l)+log2(m)+5), but also depends on
c                    input parameters lbdcnd and mbdcnd.  some typical
c                    values are listed in the table below.
c                       the solution process employed results in a loss
c                    of no more than three significant digits for l,m
c                    and n as large as 32.  more detailed information
c                    about accuracy can be found in the documentation
c                    for subroutine pois3d which is the routine that
c                    actually solves the finite difference equations.
c
c
c                       l(=m=n)     lbdcnd(=mbdcnd=nbdcnd)      t(msecs)
c                       -------     ----------------------      --------
c
c                         16                  0                    300
c                         16                  1                    302
c                         16                  3                    348
c                         32                  0                   1925
c                         32                  1                   1929
c                         32                  3                   2109
c
c     portability    american national standards institute fortran.
c                    the machine dependent constant pi is defined in
c                    function pimach.
c
c     required       cos,sin,atan
c     resident
c     routines
c
c     reference      none
c
c    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c***references  (none)
c***routines called  pois3d
c***revision history  (yymmdd)
c   801001  date written
c   890531  changed all specific intrinsics to generic.  (wrb)
c   890531  revision date from version 3.2
c   891214  prologue converted to version 4.0 format.  (bab)
c***end prologue  hw3crt
c
c
      dimension       bdxs(mdimf,*)          ,bdxf(mdimf,*)          ,
     1                bdys(ldimf,*)          ,bdyf(ldimf,*)          ,
     2                bdzs(ldimf,*)          ,bdzf(ldimf,*)          ,
     3                f(ldimf,mdimf,*)       ,w(*)
c***first executable statement  hw3crt
      ierror = 0
      if (xf .le. xs) ierror = 1
      if (l .lt. 5) ierror = 2
      if (lbdcnd.lt.0 .or. lbdcnd.gt.4) ierror = 3
      if (yf .le. ys) ierror = 4
      if (m .lt. 5) ierror = 5
      if (mbdcnd.lt.0 .or. mbdcnd.gt.4) ierror = 6
      if (zf .le. zs) ierror = 7
      if (n .lt. 5) ierror = 8
      if (nbdcnd.lt.0 .or. nbdcnd.gt.4) ierror = 9
      if (ldimf .lt. l+1) ierror = 10
      if (mdimf .lt. m+1) ierror = 11
      if (ierror .ne. 0) go to 188
      dy = (yf-ys)/m
      twbydy = 2./dy
      c2 = 1./(dy**2)
      mstart = 1
      mstop = m
      mp1 = m+1
      mp = mbdcnd+1
      go to (104,101,101,102,102),mp
  101 mstart = 2
  102 go to (104,104,103,103,104),mp
  103 mstop = mp1
  104 munk = mstop-mstart+1
      dz = (zf-zs)/n
      twbydz = 2./dz
      np = nbdcnd+1
      c3 = 1./(dz**2)
      np1 = n+1
      nstart = 1
      nstop = n
      go to (108,105,105,106,106),np
  105 nstart = 2
  106 go to (108,108,107,107,108),np
  107 nstop = np1
  108 nunk = nstop-nstart+1
      lp1 = l+1
      dx = (xf-xs)/l
      c1 = 1./(dx**2)
      twbydx = 2./dx
      lp = lbdcnd+1
      lstart = 1
      lstop = l
c
c     enter boundary data for x-boundaries.
c
      go to (122,109,109,112,112),lp
  109 lstart = 2
      do 111 j=mstart,mstop
         do 110 k=nstart,nstop
            f(2,j,k) = f(2,j,k)-c1*f(1,j,k)
  110    continue
  111 continue
      go to 115
  112 do 114 j=mstart,mstop
         do 113 k=nstart,nstop
            f(1,j,k) = f(1,j,k)+twbydx*bdxs(j,k)
  113    continue
  114 continue
  115 go to (122,116,119,119,116),lp
  116 do 118 j=mstart,mstop
         do 117 k=nstart,nstop
            f(l,j,k) = f(l,j,k)-c1*f(lp1,j,k)
  117    continue
  118 continue
      go to 122
  119 lstop = lp1
      do 121 j=mstart,mstop
         do 120 k=nstart,nstop
            f(lp1,j,k) = f(lp1,j,k)-twbydx*bdxf(j,k)
  120    continue
  121 continue
  122 lunk = lstop-lstart+1
c
c     enter boundary data for y-boundaries.
c
      go to (136,123,123,126,126),mp
  123 do 125 i=lstart,lstop
         do 124 k=nstart,nstop
            f(i,2,k) = f(i,2,k)-c2*f(i,1,k)
  124    continue
  125 continue
      go to 129
  126 do 128 i=lstart,lstop
         do 127 k=nstart,nstop
            f(i,1,k) = f(i,1,k)+twbydy*bdys(i,k)
  127    continue
  128 continue
  129 go to (136,130,133,133,130),mp
  130 do 132 i=lstart,lstop
         do 131 k=nstart,nstop
            f(i,m,k) = f(i,m,k)-c2*f(i,mp1,k)
  131    continue
  132 continue
      go to 136
  133 do 135 i=lstart,lstop
         do 134 k=nstart,nstop
            f(i,mp1,k) = f(i,mp1,k)-twbydy*bdyf(i,k)
  134    continue
  135 continue
  136 continue
c
c     enter boundary data for z-boundaries.
c
      go to (150,137,137,140,140),np
  137 do 139 i=lstart,lstop
         do 138 j=mstart,mstop
            f(i,j,2) = f(i,j,2)-c3*f(i,j,1)
  138    continue
  139 continue
      go to 143
  140 do 142 i=lstart,lstop
         do 141 j=mstart,mstop
            f(i,j,1) = f(i,j,1)+twbydz*bdzs(i,j)
  141    continue
  142 continue
  143 go to (150,144,147,147,144),np
  144 do 146 i=lstart,lstop
         do 145 j=mstart,mstop
            f(i,j,n) = f(i,j,n)-c3*f(i,j,np1)
  145    continue
  146 continue
      go to 150
  147 do 149 i=lstart,lstop
         do 148 j=mstart,mstop
            f(i,j,np1) = f(i,j,np1)-twbydz*bdzf(i,j)
  148    continue
  149 continue
c
c     define a,b,c coefficients in w-array.
c
  150 continue
      iwb = nunk+1
      iwc = iwb+nunk
      iww = iwc+nunk
      do 151 k=1,nunk
         i = iwc+k-1
         w(k) = c3
         w(i) = c3
         i = iwb+k-1
         w(i) = -2.*c3+elmbda
  151 continue
      go to (155,155,153,152,152),np
  152 w(iwc) = 2.*c3
  153 go to (155,155,154,154,155),np
  154 w(iwb-1) = 2.*c3
  155 continue
      pertrb = 0.
c
c     for singular problems adjust data to insure a solution will exist.
c
      go to (156,172,172,156,172),lp
  156 go to (157,172,172,157,172),mp
  157 go to (158,172,172,158,172),np
  158 if (elmbda) 172,160,159
  159 ierror = 12
      go to 172
  160 continue
      mstpm1 = mstop-1
      lstpm1 = lstop-1
      nstpm1 = nstop-1
      xlp = (2+lp)/3
      ylp = (2+mp)/3
      zlp = (2+np)/3
      s1 = 0.
      do 164 k=2,nstpm1
         do 162 j=2,mstpm1
            do 161 i=2,lstpm1
               s1 = s1+f(i,j,k)
  161       continue
            s1 = s1+(f(1,j,k)+f(lstop,j,k))/xlp
  162    continue
         s2 = 0.
         do 163 i=2,lstpm1
            s2 = s2+f(i,1,k)+f(i,mstop,k)
  163    continue
         s2 = (s2+(f(1,1,k)+f(1,mstop,k)+f(lstop,1,k)+f(lstop,mstop,k))/
     1                                                          xlp)/ylp
         s1 = s1+s2
  164 continue
      s = (f(1,1,1)+f(lstop,1,1)+f(1,1,nstop)+f(lstop,1,nstop)+
     1    f(1,mstop,1)+f(lstop,mstop,1)+f(1,mstop,nstop)+
     2                                   f(lstop,mstop,nstop))/(xlp*ylp)
      do 166 j=2,mstpm1
         do 165 i=2,lstpm1
            s = s+f(i,j,1)+f(i,j,nstop)
  165    continue
  166 continue
      s2 = 0.
      do 167 i=2,lstpm1
         s2 = s2+f(i,1,1)+f(i,1,nstop)+f(i,mstop,1)+f(i,mstop,nstop)
  167 continue
      s = s2/ylp+s
      s2 = 0.
      do 168 j=2,mstpm1
         s2 = s2+f(1,j,1)+f(1,j,nstop)+f(lstop,j,1)+f(lstop,j,nstop)
  168 continue
      s = s2/xlp+s
      pertrb = (s/zlp+s1)/((lunk+1.-xlp)*(munk+1.-ylp)*
     1                                              (nunk+1.-zlp))
      do 171 i=1,lunk
         do 170 j=1,munk
            do 169 k=1,nunk
               f(i,j,k) = f(i,j,k)-pertrb
  169       continue
  170    continue
  171 continue
  172 continue
      nperod = 0
      if (nbdcnd .eq. 0) go to 173
      nperod = 1
      w(1) = 0.
      w(iww-1) = 0.
  173 continue
      call pois3d (lbdcnd,lunk,c1,mbdcnd,munk,c2,nperod,nunk,w,w(iwb),
     1             w(iwc),ldimf,mdimf,f(lstart,mstart,nstart),ir,w(iww))
c
c     fill in sides for periodic boundary conditions.
c
      if (lp .ne. 1) go to 180
      if (mp .ne. 1) go to 175
      do 174 k=nstart,nstop
         f(1,mp1,k) = f(1,1,k)
  174 continue
      mstop = mp1
  175 if (np .ne. 1) go to 177
      do 176 j=mstart,mstop
         f(1,j,np1) = f(1,j,1)
  176 continue
      nstop = np1
  177 do 179 j=mstart,mstop
         do 178 k=nstart,nstop
            f(lp1,j,k) = f(1,j,k)
  178    continue
  179 continue
  180 continue
      if (mp .ne. 1) go to 185
      if (np .ne. 1) go to 182
      do 181 i=lstart,lstop
         f(i,1,np1) = f(i,1,1)
  181 continue
      nstop = np1
  182 do 184 i=lstart,lstop
         do 183 k=nstart,nstop
            f(i,mp1,k) = f(i,1,k)
  183    continue
  184 continue
  185 continue
      if (np .ne. 1) go to 188
      do 187 i=lstart,lstop
         do 186 j=mstart,mstop
            f(i,j,np1) = f(i,j,1)
  186    continue
  187 continue
  188 continue
      return
      end
