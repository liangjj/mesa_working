c---------------------------------------------------------------------c
c           bspline routines no longer in cftmath                     c
c---------------------------------------------------------------------c
*deck bsplin
      subroutine bsplin(n,xdata,ydata,k,nbreak,break,c,iflag,sc)
c
c  bsplin uses b-splines to calculate the interpolating spline approxi-
c  mation to the data (xdata(i),ydata(i)),i=1,n.
c  reference:  carl de boor,"a practicle guide to splines",
c  (springer-verlag, ny, 1978)
c
c
c  input:
c
c       n      - number of data points.
c       xdata  - array of floating point values of independent
c                variable.  the components of xdata must be
c                distinct.
c       ydata  - array of floating point values of dependent variable.
c       k      - order(degree plus 1)of spline approximation,
c                2.le.k.10.
c
c  output:
c
c       nbreak - number of polynomial pieces in spline pp(piecewise
c                polynomial) representation.
c       break  - array of nbreak+1 breakpoints for pp-representation.
c       c      - matrix of right derivatives at the first nbreak
c                breakpoints for pp-representation.  c has k rows
c                and nbreak columns.
c       iflag
c           =1 problem has unique solution.
c           =2 insufficient number of data points to properly
c              construct fit or singular matrix is encountered.  no
c              solution is returned.
c       sc     - working storage array.  it  must have dimension
c                in the calling program of at least (4n+1)*k
c
c
c  to evaluate the fit from its pp-representation, use the routine
c  ppvalu.
c
c  the spline can also be characterized by its b-spline representation
c  and evaluated in this form using the function bvalue.  the
c  necessary knot array is contained in sc(1) - sc(n+k) and the
c  coefficient array in sc(n+k+1) - sc(n+k+n).
c
c  if the fit is needed at many points, the pp-representation
c  should be used.
c
c  for a discussion on how the b-spline knots are chosen from the
c  input data points,  see the writeup.
c
c
      implicit real*8(a-h,o-z)
      dimension xdata(n),ydata(n),break(1),c(1),sc(1)
      save
c
c  check for inconsistent or non-distinct data.
c
      if(n-k+1.le.0)go to 99
      nm1=n-1
      do 1 i=1,nm1
    1     if(xdata(i+1).eq.xdata(i))go to 99
c
c  define constants, allocate storage in scratch array and initialize
c  scratch array.
c
      np2mk=n+2-k
      km1=k-1
      ia=n+k
      iq=ia+n
      idummy=iq+n*(3*k-2)
      iscr=idummy+k
      limit=(4*n+1)*k
      do 2 i=1,limit
    2     sc(i)=0.d0
c
c  construct knot sequence for b-splines.
c
      do 4 i=1,k
    4     sc(i)=xdata(1)
      kdel=km1/2
      kdp2=kdel+2
      kd1=k-kdel-1
      kstop=n-kdel-1
      if((k/2)*2.lt.k)go to 6
      do 5 i=kdp2,kstop
    5     sc(kd1+i)=xdata(i)
      go to 8
    6 kstop=kstop+1
      do 7 i=kdp2,kstop
    7     sc(kd1+i)=.5d0*(xdata(i)+xdata(i-1))
    8 do 9 i=1,k
    9     sc(n+i)=xdata(n)
c
c  set up matrix to solve for coefficients in b-spline representation.
c
      do 30 i=1,n
          call interv(sc(k),np2mk,xdata(i),ileft,mflag)
          ileft=ileft+km1
          if(mflag)99,15,14
   14     if(i .lt. n)go to 99
          ileft=n
   15     call bsplvb(sc,k,1,xdata(i),ileft,sc(idummy+1))
          l=ileft-i
          do 16 j=1,k
              l=l+1
   16         sc(iq+i+n*(l-1))=sc(idummy+j)
          if(sc(iq+i+n*(k-1)).eq.0.d0)go to 99
   30     sc(ia+i)=ydata(i)
c
c  solve for coefficients in b-spline representation.
c
      call banmat(n,km1,km1,1,1,sc(iq+1),n,sc(ia+1),n,det,sc(idummy+1))
c
c  convert b-spline representation to ppiecewise polynomial
c  representation.
c
      call bsplpp(sc,sc(ia+1),n,k,sc(idummy+1),break,c,nbreak)
      iflag=1
      return
   99 iflag=2
      return
      end
