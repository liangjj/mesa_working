*deck @(#)prespl.f	1.2  10/27/94
      subroutine prespl(n,xdata,k,sc)
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
c       k      - order(degree plus 1)of spline approximation,
c                2.le.k.10.
c
c  output:
c
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
      dimension xdata(n), sc((4*n+1)*k)
      common /io/ inp,iout
      save
c
c  check for inconsistent or non-distinct data.
c
      if(n-k+1.le.0) then
         call lnkerr('inconsistent point/order data in spline')
      endif
      do 1 i=1,n-1
         if(xdata(i+1).eq.xdata(i)) then
            call lnkerr('inconsistent data in x variable for spline')
         endif
    1 continue       

c
c  define constants, allocate storage in scratch array and initialize
c  scratch array.
c
      limit=(4*n+1)*k
      call rzero(sc,limit)
c
c  construct knot sequence for b-splines.
c
      do 4 i=1,k
         sc(i)=xdata(1)
    4 continue      
      kdel=(k-1)/2
      kdp2=kdel+2
      kd1=k-kdel-1
      kstop=n-kdel-1
      if ((k/2)*2.d0.ge.k) then
          do 5 i=kdp2,kstop
             sc(kd1+i)=xdata(i)
    5     continue      
      else
          kstop=kstop+1
          do 7 i=kdp2,kstop
            sc(kd1+i)=.5d0*(xdata(i)+xdata(i-1))
    7     continue
      endif     
      do 9 i=1,k
         sc(n+i)=xdata(n)
    9 continue      
      return
      end
