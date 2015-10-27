*deck %W%  %G%
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
*deck %W%  %G%
      subroutine splmat(n,xdata,k,sc)
      implicit real*8(a-h,o-z)
      dimension sc((4*n+1)*k), xdata(n)
      common /io/ inp, iout
c
c  set up matrix to solve for coefficients in b-spline representation.
c
      ia=n+k
      iq=ia+n
      idummy=iq+n*(3*k-2)
      iscr=idummy+k
      do 30 i=1,n
          call interv(sc(k),n+2-k,xdata(i),ileft,mflag)
          ileft=ileft+k-1
          if (mflag.lt.0) then
              call lnkerr('error in splmat')
          else
              if (mflag.gt.0) then
                  if(i .lt. n) then
                     call lnkerr('error in splmat')
                  endif
                  ileft=n
              endif
          endif
          call bsplvb(sc,k,1,xdata(i),ileft,sc(idummy+1))
          l=ileft-i
          do 16 j=1,k
              l=l+1
              sc(iq+i+n*(l-1))=sc(idummy+j)
   16     continue      
          if(sc(iq+i+n*(k-1)).eq.0.d0) then
             call lnkerr('error in splmat')
          endif
   30 continue      
      return
      end      
*deck %W%  %G%
      subroutine interv(xt,lxt,x,ileft,mflag)
c***begin prologue  interv
c***refer to  fc
c
c computes largest ileft in (1,lxt) such that xt(ileft) .le. x
c***routines called  (none)
c***end prologue  interv
      implicit real *8 (a-h,o-z)
      dimension xt(lxt)
      data ilo /1/
c***first executable statement  interv
      ihi = ilo + 1
      if (ihi .lt. lxt)                go to 20
         if (x .ge. xt(lxt))           go to 110
         if (lxt .le. 1)               go to 90
         ilo = lxt - 1
                                       go to 21
   20 if (x .ge. xt(ihi))              go to 40
   21 if (x .ge. xt(ilo))              go to 100
c *** now x .lt. xt(ihi) . find lower bound
   30 istep = 1
   31 ihi = ilo
      ilo = ihi - istep
      if (ilo .le. 1)                  go to 35
      if (x .ge. xt(ilo))              go to 50
      istep = istep*2
                                       go to 31
   35 ilo = 1
      if (x .lt. xt(1))                go to 90
                                       go to 50
c *** now x .ge. xt(ilo) . find upper bound
   40 istep = 1
   41 ilo = ihi
      ihi = ilo + istep
      if (ihi .ge. lxt)                go to 45
      if (x .lt. xt(ihi))              go to 50
      istep = istep*2
                                       go to 41
   45 if (x .ge. xt(lxt))              go to 110
      ihi = lxt
c *** now xt(ilo) .le. x .lt. xt(ihi) . narrow the interval
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)             go to 100
c     note. it is assumed that middle = ilo in case ihi = ilo+1
      if (x .lt. xt(middle))           go to 53
         ilo = middle
                                       go to 50
   53    ihi = middle
                                       go to 50
c *** set output and return
   90 mflag = -1
      ileft = 1
                                       return
  100 mflag = 0
      ileft = ilo
                                       return
  110 mflag = 1
      ileft = lxt
                                       return
      end
*deck %W%  %G%
      subroutine bsplvb ( t, jhigh, index, x, ileft, biatx )
      implicit real *8 (a-h,o-z)
calculates the value of all possibly nonzero b-splines at  x  of
c  order max(jhigh,(j+1)(index-1)) on  t .
c     dimension (t(ileft+jhigh)
      save
      dimension t(1),biatx(jhigh), deltam(20),deltap(20)
      data j/1/
content of j, deltam, deltap is expected unchanged between calls.
                                        go to (10,20), index
   10 j = 1
      biatx(1) = 1.
      if (j .ge. jhigh)                 go to 99
c
   20    jp1 = j + 1
         deltap(j) = t(ileft+j) - x
         deltam(j) = x - t(ileft-j+1)
         saved = 0.d0
         do 26 i=1,j
            term = biatx(i)/(deltap(i) + deltam(jp1-i))
            biatx(i) = saved + deltap(i)*term
   26       saved = deltam(jp1-i)*term
         biatx(jp1) = saved
         j = jp1
         if (j .lt. jhigh)              go to 20
c
   99                                   return
      end
*deck %W%  %G%
      subroutine bsplpp ( t, bcoef, n, k, scrtch, break, coef, l )
converts b-representation to pp-representation.
c     dimension t(n+k), break(*+1), coef(k,*)
c  here, * = the final value of the output parameter l .
      implicit real *8 (a-h,o-z)
      save
      dimension t(1),bcoef(n), scrtch(k,k), break(1),coef(k,1)
      dimension biatx(20)
      l = 0
      break(1) = t(k)
      do 50 ileft=k,n
c        find the next nontrivial knot interval.
         if (t(ileft+1) .eq. t(ileft))  go to 50
         l = l + 1
         break(l+1) = t(ileft+1)
         if (k .gt. 1)                  go to 9
         coef(1,l) = bcoef(ileft)
                                        go to 50
c        store the k b-spline coeff.s relevant to current knot interval
c        in  scrtch(.,1) .
    9    do 10 i=1,k
   10       scrtch(i,1) = bcoef(ileft-k+i)
c        for j=1,...,k-1, compute the k-j b-spline coeff.s relevant to
c        current knot interval for the j-th derivative by differencing
c        those for the (j-1)st derivative, and store in scrtch(.,j+1) .
         do 20 jp1=2,k
            j = jp1 - 1
            kmj = k - j
            fkmj = float(kmj)
            do 20 i=1,kmj
               diff = t(ileft+i) - t(ileft+i - kmj)
               if (diff .gt. 0.d0)  scrtch(i,jp1) =
     *                       ((scrtch(i+1,j)-scrtch(i,j))/diff)*fkmj
   20          continue
c        starting with the one b-spline of order 1 not zero at t(ileft),
c        find the values at t(ileft) of the j+1 b-splines of order j+1
c        not identically zero there from those of order j, then combine
c        with the b-spline coeff.s found earlier to compute the (k-j)-
c        th derivative at t(ileft) of the given spline.
         call bsplvb ( t, 1, 1, t(ileft), ileft, biatx )
         coef(k,l) = scrtch(1,k)
         do 30 jp1=2,k
            call bsplvb ( t, jp1, 2, t(ileft), ileft, biatx )
            kmj = k+1 - jp1
            sum = 0.d0
            do 28 i=1,jp1
   28          sum = biatx(i)*scrtch(i,kmj) + sum
   30       coef(kmj,l) = sum
   50    continue
                                        return
      end
*deck %W%  %G%
      subroutine ppval (x,y,break,coef,n,nbreak,order,ind,jderiv)
calculates value at  x  of  jderiv-th derivative of spline from pp-repr
      implicit real *8 (a-h,o-z)
      integer order
      common /io/ inp, iout
      save
      dimension break(nbreak), coef(order,nbreak), x(n), y(n), ind(n)
      call rzero(y,n)
      fmmjdr = order - jderiv
c  derivatives of order  k  or higher are identically zero.
      if (fmmjdr .gt. 0.d0) then
c  index of largest breakpoint to the left of  x is in ind.
c  evaluate  jderiv -th derivative of  i -th polynomial piece at  x .
          lim=order-jderiv
          do 10 mm=1,lim
             m=order-mm+1
             do 20 i=1,n
                h = x(i) - break(ind(i))
                y(i) = (y(i)/fmmjdr)*h + coef(m,ind(i))
   20        continue
             fmmjdr = fmmjdr - 1.d0
   10     continue     
      endif
      return
      end
*deck %W%  %G%
      subroutine banmat (in,l1,l2,nt,im,a,ia,y,iy,det,int)
      implicit real *8 (a-h,o-z)
      save
      dimension a(ia,ia), y(iy,iy), int(in)
      n=in
      l=im
      ml=l1
      mu=l2
      ll=ml+mu+1
      if (n.le.1.or.ia.lt.n.or.l1.le.0.or.l2.le.0) go to 190
      lu=ll+1
      lp=ll+ml
      n1=n-1
      if (nt.ne.1) go to 100
      sn=1.d0
      do 30 i=1,ml
      ii=mu+i
      k=ml+1-i
      do 10 j=1,ii
   10 a(i,j)=a(i,j+k)
      k=ii+1
      do 20 j=k,ll
   20 a(i,j)=0.d0
   30 continue
      do 90 nr=1,n1
      np=nr+1
      lr=nr+ml
      if (lr.gt.n) lr=n
      lc=nr+ll-1
      if (lc.gt.n) lc=n
      kk=lc-nr
      mx=nr
      xm=abs(a(nr,1))
      do 40 i=np,lr
      if (abs(a(i,1)).le.xm) go to 40
      mx=i
      xm=abs(a(i,1))
   40 continue
      int(nr)=mx
      if (mx.eq.nr) go to 60
      do 50 i=1,ll
      xx=a(nr,i)
      a(nr,i)=a(mx,i)
   50 a(mx,i)=xx
      sn=-sn
   60 xm=a(nr,1)
      if (xm.eq.0.d0) go to 190
      xm=-1.d0/xm
      do 80 i=np,lr
      j=ll+i-nr
      xx=a(i,1)*xm
      a(nr,j)=xx
      do 70 kt=1,kk
      kt1=(kt-1)*ia
      kti=kt1+i
      ktr=kt1+nr
      a(kti,1)=xx*a(ktr,2)+a(kti,2)
   70 continue
   80 a(i,ll)=0.d0
   90 continue
      if (a(n,1).eq.0.d0) sn=0.d0
  100 det=sn
      if (im.eq.0) return
      do 150 nr=1,n1
      np=nr+1
      lr=nr+ml
      if (lr.gt.n) lr=n
      kk=lr-nr
      if (a(nr,1).eq.0.d0) go to 190
      if (int(nr).eq.nr) go to 120
      j=int(nr)
      do 110 i=1,l
      xx=y(nr,i)
      y(nr,i)=y(j,i)
  110 y(j,i)=xx
  120 do 140 i=1,l
      do 130 kt=1,kk
      ktp=(kt-1)+np
      ktr=(kt-1)*ia+nr
      y(ktp,i)=y(nr,i)*a(ktr,lu)+y(ktp,i)
  130 continue
  140 continue
  150 continue
      xm=a(n,1)
      nr=n
      if (xm.eq.0.d0) go to 190
      xm=1.d0/xm
      do 160 i=1,l
  160 y(n,i)=y(n,i)*xm
      ns=1
      do 180 nb=ns,n1
      nr=n-nb
      np=nr+1
      lc=nr+ll-1
      if (lc.gt.n) lc=n
      kk=lc-nr
      do 170 i=1,l
  170 y(nr,i)=(y(nr,i)-sdot(kk,a(nr,2),ia,y(np,i),1))/a(nr,1)
  180 continue
      return
  190 det=0.d00
      return
      end
*deck %W%  %G%
      subroutine splcof(n,ydata,k,nbreak,break,c,sc)
      implicit real*8(a-h,o-z)
      dimension break(nbreak+1), ydata(n), c(k,nbreak), sc((4*n+1)*k)   
      common /io/ inp, iout       
c
c  solve for coefficients in b-spline representation.
c
      ia=n+k
      call copy(ydata,sc(ia+1),n)
      iq=ia+n
      idummy=iq+n*(3*k-2)
      call banmat(n,k-1,k-1,1,1,sc(iq+1),n,sc(ia+1),n,det,sc(idummy+1))
c
c  convert b-spline representation to ppiecewise polynomial
c  representation.
c
      call bsplpp(sc,sc(ia+1),n,k,sc(idummy+1),break,c,nbreak)
      return
      end
