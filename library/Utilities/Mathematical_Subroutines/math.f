*deck banmat
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
*deck bffgh
      subroutine bffgh (z,l,ba,bjp,bjpp,bb,byp,bypp,qsq,ikq,y)
c     c3 ucsd bffgh
      implicit real*8 (a-h,o-z)
      dimension bj(5), by(5), r(3)
c
c     revised 1-22-68
c
      if (l.gt.0) go to 20
      if (ikq.ge.2.and.qsq.lt.0.d0) go to 10
      ba=sin(z)/z
      bb=-cos(z)/z
      uta=ba/z
      utb=bb/z
      bjp=-bb-uta
      byp=ba-utb
      ujtemp=(2.0d0/z**2)-1.0d0
      bjpp=ba*ujtemp+2.0d0*utb
      bypp=bb*ujtemp-2.0d0*uta
      return
   10 exx=exp(z)
      exxi=1.d0/exx
      ba=(exx-exxi)*0.5d0
      ba=ba/z
      bb=exxi/z
      bbb=(exx+exxi)*0.5d0
      bbb=bbb/z
      bjp=bbb-ba/z
      byp=-bb-bb/z
      bjpp=ba-2.d0*bjp/z
      bypp=bb-2.d0*byp/z
      return
   20 ysave=y
      if (y-1.d0) 40,40,30
   30 y=1.
   40 fl=l
      l1=l+1
      zb=z
      ze=0.d0
      fl1=fl+1.d0
      zsq=z*z
      flr=sqrt(fl*fl1)
      sine=sin(z)
      cose=cos(z)
      if (ikq-1) 60,60,70
   50 if (z-flr) 410,410,570
   60 if (z-1.d0) 200,90,90
   70 if (qsq) 80,80,60
   80 if (z-1.d0) 330,320,320
   90 bj(1)=sine/zb
      sign=1.d0
      er=-y
      by(1)=er*cose
  100 by(1)=by(1)/zb
      er=bj(1)/zb
      if (sign) 110,120,120
  110 er=-er
  120 bj(2)=(by(1)/y+er)/y
      by(2)=(by(1)/z-y*bj(1))*y
      flx=2.d0
  130 temp=(2.d0*flx-1.d0)/z
      bj(3)=(temp*bj(2)-bj(1)/y)/y
      if (sign) 140,150,150
  140 bj(3)=-bj(3)
  150 er=temp*by(2)
      by(3)=y*by(1)
      if (sign) 160,170,170
  160 by(3)=-by(3)
  170 by(3)=((er-by(3)))*y
      if (flx-fl1) 190,180,180
  180 if (sign) 340,340,50
  190 bj(1)=bj(2)
      by(1)=by(2)
      bj(2)=bj(3)
      by(2)=by(3)
      flx=flx+1.d0
      go to 130
  200 lx=l-1
      sexit=1.d0
  210 do 310 i=1,3
      sum2=1.d0
      sum1=1.d0
      term1=1.d0
      term2=1.d0
      fm=1.d0
      tl=2*lx
  220 tm=2.d0*fm
      term1=-term1*zsq/(tm*(tl+tm+1.d0))*sexit
      term2=term2*zsq/(tm*(tl-tm+1.d0))*sexit
      old1=sum1
      old2=sum2
      sum1=sum1+term1
      sum2=sum2+term2
      diff=abs((old1-sum1)/sum1)+abs((old2-sum2)/sum2)
      if (diff-.000001d0) 240,240,230
  230 fm=fm+1.d0
      go to 220
  240 ix=2*lx-1
      asm=-1.d0
      if (lx) 270,270,250
  250 do 260 j=1,ix,2
      aj=j
  260 asm=asm*aj
  270 zyl=(z/y)**lx
      zyl1=zyl*(z/y)
      if (lx) 280,280,290
  280 asm1=1.d0
      go to 300
  290 alx=lx
      asm1=-asm*(2.d0*alx+1.d0)
  300 bj(i)=zyl*sum1/asm1
      by(i)=asm*sum2/zyl1
  310 lx=lx+1
      if (sexit) 340,340,570
  320 sign=-1.d0
      er=exp(zb)
      erl=1.d0/er
      snh=.5d0*(er-erl)
      csh=.5d0*(er+erl)
      bj(1)=snh/zb
      by(1)=y*csh
      go to 100
  330 lx=l-1
      sexit=-1.d0
      go to 210
  340 knt=2
      by(1)=y/zb
      by(2)=y*y*(1.d0+zb)/(zb*zb)
  350 by(3)=(2.d0*float(knt)-1.d0)*y*by(2)/zb+y*y*by(1)
      if (knt-l-1) 360,370,370
  360 by(1)=by(2)
      by(2)=by(3)
      knt=knt+1
      go to 350
  370 exz=exp(-zb)
      do 380 n=1,3
  380 by(n)=by(n)*exz
      if (zb-1.d0) 400,390,390
  390 if (z-4.d0*flr) 500,400,400
  400 er=2.d0*fl+1.d0
      bjp=(fl*bj(1)/y+fl1*bj(3)*y)/er
      byp=-(fl*by(1)*y+fl1*by(3)/y)/er
      bjpp=2.d0*bjp/z-(1.d0+fl*fl1/zsq)*bj(2)
      bypp=2.d0*byp/z-(1.d0+fl*fl1/zsq)*by(2)
      go to 580
  410 a=0.1d0
      ize=0
      b=0.35
      l1=l+1
      fn1=l
      idx=-1
      fn2=z-0.5d0+sqrt(30.0d0*b*z)
      ze=0.d0
      u1=2.d0*z/(2.d0*fn1+1.d0)
      u2=2.d0*z/(2.d0*fn2+1.d0)
      m1=fn1+30.d0*(a+b*u1*(2.d0-u1*u1)/(2.d0*(1.d0-u1*u1)))
      m2=fn2+30.d0*(a+b*u2*(2.d0-u2*u2)/(2.d0*(1.d0-u2*u2)))
      al=l
      if (fn2-al) 420,420,430
  420 am=m1+1
      m=m1
      go to 450
  430 if (m1-m2) 420,420,440
  440 am=m2+1
      m=m2
  450 r(3)=z/(2.d0*am+3.d0)
      n=m
  460 an=n
      r(2)=z/(2.d0*an+3.d0-z*r(3))
      n=n-1
      if (l-n) 470,470,480
  470 r(3)=r(2)
      go to 460
  480 bj(3)=r(2)
      bj(2)=1.d0
      al=2*l+1
      bj(1)=al/z-r(2)
      la=-l
      la1=-l-1
      alpha=z*z*(bj(2)*by(1)*y**la-bj(1)*by(2)*y**la1)
      na=l-1
      do 490 n=1,3
      bj(n)=(1.d0/(y**na*alpha))*bj(n)
  490 na=na+1
      go to 570
  500 a=0.1
      b=.35
      fn=-.5d0+sqrt(30.d0*b*zb)
      if (fn-fl) 510,520,520
  510 fn=fl
  520 u=2.d0*zb/(2.d0*fn+1.d0)
      m=fn+30.d0*(a+b*u)+1.d0
      r(3)=zb/2.d0*float(m)+1.d0
      m=m-1
  530 an=m
      r(2)=zb/(2.d0*an+1.d0+zb*r(3))
      m=m-1
      if (m-l-1) 550,540,540
  540 r(3)=r(2)
      go to 530
  550 bj(3)=r(2)
      bj(2)=1.d0
      bj(1)=(2.d0*fl+1.d0)/zb+r(2)
      alpha=zb*zb*(bj(2)*by(1)/y**l+bj(1)*by(2)/y**l1)
      na=fl-1.d0
      do 560 n=1,3
      bj(n)=(1.d0/(y**na*alpha))*bj(n)
  560 na=na+1
      go to 400
  570 er=1.d0/(2.d0*fl+1.d0)
      bjp=er*(fl*bj(1)/y-fl1*bj(3)*y)
      byp=er*(fl*by(1)*y-fl1*by(3)/y)
      bjpp=(fl*fl1/zsq-1.d0)*bj(2)-2.d0*bjp/z
      bypp=(fl*fl1/zsq-1.d0)*by(2)-2.d0*byp/z
  580 ba=bj(2)
      bb=by(2)
      y=ysave
      return
      end
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
      save
      dimension xdata(n),ydata(n),break(1),c(1),sc(1)
c
c  check for inconsistent or non-distinct data.
c
      implicit real *8 (a-h,o-z)
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
*deck bsplpp
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
*deck bsplvb
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
*deck fact
c***begin prologue     fact
c***date written       880721   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           fact, link 2702, factorials
c***author             schneider, barry (lanl)
c***source             m2702
c***purpose            factorials
c***description        calculation of factorials
c***references         none
c
c***routines called
c***end prologue       fact
      subroutine fact(dfct,ddfct,maxfac)
      implicit integer (a-z)
      real *8 dfct, ddfct
      dimension dfct(0:maxfac), ddfct(0:maxfac)
c----------------------------------------------------------------------c
c               calculate factorials                                   c
c----------------------------------------------------------------------c
      dfct(0)=1.d+00
      dfct(1)=1.d+00
      if (maxfac.gt.1) then
          do 10 i=2,maxfac
             dfct(i)=i*dfct(i-1)
   10     continue
      endif
c----------------------------------------------------------------------c
c           calculate (2*m-1) double factorial                         c
c----------------------------------------------------------------------c
      ddfct(0)=1.d+00
      ddfct(1)=1.d+00
      ddfct(2)=3.d+00
      if (maxfac.gt.2) then
         do 20 i=3,maxfac
            ddfct(i)=(i+i-1)*ddfct(i-1)
   20    continue  
      endif
      return
      end
*deck legend
c***begin prologue     legend
c***date written       880721   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           legend, link 2702, legendre functions
c***author             schneider, barry (lanl)
c***source             m2702
c***purpose            legendre functions
c***description        calculation of p(l,m) functions
c***references         none
c
c***routines called
c***end prologue       legend
      subroutine legend (plm,x,dfct,ddfct,npt,lmax,m,maxfac)
      implicit integer (a-z)
      real *8 plm, x, dfct, ddfct, fm, facx, f1
      dimension plm(npt,0:lmax), x(npt)
      dimension  dfct(0:maxfac), ddfct(0:maxfac)
c----------------------------------------------------------------------c
c           start recursion with plm(m,m) and plm(m+1,m)               c
c                      and recur upward                                c
c----------------------------------------------------------------------c
      do 10 i=m,lmax
         do 10 j=1,npt
            plm(j,i)=0.d+00
   10 continue
      fm=.5d+00*m
      do 20 i=1,npt
         facx=1.d+00
         if (fm.ne.0.d+00) then
             facx=(1.d+00-x(i)*x(i))**fm
         endif
         plm(i,m)=ddfct(m)*facx
   20 continue
      if (lmax.ne.m) then
          mm=m+m+1
          mpls1=m+1
          do 30 i=1,npt
             plm(i,mpls1)=mm*x(i)*plm(i,m)
   30     continue
          if (lmax.ne.mpls1) then
              lind=m+2
              n1=2
              n2=m+m+3
              n3=n2-2
              do 50 i=lind,lmax
                 ii=i-1
                 jj=i-2
                 do 40 j=1,npt
                    f1=n2*x(j)*plm(j,ii)-n3*plm(j,jj)
                    f1=f1/n1
                    plm(j,i)=f1
   40            continue
                 n1=n1+1
                 n2=n2+2
                 n3=n3+1
   50       continue
          endif
      endif
      return
c
      end
*deck lgndx2
c***begin prologue     lgndx2
c***date written       861117   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           gauss-legendre, quadrature
c***author             collins lee (lanl)
c***source             mylib
c***purpose            quadrature
c***description        points and weights for gauss-legendre quadrature
c***                   for 1, 2, 3 points
c***references         none
c
c***routines called    none
      subroutine lgndx2 (n,j,wtg,rtg)
      implicit real *8 (a-h,o-z)
*
      data x11, w11 / .5d0 , 1.d0/
      data x21, w21 /0.5773502692d0,1.d0/
      data x31, x32 /0.7745966692d0,0.0d0/
      data w31, w32 /0.55555556d0,0.88888888d0/
      data half, one /0.5d0,1.0d0/
*
      if (n.eq.1) then
          rtg=x11
          wtg=w11
      else
      if (n.eq.3) go to 20
*
      wtg=w21
      if (j.eq.2) go to 10
      rtg=-x21
      go to 60
   10 rtg=x21
      go to 60
   20 go to (30,40,50), j
   30 wtg=w31
      rtg=-x31
      go to 60
   40 wtg=w32
      rtg=x32
      go to 60
   50 wtg=w31
      rtg=x31
   60 wtg=wtg*half
      rtg=half*(rtg+one)
      endif
      return
      end
*deck modbes
c***begin prologue     modbes
c***date written       890801   (yymmdd)
c***revision date               (yymmdd)
c***keywords           m6004, link 6004, modified bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            modified spherical bessel funcions
c***description        the analog of ricatti-bessel functions, the
c***                   modified spherical bessel functions multiplied
c***                   by x are computed by upward recursion. this
c***                   routine is for the exponentially decaying form
c***                   needed for closed channel asymptotic solutions.
c***references       
c
c***routines called
c***end prologue       modbes
      subroutine modbes (x,xinv,mn,dmn,ptdim,lmax,prnt)
      implicit integer (a-z)
      real *8 x, xinv, mn, dmn, ex, pi2
      logical prnt
      dimension x(ptdim), xinv(ptdim), mn(ptdim,0:lmax)
      dimension dmn(ptdim,0:lmax)
      common /io/ inp, iout
      data pi2 / 1.570796327d+00 /
c----------------------------------------------------------------------c
c       calculate first two members to start upward recursion          c
c----------------------------------------------------------------------c
      do 10 i=1,ptdim
         ex=exp(-x(i))
         mn(i,0)=pi2*ex
         dmn(i,0)=-mn(i,0)
   10 continue
      if (lmax.ge.1) then
          do 20 i=1,ptdim
             ex=exp(-x(i))
             mn(i,1)=pi2*(1.d+00+xinv(i))*ex 
             dmn(i,1)=-pi2*(1.d+00+xinv(i)+xinv(i)*xinv(i))*ex
   20     continue
      endif
      if (lmax.gt.1) then
          do 30 l=1,lmax-1
             lp=l+1
             lm=l-1
             lfac=l+l+1
             do 40 i=1,ptdim
                mn(i,lp)=mn(i,lm)+lfac*xinv(i)*mn(i,l)
                dmn(i,lp)=dmn(i,lm)+lfac*xinv(i)*dmn(i,l)
   40        continue
   30     continue
      endif
      if (prnt) then
          do 50 l=0,lmax
             write (iout,900) l
             write (iout,1000) (mn(i,l),i=1,ptdim)
             write (iout,1010) l
             write (iout,1000) (dmn(i,l),i=1,ptdim)
   50     continue
      endif
  900 format(/,5x,'regular modified function l = ',1x,i3)
 1000 format( (/,5x,5(e15.8,1x) ) )
 1010 format(/,5x,'derivative regular modified function l = ',1x,i3)
      return
      end
c
c                       least squares polynomial fitting
*deck plyfit
c***begin prologue     plyfit
c***date written       850601  (yymmdd)
c***revision date      yymmdd  (yymmdd)
c***keywords           least squares,
c***author             schneider, barry (lanl)
c***source             mylib
c***purpose
c                      least squares polynomial fit of input function
c***description
c                      input function f is fit in a least squares sense
c                      to a power series
c
c                      f is input function
c                      power are the powers used in fitting
c                      pt are the points
c                      rhs is the output set of coefficients
c                      the rest of the variables are used only internally
c***references
c***routines called    sgemm(clams)
c***end prologue       plyfit
      subroutine plyfit(f,coef,rhs,xn,pt,power,ipvt,npwr,npnts,prnt,
     1                  first)
      implicit integer (a-z)
      real *8 f, coef, rhs, xn, pt
      logical prnt
      character *80 title
      character *(*) first
      dimension f(npnts), power(npwr), xn(npwr,npnts), pt(npnts)
      dimension coef(npwr,npwr), rhs(npwr), ipvt(npwr)
      common /io/ inp, iout
      if (first.eq.'first') then
c----------------------------------------------------------------------c
c                     print out input information                      c
c----------------------------------------------------------------------c
      write (iout,100) npwr
  100 format(/,5x,'least squares polynomial fitting subroutine',//,5x,
     1            'order of fit',1x,i3)
      write (iout,200) (power(i),i=1,npwr)
  200 format (//,5x,'powers',(/,15x,10(i2,1x)))
      call rzero(xn,npwr*npnts)
c----------------------------------------------------------------------c
c                    set up coefficient matrix                         c
c----------------------------------------------------------------------c
      do 10 i=1,npwr
         do 20 j=1,npnts
            xn(i,j)=pt(j)**power(i)
   20    continue
   10 continue      
      call rzero(coef,npwr*npwr)
      call ebct(coef,xn,xn,npwr,npnts,npwr)
c----------------------------------------------------------------------c
c               factor coefficient matrix                              c
c----------------------------------------------------------------------c
      call sgefa(coef,npwr,npwr,ipvt,info)
      else
c----------------------------------------------------------------------c
c                get right hand side                                   c
c----------------------------------------------------------------------c
      call ebc(rhs,xn,f,npwr,npnts,1)
c----------------------------------------------------------------------c
c               solve equations for coefficients                       c
c----------------------------------------------------------------------c
      call rzero(rhs,npwr)
      call sgesl(coef,npwr,npwr,ipvt,rhs,0)
      if (prnt) then
         title='least squares coefficients'
         call prntrm(title,rhs,npwr,1,npwr,1,iout)
      endif
      endif
      return
      end
*deck ppvalu
      function ppvalu (break, coef, l, k, x, jderiv )
calculates value at  x  of  jderiv-th derivative of spline from pp-repr
      implicit real *8 (a-h,o-z)
      save
      dimension break(l), coef(k,l)
      ppvalu = 0.d0
      fmmjdr = k - jderiv
c  derivatives of order  k  or higher are identically zero.
      if (fmmjdr .le. 0.d0)               go to 99
c  find index i of largest breakpoint to the left of  x .
      call interv ( break, l, x, i, ndummy )
c  evaluate  jderiv -th derivative of  i -th polynomial piece at  x .
      h = x - break(i)
      lim=k-jderiv
      do 10 mm=1,lim
         m=k-mm+1
         ppvalu = (ppvalu/fmmjdr)*h + coef(m,i)
   10    fmmjdr = fmmjdr - 1.d0
   99                                   return
      end
c----------------------------------------------------------------------c
c                   ricatti-bessel or bessel function                  c
c                           programs                                   c
c----------------------------------------------------------------------c
*deck rbes
      subroutine rbes (type,l,z,b,bp,bpp,y,yp,ypp)
      implicit real*8 (a-h,o-z)
      character *(*) type
      call bffgh (z,l,b,bp,bpp,y,yp,ypp,z*z,2,1.d0)
c----------------------------------------------------------------------c
c                make ricatti-bessels if you want                      c
c----------------------------------------------------------------------c   
      if (type.eq.'ricatti-bessel') then
          bpp=z*bpp+2.d+00*bp
          ypp=z*ypp+2.d+00*yp
          bp=z*bp+b
          yp=z*yp+y
          b=z*b
          y=z*y
      endif
      return
      end
*deck rcbesb
c***begin prologue     rcbesb
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate bessel functions
c***description        bessel functions calculated backward
c***                   recursion. if greater accuracy needed change
c***                   parameter statement.
c***references       
c
c***routines called
c***end prologue       rcbesb
      subroutine rcbesb(x,xinv,j,jp,y,yp,norm,np,ptdim,lmax,ltop,prnt)
      implicit integer (a-z)
      parameter (acc=30)
      common /io/ inp, iout
      real *8 one, two, x, j, jp, y, yp, xinv, norm, rl
      logical prnt
      dimension x(ptdim), xinv(ptdim), j(ptdim,0:ltop), jp(ptdim,0:ltop)
      dimension y(ptdim,0:ltop), yp(ptdim,0:ltop), norm(ptdim)
      data one, two/ 1.d+00, 2.d+00 /
c----------------------------------------------------------------------c
c            estimate starting l                                       c
c----------------------------------------------------------------------c
      rl=lmax*acc
      strtl=lmax+sqrt(rl)
      strtl=max(strtl,lfinal)
      if (strtl.gt.ltop) then
          call lnkerr('starting l bigger than ltop')
      endif
c----------------------------------------------------------------------c
c               make first two bessel functions                        c
c----------------------------------------------------------------------c
      do 10 i=1,np
         j(i,0)=sin(x(i))
         norm(i)=j(i,0)
         y(i,0)=-cos(x(i))
         jp(i,0)=-y(i,0)
         yp(i,0)=j(i,0)
   10 continue
      if (lmax.eq.0) then
          return
      endif
      do 20 i=1,np
         j(i,1)=j(i,0)*xinv(i)+y(i,0)
         y(i,1)=y(i,0)*xinv(i)-j(i,0)
         jp(i,1)=j(i,0)-one*j(i,1)*xinv(i)
         yp(i,1)=y(i,0)-one*y(i,1)*xinv(i)    
   20 continue
      if (lmax.eq.1) then
          return
      endif       
c----------------------------------------------------------------------c
c                   lmax is greater than one                           c
c----------------------------------------------------------------------c
c----------------------------------------------------------------------c
c              calculate y by upward recursion                         c
c----------------------------------------------------------------------c
       lfinal=min(ltop,lmax+1)     
       do 30 l=1,lfinal-1
          ll1=l+l+1
          lm1=l-1
          lp1=l+1
          do 40 i=1,np
             y(i,lp1)=ll1*y(i,l)*xinv(i)-y(i,lm1)
   40     continue
   30  continue
c----------------------------------------------------------------------c
c             calculate j by downward recursion                        c
c----------------------------------------------------------------------c
       call rzero(j(1,strtl),np)
       onelss=strtl-1
       do 50 i=1,np
          j(i,onelss)=one
   50  continue
       do 60 l=onelss-1,0,-1
          ll3=l+l+3
          lp1=l+1
          lp2=l+2
          do 70 i=1,np
             j(i,l)=ll3*j(i,lp1)*xinv(i)-j(i,lp2)
   70     continue
   60  continue
c----------------------------------------------------------------------c
c                normalize the j                                       c
c----------------------------------------------------------------------c
       do 80 i=1,np
          norm(i)=norm(i)/j(i,0)
  80   continue
       do 90 l=0,lfinal
          do 100 i=1,np
             j(i,l)=j(i,l)*norm(i)
  100     continue
   90  continue
c---------------------------------------------------------------------c
c             finish calculation by getting derivatives               c
c---------------------------------------------------------------------c
       do 200 l=0,lmax
          lp1=l+1
          do 300 i=1,np
             jp(i,l)=lp1*j(i,l)*xinv(i)-j(i,lp1)            
             yp(i,l)=lp1*y(i,l)*xinv(i)-y(i,lp1) 
  300    continue
  200 continue
      if (prnt) then
          do 800 l=0,lmax
             write (iout,900) l
             write (iout,1000) (j(i,l),i=1,np)
             write (iout,1010) l
             write (iout,1000) (jp(i,l),i=1,np)
             write (iout,1020) l
             write (iout,1000) (y(i,l),i=1,np)
             write (iout,1030) l
             write (iout,1000) (yp(i,l),i=1,np)
  800     continue
      endif
  900 format(/,5x,'regular function l = ',1x,i3)
 1000 format( (/,5x,5(e15.8,1x) ) )
 1010 format(/,5x,'derivative regular function l = ',1x,i3)
 1020 format(/,5x,'irregular function l = ',1x,i3)
 1030 format(/,5x,'derivative irregular function l = ',1x,i3)
      return
      end
*deck rcbesf
c***begin prologue     rcbesf
c***date written       xxxxxx   (yymmdd)
c***revision date      890422   (yymmdd)
c***keywords           m6004, link 6004, bessel, spline
c***author             schneider, barry (lanl)
c***source             m6004
c***purpose            generate bessel functions
c***description        bessel functions calculated by forward
c***                   recursion. use only when lmax lt arg
c***references       
c
c***routines called
c***end prologue       rcbesf
      subroutine rcbesf(x,xinv,j,jp,y,yp,nr,ptdim,lmax,ltop,prnt) 
      implicit integer (a-z)
      common /io/ inp, iout
      logical prnt
      real *8 x, j, jp, y, yp, xinv
      dimension j(ptdim,0:lmax), jp(ptdim,0:lmax), x(*), xinv(*)
      dimension y(ptdim,0:lmax), yp(ptdim,0:lmax)
      do 10 i=1,nr
         j(i,0)=sin(x(i))
         y(i,0)=-cos(x(i))
         jp(i,0)=-y(i,0)
         yp(i,0)=j(i,0)
   10 continue
      if (lmax.eq.0) then
          return
      else
          do 20 i=1,nr
             y(i,1)=y(i,0)*xinv(i)-j(i,0)
             j(i,1)=j(i,0)*xinv(i)+y(i,0)
   20     continue
          if (lmax.ne.1) then
              do 30 l=1,lmax-1
                 lp1=l+1
                 ll1=l+l+1
                 lm1=l-1
                 do 40 i=1,nr
                    j(i,lp1)=ll1*j(i,l)*xinv(i)-j(i,lm1)
                    y(i,lp1)=ll1*y(i,l)*xinv(i)-y(i,lm1)
   40            continue
   30         continue
          endif
          do  50 l=1,lmax
              lm1=l-1
              do 60 i=1,nr
                 jp(i,l)=j(i,lm1)-l*j(i,l)*xinv(i)
                 yp(i,l)=y(i,lm1)-l*y(i,l)*xinv(i)
   60         continue
   50     continue
      endif
      if (prnt) then
          do 800 l=0,lmax
             write (iout,900) l
             write (iout,1000) (j(i,l),i=1,nr)
             write (iout,1010) l
             write (iout,1000) (jp(i,l),i=1,nr)
             write (iout,1020) l
             write (iout,1000) (y(i,l),i=1,nr)
             write (iout,1030) l
             write (iout,1000) (yp(i,l),i=1,nr)
  800     continue
      endif
  900 format(/,5x,'regular function l = ',1x,i3)
 1000 format( (/,5x,5(e15.8,1x) ) )
 1010 format(/,5x,'derivative regular function l = ',1x,i3)
 1020 format(/,5x,'irregular function l = ',1x,i3)
 1030 format(/,5x,'derivative irregular function l = ',1x,i3)
      return
      end
*deck snorm
c***begin prologue     snorm
c***date written       861108   (yymmdd)
c***revision date      yymmdd   (yymmdd)
c***keywords           m1104, link 1104, wighted norm
c***author             schneider, barry (lanl)
c***source             m1104
c***purpose            scalar product in r space.
c***description        norm of two vectors including integration weight
c***                   in definition of scalar product.
c
c***references         none
c
c***routines called
c***end prologue       snorm
      function snorm (f1,f2,wt,n,nowgt)
      real *8 f1, f2, wt, sum, snorm
      logical nowgt
      dimension f1(n), f2(n), wt(n)
      sum=0.d+00
      if (.not.nowgt) then
         do 10 i=1,n
            sum=sum+f1(i)*wt(i)*f2(i)
   10    continue
      else
          do 20 i=1,n
             sum=sum+f1(i)*f2(i)
   20     continue
      endif
      snorm=sum
      return
      end
