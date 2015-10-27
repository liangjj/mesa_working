        program mpbwa
c
c       ============================================================
c       purpose: this program computes the parabolic cylinder 
c                functions w(a,ñx) and their derivatives using
c                subroutine pbwa
c       input  : a --- parameter  ( 0 ó |a| ó 5 )
c                x --- argument of w(a,ñx)  ( 0 ó |x| ó 5 )
c       output : w1f --- w(a,x)
c                w1d --- w'(a,x)
c                w2f --- w(a,-x)
c                w2d --- w'(a,-x)
c       example: x = 5.0
c                 a      w(a,x)     w'(a,x)    w(a,-x)   w'(a,-x)
c              ----------------------------------------------------
c                0.5   .1871153    .1915744  -.8556585   4.4682493
c                1.5  -.0215853    .0899870 -8.8586002  -9.3971967
c                0.0   .3009549   -.7148233   .6599634   1.7552224
c               -0.5  -.1934088  -1.3474400   .6448148   -.6781011
c               -1.5  -.5266539    .8219516  -.2822774  -1.4582283
c               -5.0   .0893618  -1.8118641   .5386084    .2698553
c       ============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter a and x '
        read(*,*)a,x
        write(*,10)a,x
        write(*,*)
        write(*,*)'   a       w(a,x)          w''(a,x)',
     &            '         w(a,-x)         w''(a,-x)'
        write(*,*)' -----------------------------------------',
     &            '----------------------------'
        call pbwa(a,x,w1f,w1d,w2f,w2d)
        write(*,20)a,w1f,w1d,w2f,w2d
10      format(1x,'a=',f5.1,3x,'x=',f5.1)
20      format(1x,f5.1,4d16.8)
        end


        subroutine pbwa(a,x,w1f,w1d,w2f,w2d)
c
c       ======================================================
c       purpose: compute parabolic cylinder functions w(a,ñx)
c                and their derivatives
c       input  : a --- parameter  ( 0 ó |a| ó 5 )
c                x --- argument of w(a,ñx)  ( 0 ó |x| ó 5 )
c       output : w1f --- w(a,x)
c                w1d --- w'(a,x)
c                w2f --- w(a,-x)
c                w2d --- w'(a,-x)
c       routine called:
c               cgama for computing complex gamma function
c       ======================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex *16 (c,z)
        dimension h(100),d(100)
        eps=1.0d-15
        p0=0.59460355750136d0
        if (a.eq.0.0d0) then
           g1=3.625609908222d0
           g2=1.225416702465d0
        else
           x1=0.25d0
           y1=0.5d0*a
           call cgama(x1,y1,1,ugr,ugi)
           g1=dsqrt(ugr*ugr+ugi*ugi)
           x2=0.75d0
           call cgama(x2,y1,1,vgr,vgi)
           g2=dsqrt(vgr*vgr+vgi*vgi)
        endif
        f1=dsqrt(g1/g2)
        f2=dsqrt(2.0d0*g2/g1)
        h0=1.0d0
        h1=a
        h(1)=a
        do 10 l1=4,200,2
           m=l1/2
           hl=a*h1-0.25d0*(l1-2.0d0)*(l1-3.0d0)*h0
           h(m)=hl
           h0=h1
10         h1=hl
        y1f=1.0d0
        r=1.0d0
        do 15 k=1,100
           r=0.5d0*r*x*x/(k*(2.0d0*k-1.0d0))
           r1=h(k)*r
           y1f=y1f+r1
           if (dabs(r1/y1f).le.eps.and.k.gt.30) go to 20
15      continue
20      y1d=a
        r=1.0d0
        do 25 k=1,100
           r=0.5d0*r*x*x/(k*(2.0d0*k+1.0d0))
           r1=h(k+1)*r
           y1d=y1d+r1
           if (dabs(r1/y1d).le.eps.and.k.gt.30) go to 30
25      continue
30      y1d=x*y1d
        d1=1.0d0
        d2=a
        d(1)=1.0d0
        d(2)=a
        do 40 l2=5,160,2
           m=(l2+1)/2
           dl=a*d2-0.25d0*(l2-2.0d0)*(l2-3.0d0)*d1
           d(m)=dl
           d1=d2
40         d2=dl
        y2f=1.0d0
        r=1.0d0
        do 45 k=1,100
           r=0.5d0*r*x*x/(k*(2.0d0*k+1.0d0))
           r1=d(k+1)*r
           y2f=y2f+r1
           if (dabs(r1/y2f).le.eps.and.k.gt.30) go to 50
45      continue
50      y2f=x*y2f
        y2d=1.0d0
        r=1.0d0
        do 55 k=1,100
           r=0.5d0*r*x*x/(k*(2.0d0*k-1.0d0))
           r1=d(k+1)*r
           y2d=y2d+r1
           if (dabs(r1/y2d).le.eps.and.k.gt.30) go to 60
55      continue
60      w1f=p0*(f1*y1f-f2*y2f)
        w2f=p0*(f1*y1f+f2*y2f)
        w1d=p0*(f1*y1d-f2*y2d)
        w2d=p0*(f1*y1d+f2*y2d)
        return
        end


        subroutine cgama(x,y,kf,gr,gi)
c
c       =========================================================
c       purpose: compute complex gamma function â(z) or ln[â(z)]
c       input :  x  --- real part of z
c                y  --- imaginary part of z
c                kf --- function code
c                       kf=0 for ln[â(z)]
c                       kf=1 for â(z)
c       output:  gr --- real part of ln[â(z)] or â(z)
c                gi --- imaginary part of ln[â(z)] or â(z)
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension a(10)
        pi=3.141592653589793d0
        data a/8.333333333333333d-02,-2.777777777777778d-03,
     &         7.936507936507937d-04,-5.952380952380952d-04,
     &         8.417508417508418d-04,-1.917526917526918d-03,
     &         6.410256410256410d-03,-2.955065359477124d-02,
     &         1.796443723688307d-01,-1.39243221690590d+00/
        if (y.eq.0.0d0.and.x.eq.int(x).and.x.le.0.0d0) then
           gr=1.0d+300
           gi=0.0d0
           return
        else if (x.lt.0.0d0) then
           x1=x
           y1=y
           x=-x
           y=-y
        endif
        x0=x
        if (x.le.7.0) then
           na=int(7-x)
           x0=x+na
        endif
        z1=dsqrt(x0*x0+y*y)
        th=datan(y/x0)
        gr=(x0-.5d0)*dlog(z1)-th*y-x0+0.5d0*dlog(2.0d0*pi)
        gi=th*(x0-0.5d0)+y*dlog(z1)-y
        do 10 k=1,10
           t=z1**(1-2*k)
           gr=gr+a(k)*t*dcos((2.0d0*k-1.0d0)*th)
10         gi=gi-a(k)*t*dsin((2.0d0*k-1.0d0)*th)
        if (x.le.7.0) then
           gr1=0.0d0
           gi1=0.0d0
           do 15 j=0,na-1
              gr1=gr1+.5d0*dlog((x+j)**2+y*y)
15            gi1=gi1+datan(y/(x+j))
           gr=gr-gr1
           gi=gi-gi1
        endif
        if (x1.lt.0.0d0) then
           z1=dsqrt(x*x+y*y)
           th1=datan(y/x)
           sr=-dsin(pi*x)*dcosh(pi*y)
           si=-dcos(pi*x)*dsinh(pi*y)
           z2=dsqrt(sr*sr+si*si)
           th2=datan(si/sr)
           if (sr.lt.0.0d0) th2=pi+th2
           gr=dlog(pi/(z1*z2))-gr
           gi=-th1-th2-gi
           x=x1
           y=y1
        endif
        if (kf.eq.1) then
           g0=dexp(gr)
           gr=g0*dcos(gi)
           gi=g0*dsin(gi)
        endif
        return
        end
