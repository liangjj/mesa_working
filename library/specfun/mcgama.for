        program mcgama
c
c       ==========================================================
c       purpose: this program computes the gamma function â(z)  
c                or ln[â(z)] for a complex argument using 
c                subroutine cgama
c       input :  x  --- real part of z
c                y  --- imaginary part of z
c                kf --- function code
c                       kf=0 for ln[â(z)]
c                       kf=1 for â(z)
c       output:  gr --- real part of ln[â(z)] or â(z)
c                gi --- imaginary part of ln[â(z)] or â(z)
c       examples:
c
c         x         y           re[â(z)]           im[â(z)]
c       --------------------------------------------------------
c        2.50      5.00     .2267360319d-01    -.1172284404d-01
c        5.00     10.00     .1327696517d-01     .3639011746d-02
c        2.50     -5.00     .2267360319d-01     .1172284404d-01
c        5.00    -10.00     .1327696517d-01    -.3639011746d-02
c
c         x         y          re[lnâ(z)]         im[lnâ(z)]
c      ---------------------------------------------------------
c        2.50      5.00    -.3668103262d+01     .5806009801d+01
c        5.00     10.00    -.4285507444d+01     .1911707090d+02
c        2.50     -5.00    -.3668103262d+01    -.5806009801d+01
c        5.00    -10.00    -.4285507444d+01    -.1911707090d+02
c       ==========================================================
c
        double precision x,y,gr,gi
        write(*,*)'  please enter kf, x and y'
        read(*,*)kf,x,y
        write(*,*)
        if (kf.eq.1) then
            write(*,*)'       x         y           re[â(z)]',
     &                '           im[â(z)]'
        else
            write(*,*)'       x         y          re[lnâ(z)]',
     &                '         im[lnâ(z)]'
        endif
        write(*,*)'    ------------------------------------',
     &            '---------------------'
        call cgama(x,y,kf,gr,gi)
        write(*,10)x,y,gr,gi
10      format(1x,2f10.2,2d20.10)
        end


        subroutine cgama(x,y,kf,gr,gi)
c
c       =========================================================
c       purpose: compute the gamma function â(z) or ln[â(z)]
c                for a complex argument
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
