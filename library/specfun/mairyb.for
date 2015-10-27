        program mairyb
c
c       ============================================================
c       purpose: this program computes airy functions and their 
c                derivatives using subroutine airyb
c       input:   x  --- argument of airy function
c       output:  ai --- ai(x)
c                bi --- bi(x)
c                ad --- ai'(x)
c                bd --- bi'(x)
c       example:
c
c   x       ai(x)          bi(x)          ai'(x)         bi'(x)
c  ----------------------------------------------------------------
c   0   .35502805d+00  .61492663d+00 -.25881940d+00  .44828836d+00
c  10   .11047533d-09  .45564115d+09 -.35206337d-09  .14292361d+10
c  20   .16916729d-26  .21037650d+26 -.75863916d-26  .93818393d+26
c  30   .32082176d-48  .90572885d+47 -.17598766d-47  .49533045d+48
c
c   x       ai(-x)         bi(-x)         ai'(-x)        bi'(-x)
c  ----------------------------------------------------------------
c   0       .35502805      .61492663     -.25881940      .44828836
c  10       .04024124     -.31467983      .99626504      .11941411
c  20      -.17640613     -.20013931      .89286286     -.79142903
c  30      -.08796819     -.22444694     1.22862060     -.48369473
c       ============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        call airyb(x,ai,bi,ad,bd)
        write(*,30)
        write(*,40)
        write(*,10)x,ai,bi,ad,bd
        write(*,*)
        call airyb(-x,ai,bi,ad,bd)
        write(*,50)
        write(*,40)
        write(*,20)x,ai,bi,ad,bd
10      format(1x,f5.1,4d16.8)
20      format(1x,f5.1,4d16.8)
30      format(4x,'x',8x,'ai(x)',11x,'bi(x)',11x,'ai''(x)',
     &         10x,'bi''(x)')
40      format(2x,'----------------------------------',
     &        '-----------------------------------')
50      format(4x,'x',8x,'ai(-x)',10x,'bi(-x)',10x,
     &        'ai''(-x)',9x,'bi''(-x)')
        end


        subroutine airyb(x,ai,bi,ad,bd)
c
c       =======================================================
c       purpose: compute airy functions and their derivatives
c       input:   x  --- argument of airy function
c       output:  ai --- ai(x)
c                bi --- bi(x)
c                ad --- ai'(x)
c                bd --- bi'(x)
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension ck(41),dk(41)
        eps=1.0d-15
        pi=3.141592653589793d0
        c1=0.355028053887817d0
        c2=0.258819403792807d0
        sr3=1.732050807568877d0
        xa=dabs(x)
        xq=dsqrt(xa)
        if (x.gt.0.0d0) xm=5.0
        if (x.le.0.0d0) xm=8.0
        if (x.eq.0.0d0) then
           ai=c1
           bi=sr3*c1
           ad=-c2
           bd=sr3*c2
           return
        endif
        if (xa.le.xm) then
           fx=1.0d0
           r=1.0d0
           do 10 k=1,40
              r=r*x/(3.0d0*k)*x/(3.0d0*k-1.0d0)*x
              fx=fx+r
              if (dabs(r).lt.dabs(fx)*eps) go to 15
10         continue
15         gx=x
           r=x
           do 20 k=1,40
              r=r*x/(3.0d0*k)*x/(3.0d0*k+1.0d0)*x
              gx=gx+r
              if (dabs(r).lt.dabs(gx)*eps) go to 25
20         continue
25         ai=c1*fx-c2*gx
           bi=sr3*(c1*fx+c2*gx)
           df=0.5d0*x*x
           r=df
           do 30 k=1,40
              r=r*x/(3.0d0*k)*x/(3.0d0*k+2.0d0)*x
              df=df+r
              if (dabs(r).lt.dabs(df)*eps) go to 35
30         continue
35         dg=1.0d0
           r=1.0d0
           do 40 k=1,40
              r=r*x/(3.0d0*k)*x/(3.0d0*k-2.0d0)*x
              dg=dg+r
              if (dabs(r).lt.dabs(dg)*eps) go to 45
40         continue
45         ad=c1*df-c2*dg
           bd=sr3*(c1*df+c2*dg)
        else
           xe=xa*xq/1.5d0
           xr1=1.0d0/xe
           xar=1.0d0/xq
           xf=dsqrt(xar)
           rp=0.5641895835477563d0
           r=1.0d0
           do 50 k=1,40
              r=r*(6.0d0*k-1.0d0)/216.0d0*(6.0d0*k-3.0d0)
     &          /k*(6.0d0*k-5.0d0)/(2.0d0*k-1.0d0)
              ck(k)=r
50            dk(k)=-(6.0d0*k+1.0d0)/(6.0d0*k-1.0d0)*ck(k)
           km=int(24.5-xa)
           if (xa.lt.6.0) km=14
           if (xa.gt.15.0) km=10
           if (x.gt.0.0d0) then
              sai=1.0d0
              sad=1.0d0
              r=1.0d0
              do 55 k=1,km
                 r=-r*xr1
                 sai=sai+ck(k)*r
55               sad=sad+dk(k)*r
              sbi=1.0d0
              sbd=1.0d0
              r=1.0d0
              do 60 k=1,km
                 r=r*xr1
                 sbi=sbi+ck(k)*r
60               sbd=sbd+dk(k)*r
              xp1=dexp(-xe)
              ai=0.5d0*rp*xf*xp1*sai
              bi=rp*xf/xp1*sbi
              ad=-.5d0*rp/xf*xp1*sad
              bd=rp/xf/xp1*sbd
           else
              xcs=dcos(xe+pi/4.0d0)
              xss=dsin(xe+pi/4.0d0)
              ssa=1.0d0
              sda=1.0d0
              r=1.0d0
              xr2=1.0d0/(xe*xe)
              do 65 k=1,km
                 r=-r*xr2
                 ssa=ssa+ck(2*k)*r
65               sda=sda+dk(2*k)*r
              ssb=ck(1)*xr1
              sdb=dk(1)*xr1
              r=xr1
              do 70 k=1,km
                 r=-r*xr2
                 ssb=ssb+ck(2*k+1)*r
70               sdb=sdb+dk(2*k+1)*r
              ai=rp*xf*(xss*ssa-xcs*ssb)
              bi=rp*xf*(xcs*ssa+xss*ssb)
              ad=-rp/xf*(xcs*sda+xss*sdb)
              bd=rp/xf*(xss*sda-xcs*sdb)
           endif
        endif
        return
        end
