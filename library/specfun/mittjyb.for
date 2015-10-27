        program mittjyb
c
c       ===========================================================
c       purpose: this program computes the integral of [1-j0(t)]/t 
c                with respect to t from 0 to x and y0(t)/t with 
c                respect to t from x to ì using subroutine ittjyb
c       input :  x   --- variable in the limits  ( x ò 0 )
c       output:  ttj --- integration of [1-j0(t)]/t from 0 to x
c                tty --- integration of y0(t)/t from x to ì
c       example:
c                  x      [1-j0(t)]/tdt       y0(t)/tdt
c                ----------------------------------------
c                 5.0     .1540347d+01    -.4632208d-01
c                10.0     .2177866d+01    -.2298791d-01
c                15.0     .2578551d+01     .3857453d-03
c                20.0     .2877311d+01     .8503154d-02
c                25.0     .3108231d+01     .3526339d-02
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x     [1-j0(t)]/tdt      y0(t)/tdt'
        write(*,*)'----------------------------------------'
        call ittjyb(x,ttj,tty)
        write(*,10)x,ttj,tty
10      format(1x,f5.1,2d17.7)
        end


        subroutine ittjyb(x,ttj,tty)
c
c       ==========================================================
c       purpose: integrate [1-j0(t)]/t with respect to t from 0
c                to x, and y0(t)/t with respect to t from x to ì
c       input :  x   --- variable in the limits  ( x ò 0 )
c       output:  ttj --- integration of [1-j0(t)]/t from 0 to x
c                tty --- integration of y0(t)/t from x to ì
c       ==========================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        if (x.eq.0.0d0) then
           ttj=0.0d0
           tty=-1.0d+300
        else if (x.le.4.0d0) then
           x1=x/4.0d0
           t=x1*x1
           ttj=((((((.35817d-4*t-.639765d-3)*t+.7092535d-2)*t
     &         -.055544803d0)*t+.296292677d0)*t-.999999326d0)
     &         *t+1.999999936d0)*t
           tty=(((((((-.3546d-5*t+.76217d-4)*t-.1059499d-2)*t
     &         +.010787555d0)*t-.07810271d0)*t+.377255736d0)
     &         *t-1.114084491d0)*t+1.909859297d0)*t
           e0=el+dlog(x/2.0d0)
           tty=pi/6.0d0+e0/pi*(2.0d0*ttj-e0)-tty
        else if (x.le.8.0d0) then
           xt=x+.25d0*pi
           t1=4.0d0/x
           t=t1*t1
           f0=(((((.0145369d0*t-.0666297d0)*t+.1341551d0)*t
     &        -.1647797d0)*t+.1608874d0)*t-.2021547d0)*t
     &        +.7977506d0
           g0=((((((.0160672d0*t-.0759339d0)*t+.1576116d0)*t
     &        -.1960154d0)*t+.1797457d0)*t-.1702778d0)*t
     &        +.3235819d0)*t1
           ttj=(f0*dcos(xt)+g0*dsin(xt))/(dsqrt(x)*x)
           ttj=ttj+el+dlog(x/2.0d0)
           tty=(f0*dsin(xt)-g0*dcos(xt))/(dsqrt(x)*x)
        else
           t=8.0d0/x
           xt=x+.25d0*pi
           f0=(((((.18118d-2*t-.91909d-2)*t+.017033d0)*t
     &        -.9394d-3)*t-.051445d0)*t-.11d-5)*t+.7978846d0
           g0=(((((-.23731d-2*t+.59842d-2)*t+.24437d-2)*t
     &      -.0233178d0)*t+.595d-4)*t+.1620695d0)*t
           ttj=(f0*dcos(xt)+g0*dsin(xt))/(dsqrt(x)*x)
     &         +el+dlog(x/2.0d0)
           tty=(f0*dsin(xt)-g0*dcos(xt))/(dsqrt(x)*x)
        endif
        return
        end
