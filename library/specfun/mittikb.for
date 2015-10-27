        program mittikb
c
c       ============================================================
c       purpose: this program computes the integral of [i0(t)-1]/t
c                with respect to t from 0 to x and k0(t)/t with 
c                respect to t from x to ì using subroutine ittikb
c       input :  x   --- upper limit of the integral
c       output:  tti --- integration of [i0(t)-1]/t from 0 to x
c                ttk --- integration of k0(t)/t from x to ì
c       example:
c                   x     [1-i0(t)]/tdt      k0(t)/tdt
c                ---------------------------------------
c                  5.0     .710478d+01     .586361d-03
c                 10.0     .340811d+03     .156293d-05
c                 15.0     .254373d+05     .598363d-08
c                 20.0     .236735d+07     .267906d-10
c                 25.0     .246534d+09     .131007d-12
c       ============================================================
c
        double precision x,tti,ttk
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x     [1-i0(t)]/tdt      k0(t)/tdt'
        write(*,*)'---------------------------------------'
        call ittikb(x,tti,ttk)
        write(*,10)x,tti,ttk
10      format(1x,f5.1,2d16.6)
        end


        subroutine ittikb(x,tti,ttk)
c
c       =========================================================
c       purpose: integrate [i0(t)-1]/t with respect to t from 0
c                to x, and k0(t)/t with respect to t from x to ì
c       input :  x   --- variable in the limits  ( x ò 0 )
c       output:  tti --- integration of [i0(t)-1]/t from 0 to x
c                ttk --- integration of k0(t)/t from x to ì
c       =========================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        if (x.eq.0.0d0) then
           tti=0.0d0
        else if (x.le.5.0d0) then
           x1=x/5.0d0
           t=x1*x1
           tti=(((((((.1263d-3*t+.96442d-3)*t+.968217d-2)*t
     &         +.06615507d0)*t+.33116853d0)*t+1.13027241d0)
     &         *t+2.44140746d0)*t+3.12499991d0)*t
        else
           t=5.0d0/x
           tti=(((((((((2.1945464d0*t-3.5195009d0)*t
     &         -11.9094395d0)*t+40.394734d0)*t-48.0524115d0)
     &         *t+28.1221478d0)*t-8.6556013d0)*t+1.4780044d0)
     &         *t-.0493843d0)*t+.1332055d0)*t+.3989314d0
           tti=tti*dexp(x)/(dsqrt(x)*x)
        endif
        if (x.eq.0.0d0) then
           ttk=1.0d+300
        else if (x.le.2.0d0) then
           t1=x/2.0d0
           t=t1*t1
           ttk=(((((.77d-6*t+.1544d-4)*t+.48077d-3)*t
     &         +.925821d-2)*t+.10937537d0)*t+.74999993d0)*t
           e0=el+dlog(x/2.0d0)
           ttk=pi*pi/24.0d0+e0*(.5d0*e0+tti)-ttk
        else if (x.le.4.0d0) then
           t=2.0d0/x
           ttk=(((.06084d0*t-.280367d0)*t+.590944d0)*t
     &         -.850013d0)*t+1.234684d0
           ttk=ttk*dexp(-x)/(dsqrt(x)*x)
        else
           t=4.0d0/x
           ttk=(((((.02724d0*t-.1110396d0)*t+.2060126d0)*t
     &         -.2621446d0)*t+.3219184d0)*t-.5091339d0)*t
     &         +1.2533141d0
           ttk=ttk*dexp(-x)/(dsqrt(x)*x)
        endif
        return
        end
