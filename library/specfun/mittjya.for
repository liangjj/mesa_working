        program mittjya
c
c       ===========================================================
c       purpose: this program computes the integral of [1-j0(t)]/t 
c                with respect to t from 0 to x and y0(t)/t with 
c                respect to t from x to ì using subroutine ittjya
c       input :  x   --- variable in the limits  ( x ò 0 )
c       output:  ttj --- integration of [1-j0(t)]/t from 0 to x
c                tty --- integration of y0(t)/t from x to ì
c       example:
c                  x       [1-j0(t)]/tdt       y0(t)/tdt
c               -------------------------------------------
c                 5.0     .15403472d+01    -.46322055d-01
c                10.0     .21778664d+01    -.22987934d-01
c                15.0     .25785507d+01     .38573574d-03
c                20.0     .28773106d+01     .85031527d-02
c                25.0     .31082313d+01     .35263393d-02
c       ===========================================================
c
        double precision x,ttj,tty
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x      [1-j0(t)]/tdt       y0(t)/tdt'
        write(*,*)'-------------------------------------------'
        call ittjya(x,ttj,tty)
        write(*,10)x,ttj,tty
10      format(1x,f5.1,2d18.8)
        end


        subroutine ittjya(x,ttj,tty)
c
c       =========================================================
c       purpose: integrate [1-j0(t)]/t with respect to t from 0
c                to x, and y0(t)/t with respect to t from x to ì
c       input :  x   --- variable in the limits  ( x ò 0 )
c       output:  ttj --- integration of [1-j0(t)]/t from 0 to x
c                tty --- integration of y0(t)/t from x to ì
c       =========================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        if (x.eq.0.0d0) then
           ttj=0.0d0
           tty=-1.0d+300
        else if (x.le.20.0d0) then
           ttj=1.0d0
           r=1.0d0
           do 10 k=2,100
              r=-.25d0*r*(k-1.0d0)/(k*k*k)*x*x
              ttj=ttj+r
              if (dabs(r).lt.dabs(ttj)*1.0d-12) go to 15
10         continue
15         ttj=ttj*.125d0*x*x
           e0=.5d0*(pi*pi/6.0d0-el*el)-(.5d0*dlog(x/2.0d0)+el)
     &        *dlog(x/2.0d0)
           b1=el+dlog(x/2.0d0)-1.5d0
           rs=1.0d0
           r=-1.0d0
           do 20 k=2,100
              r=-.25d0*r*(k-1.0d0)/(k*k*k)*x*x
              rs=rs+1.0d0/k
              r2=r*(rs+1.0d0/(2.0d0*k)-(el+dlog(x/2.0d0)))
              b1=b1+r2
              if (dabs(r2).lt.dabs(b1)*1.0d-12) go to 25
20         continue
25         tty=2.0d0/pi*(e0+.125d0*x*x*b1)
        else
           a0=dsqrt(2.0d0/(pi*x))
           do 50 l=0,1
              vt=4.0d0*l*l
              px=1.0d0
              r=1.0d0
              do 30 k=1,14
                 r=-.0078125d0*r*(vt-(4.0d0*k-3.0d0)**2)
     &             /(x*k)*(vt-(4.0d0*k-1.0d0)**2)
     &             /((2.0d0*k-1.0d0)*x)
                 px=px+r
                 if (dabs(r).lt.dabs(px)*1.0d-12) go to 35
30            continue
35            qx=1.0d0
              r=1.0d0
              do 40 k=1,14
                 r=-.0078125d0*r*(vt-(4.0d0*k-1.0d0)**2)
     &             /(x*k)*(vt-(4.0d0*k+1.0d0)**2)
     &             /(2.0d0*k+1.0d0)/x
                 qx=qx+r
                 if (dabs(r).lt.dabs(qx)*1.0d-12) go to 45
40            continue
45            qx=.125d0*(vt-1.0d0)/x*qx
              xk=x-(.25d0+.5d0*l)*pi
              bj1=a0*(px*dcos(xk)-qx*dsin(xk))
              by1=a0*(px*dsin(xk)+qx*dcos(xk))
              if (l.eq.0) then
                 bj0=bj1
                 by0=by1
              endif
50         continue
           t=2.0d0/x
           g0=1.0d0
           r0=1.0d0
           do 55 k=1,10
              r0=-k*k*t*t*r0
55            g0=g0+r0
           g1=1.0d0
           r1=1.0d0
           do 60 k=1,10
              r1=-k*(k+1.0d0)*t*t*r1
60            g1=g1+r1
           ttj=2.0d0*g1*bj0/(x*x)-g0*bj1/x+el+dlog(x/2.0d0)
           tty=2.0d0*g1*by0/(x*x)-g0*by1/x
        endif
        return
        end
