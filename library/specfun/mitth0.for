        program mitth0
c
c       ===========================================================
c       purpose: this program evaluates the integral of h0(t)/t 
c                with respect to t from x to infinity using
c                subroutine itth0
c       input :  x   --- lower limit  ( x ò 0 )
c       output:  tth --- integration of h0(t)/t from x to infinity
c       example:
c                    x        h0(t)/t dt
c                 -----------------------
c                   0.0      1.57079633
c                   5.0       .07954575
c                  10.0       .04047175
c                  15.0       .04276558
c                  20.0       .04030796
c                  30.0       .01815256
c                  40.0       .01621331
c                  50.0       .01378661
c       =======================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*) '   x       h0(t)/t dt'
        write(*,*)'-----------------------'
        call itth0(x,tth)
        write(*,10)x,tth
10      format(1x,f5.1,1x,e16.8)
        end


        subroutine itth0(x,tth)
c
c       ===========================================================
c       purpose: evaluate the integral h0(t)/t with respect to t
c                from x to infinity
c       input :  x   --- lower limit  ( x ò 0 )
c       output:  tth --- integration of h0(t)/t from x to infinity
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        s=1.0d0
        r=1.0d0
        if (x.lt.24.5d0) then
           do 10 k=1,60
              r=-r*x*x*(2.0*k-1.0d0)/(2.0*k+1.0d0)**3
              s=s+r
              if (dabs(r).lt.dabs(s)*1.0d-12) go to 15
10         continue
15         tth=pi/2.0d0-2.0d0/pi*x*s
        else
           do 20 k=1,10
              r=-r*(2.0*k-1.0d0)**3/((2.0*k+1.0d0)*x*x)
              s=s+r
              if (dabs(r).lt.dabs(s)*1.0d-12) go to 25
20            continue
25         tth=2.0d0/(pi*x)*s
           t=8.0d0/x
           xt=x+.25d0*pi
           f0=(((((.18118d-2*t-.91909d-2)*t+.017033d0)*t
     &        -.9394d-3)*t-.051445d0)*t-.11d-5)*t+.7978846d0
           g0=(((((-.23731d-2*t+.59842d-2)*t+.24437d-2)*t
     &        -.0233178d0)*t+.595d-4)*t+.1620695d0)*t
           tty=(f0*dsin(xt)-g0*dcos(xt))/(dsqrt(x)*x)
           tth=tth+tty
        endif
        return
        end
