        program mitsl0
c
c       ===========================================================
c       purpose: this program evaluates the integral of modified 
c                struve function l0(t) with respect to t from 0  
c                to x using subroutine itsl0
c       input :  x   --- upper limit  ( x ò 0 )
c       output:  tl0 --- integration of l0(t) from 0 to x
c       example:
c                      x        l0(t)dt
c                   -----------------------
c                     0.0    .0000000d+00
c                     5.0    .3003079d+02
c                    10.0    .2990773d+04
c                    15.0    .3526179d+06
c                    20.0    .4475860d+08
c                    30.0    .7955389d+12
c                    40.0    .1508972d+17
c                    50.0    .2962966d+21
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x        l0(t)dt'
        write(*,*)'-----------------------'
        call itsl0(x,tl0)
        write(*,10)x,tl0
10      format(1x,f5.1,d16.7)
        end


        subroutine itsl0(x,tl0)
c
c       ===========================================================
c       purpose: evaluate the integral of modified struve function
c                l0(t) with respect to t from 0 to x
c       input :  x   --- upper limit  ( x ò 0 )
c       output:  tl0 --- integration of l0(t) from 0 to x
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        dimension a(18)
        pi=3.141592653589793d0
        r=1.0d0
        if (x.le.20.0) then
           s=0.5d0
           do 10 k=1,100
              rd=1.0d0
              if (k.eq.1) rd=0.5d0
              r=r*rd*k/(k+1.0d0)*(x/(2.0d0*k+1.0d0))**2
              s=s+r
              if (dabs(r/s).lt.1.0d-12) go to 15
10         continue
15         tl0=2.0d0/pi*x*x*s
        else
           s=1.0d0
           do 20 k=1,10
              r=r*k/(k+1.0d0)*((2.0d0*k+1.0d0)/x)**2
              s=s+r
              if (dabs(r/s).lt.1.0d-12) go to 25
20            continue
25         el=.57721566490153d0
           s0=-s/(pi*x*x)+2.0d0/pi*(dlog(2.0d0*x)+el)
           a0=1.0d0
           a1=5.0d0/8.0d0
           a(1)=a1
           do 30 k=1,10
              af=((1.5d0*(k+.50d0)*(k+5.0d0/6.0d0)*a1-.5d0*
     &            (k+.5d0)**2*(k-.5d0)*a0))/(k+1.0d0)
              a(k+1)=af
              a0=a1
30            a1=af
           ti=1.0d0
           r=1.0d0
           do 35 k=1,11
              r=r/x
35            ti=ti+a(k)*r
           tl0=ti/dsqrt(2*pi*x)*dexp(x)+s0
        endif
        return
        end
