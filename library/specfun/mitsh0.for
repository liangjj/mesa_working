        program mitsh0
c
c       ====================================================
c       purpose: this program evaluates the integral of 
c                struve function h0(t) with respect to t 
c                from 0 and x using subroutine itsh0
c       input :  x   --- upper limit  ( x ò 0 )
c       output:  th0 --- integration of h0(t) from 0 and x
c       example:
c                    x        h0(t)dt
c                 ----------------------
c                   0.0       .0000000
c                   5.0      2.0442437
c                  10.0      2.5189577
c                  15.0      2.5415824
c                  20.0      2.5484517
c                  30.0      3.0625848
c                  40.0      3.1484123
c                  50.0      3.2445168
c       ====================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x         h0(t)dt'
        write(*,*)'----------------------Ä'
        call itsh0(x,th0)
        write(*,10)x,th0
10      format(1x,f5.1,e16.7)
        end


        subroutine itsh0(x,th0)
c
c       ===================================================
c       purpose: evaluate the integral of struve function
c                h0(t) with respect to t from 0 and x
c       input :  x   --- upper limit  ( x ò 0 )
c       output:  th0 --- integration of h0(t) from 0 and x
c       ===================================================
c
        implicit double precision (a-h,o-z)
        dimension a(25)
        pi=3.141592653589793d0
        r=1.0d0            
        if (x.le.30.0) then
           s=0.5d0
           do 10 k=1,100
              rd=1.0d0
              if (k.eq.1) rd=0.5d0
              r=-r*rd*k/(k+1.0d0)*(x/(2.0d0*k+1.0d0))**2
              s=s+r
              if (dabs(r).lt.dabs(s)*1.0d-12) go to 15
10         continue
15         th0=2.0d0/pi*x*x*s
        else
           s=1.0d0
           do 20 k=1,12
              r=-r*k/(k+1.0d0)*((2.0d0*k+1.0d0)/x)**2
              s=s+r
              if (dabs(r).lt.dabs(s)*1.0d-12) go to 25
20         continue
25         el=.57721566490153d0
           s0=s/(pi*x*x)+2.0d0/pi*(dlog(2.0d0*x)+el)
           a0=1.0d0
           a1=5.0d0/8.0d0
           a(1)=a1
           do 30 k=1,20
              af=((1.5d0*(k+.5d0)*(k+5.0d0/6.0d0)*a1-.5d0
     &           *(k+.5d0)*(k+.5d0)*(k-.5d0)*a0))/(k+1.0d0)
              a(k+1)=af
              a0=a1
30            a1=af
           bf=1.0d0
           r=1.0d0
           do 35 k=1,10
              r=-r/(x*x)
35            bf=bf+a(2*k)*r
           bg=a(1)/x
           r=1.0d0/x
           do 40 k=1,10
              r=-r/(x*x)
40            bg=bg+a(2*k+1)*r
           xp=x+.25d0*pi
           ty=dsqrt(2.0d0/(pi*x))*(bg*dcos(xp)-bf*dsin(xp))
           th0=ty+s0
        endif
        return
        end
