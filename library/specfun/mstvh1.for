        program mstvh1
c
c       =====================================================
c       purpose: this program computes struve function 
c                h1(x) using subroutine stvh1
c       input :  x   --- argument of h1(x) ( x ò 0 )
c       output:  sh1 --- h1(x)
c       example:
c                   x          h1(x)
c                -----------------------
c                  0.0       .00000000
c                  5.0       .80781195
c                 10.0       .89183249
c                 15.0       .66048730
c                 20.0       .47268818
c                 25.0       .53880362
c       =====================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x          h1(x)'
        write(*,*)'-----------------------'
        call stvh1(x,sh1)
        write(*,10)x,sh1
10      format(1x,f5.1,e16.8)
        end


        subroutine stvh1(x,sh1)
c
c       =============================================
c       purpose: compute struve function h1(x)
c       input :  x   --- argument of h1(x) ( x ò 0 )
c       output:  sh1 --- h1(x)
c       =============================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        r=1.0d0
        if (x.le.20.0d0) then
           s=0.0d0
           a0=-2.0d0/pi
           do 10 k=1,60
              r=-r*x*x/(4.0d0*k*k-1.0d0)
              s=s+r
              if (dabs(r).lt.dabs(s)*1.0d-12) go to 15
10         continue
15         sh1=a0*s
        else
           s=1.0d0
           km=int(.5*x)
           if (x.gt.50.d0) km=25
           do 20 k=1,km
              r=-r*(4.0d0*k*k-1.0d0)/(x*x)
              s=s+r
              if (dabs(r).lt.dabs(s)*1.0d-12) go to 25
20         continue
25         t=4.0d0/x
           t2=t*t
           p1=((((.42414d-5*t2-.20092d-4)*t2+.580759d-4)*t2
     &        -.223203d-3)*t2+.29218256d-2)*t2+.3989422819d0
           q1=t*(((((-.36594d-5*t2+.1622d-4)*t2-.398708d-4)*
     &        t2+.1064741d-3)*t2-.63904d-3)*t2+.0374008364d0)
           ta1=x-.75d0*pi
           by1=2.0d0/dsqrt(x)*(p1*dsin(ta1)+q1*dcos(ta1))
           sh1=2.0/pi*(1.0d0+s/(x*x))+by1
        endif
        return
        end
