        program mstvl1
c
c       =====================================================
c       purpose: this program computes the modified struve 
c                function l1(x) using subroutine stvl1
c       input :  x   --- argument of l1(x) ( x ò 0 )
c       output:  sl1 --- l1(x)
c       example:
c                     x        l1(x)
c                 -----------------------
c                   0.0   .00000000d+00
c                   5.0   .23728216d+02
c                  10.0   .26703583d+04
c                  15.0   .32812429d+06
c                  20.0   .42454973d+08
c                  30.0   .76853204d+12
c                  40.0   .14707396d+17
c                  50.0   .29030786d+21
c       =====================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x        l1(x)'
        write(*,*)'-----------------------'
        call stvl1(x,sl1)
        write(*,10)x,sl1
10      format(1x,f5.1,d16.8)
        end


        subroutine stvl1(x,sl1)
c
c       ================================================
c       purpose: compute modified struve function l1(x)
c       input :  x   --- argument of l1(x) ( x ò 0 )
c       output:  sl1 --- l1(x)
c       ================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        r=1.0d0
        if (x.le.20.0d0) then
           s=0.0d0
           do 10 k=1,60
              r=r*x*x/(4.0d0*k*k-1.0d0)
              s=s+r
              if (dabs(r/s).lt.1.0d-12) go to 15
10         continue
15         sl1=2.0d0/pi*s
        else
           s=1.0d0
           km=int(.50*x)
           if (x.gt.50) km=25
           do 20 k=1,km
              r=r*(2.0d0*k+3.0d0)*(2.0d0*k+1.0d0)/(x*x)
              s=s+r
              if (dabs(r/s).lt.1.0d-12) go to 25
20            continue
25         sl1=2.0d0/pi*(-1.0d0+1.0d0/(x*x)+3.0d0*s/x**4)
           a1=dexp(x)/dsqrt(2.0d0*pi*x)
           r=1.0d0
           bi1=1.0d0
           do 30 k=1,16
              r=-0.125d0*r*(4.0d0-(2.0d0*k-1.0d0)**2)/(k*x)
              bi1=bi1+r
              if (dabs(r/bi1).lt.1.0d-12) go to 35
30         continue
35         sl1=sl1+a1*bi1
        endif
        return
        end
