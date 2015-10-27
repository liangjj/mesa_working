        program mstvl0
c
c       =====================================================
c       purpose: this program computes modified struve 
c                function l0(x) using subroutine stvl0
c       input :  x   --- argument of l0(x) ( x ò 0 )
c       output:  sl0 --- l0(x)
c       example:
c                   x        l0(x)
c               ------------------------
c                  0.0   .00000000d+00
c                  5.0   .27105917d+02
c                 10.0   .28156522d+04
c                 15.0   .33964933d+06
c                 20.0   .43558283d+08
c                 30.0   .78167230d+12
c                 40.0   .14894775d+17
c                 50.0   .29325538d+21
c       =====================================================
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x        l0(x)'
        write(*,*)'-----------------------'
        call stvl0(x,sl0)
        write(*,10)x,sl0
10      format(1x,f5.1,d16.8)
        end


        subroutine stvl0(x,sl0)
c
c       ================================================
c       purpose: compute modified struve function l0(x)
c       input :  x   --- argument of l0(x) ( x ò 0 )
c       output:  sl0 --- l0(x)
c       ================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        s=1.0d0
        r=1.0d0
        if (x.le.20.0d0) then
           a0=2.0d0*x/pi
           do 10 k=1,60
              r=r*(x/(2.0d0*k+1.0d0))**2
              s=s+r
              if (dabs(r/s).lt.1.0d-12) go to 15
10         continue
15         sl0=a0*s
        else
           km=int(.5*(x+1.0))
           if (x.ge.50.0) km=25
           do 20 k=1,km
              r=r*((2.0d0*k-1.0d0)/x)**2
              s=s+r
              if (dabs(r/s).lt.1.0d-12) go to 25
20         continue
25         a1=dexp(x)/dsqrt(2.0d0*pi*x)
           r=1.0d0
           bi0=1.0d0
           do 30 k=1,16
              r=0.125d0*r*(2.0d0*k-1.0d0)**2/(k*x)
              bi0=bi0+r
              if (dabs(r/bi0).lt.1.0d-12) go to 35
30         continue
35         bi0=a1*bi0
           sl0=-2.0d0/(pi*x)*s+bi0
        endif
        return
        end
