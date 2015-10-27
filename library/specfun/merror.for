        program merror
c
c       ===================================================
c       purpose: this program computes the error function 
c                erf(x) using subroutine error
c       input:   x   --- argument of erf(x)
c       output:  err --- erf(x)
c       example:
c                  x         erf(x)
c                ---------------------
c                 1.0       .84270079
c                 2.0       .99532227
c                 3.0       .99997791
c                 4.0       .99999998
c                 5.0      1.00000000
c       ===================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x         erf(x)'
        write(*,*)'---------------------'
        call error(x,er)
        write(*,10)x,er
10      format(1x,f5.2,f15.8)
        end


        subroutine error(x,err)
c
c       =========================================
c       purpose: compute error function erf(x)
c       input:   x   --- argument of erf(x)
c       output:  err --- erf(x)
c       =========================================
c
        implicit double precision (a-h,o-z)
        eps=1.0d-15
        pi=3.141592653589793d0
        x2=x*x
        if (dabs(x).lt.3.5d0) then
           er=1.0d0
           r=1.0d0
           do 10 k=1,50
              r=r*x2/(k+0.5d0)
              er=er+r
              if (dabs(r).le.dabs(er)*eps) go to 15
10         continue
15         c0=2.0d0/dsqrt(pi)*x*dexp(-x2)
           err=c0*er
        else
           er=1.0d0
           r=1.0d0
           do 20 k=1,12
              r=-r*(k-0.5d0)/x2
20            er=er+r
           c0=dexp(-x2)/(dabs(x)*dsqrt(pi))
           err=1.0d0-c0*er
           if (x.lt.0.0) err=-err
        endif
        return
        end
