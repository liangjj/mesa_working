        program me1xb
c
c       =========================================================
c       purpose: this program computes the exponential integral 
c                e1(x) using subroutine e1xb
c       input :  x  --- argument of e1(x)  ( x > 0 )
c       output:  e1 --- e1(x)
c       example:
c                  x          e1(x)
c                -------------------------
c                 0.0     .1000000000+301
c                 1.0     .2193839344e+00
c                 2.0     .4890051071e-01
c                 3.0     .1304838109e-01
c                 4.0     .3779352410e-02
c                 5.0     .1148295591e-02
c       =========================================================
c
        double precision e1,x
        write(*,*)'please enter x '
        read(*,*) x
        write(*,*)'   x          e1(x)'
        write(*,*)' -------------------------'
        call e1xb(x,e1)
        write(*,10)x,e1
10      format(1x,f5.1,e20.10)
        end


        subroutine e1xb(x,e1)
c
c       ============================================
c       purpose: compute exponential integral e1(x)
c       input :  x  --- argument of e1(x)
c       output:  e1 --- e1(x)  ( x > 0 )
c       ============================================
c
        implicit double precision (a-h,o-z)
        if (x.eq.0.0) then
           e1=1.0d+300
        else if (x.le.1.0) then
           e1=1.0d0
           r=1.0d0
           do 10 k=1,25
              r=-r*k*x/(k+1.0d0)**2
              e1=e1+r
              if (dabs(r).le.dabs(e1)*1.0d-15) go to 15
10         continue
15         ga=0.5772156649015328d0
           e1=-ga-dlog(x)+x*e1
        else
           m=20+int(80.0/x)
           t0=0.0d0
           do 20 k=m,1,-1
              t0=k/(1.0d0+k/(x+t0))
20         continue
           t=1.0d0/(x+t0)
           e1=dexp(-x)*t
        endif
        return
        end
