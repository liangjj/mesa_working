        program me1xa
c
c       =========================================================
c       purpose: this program computes the exponential integral 
c                e1(x) using subroutine e1xa
c       input :  x  --- argument of e1(x)  ( x > 0 )
c       output:  e1 --- e1(x)
c       example:
c                  x        e1(x)
c                ----------------------
c                 0.0     .1000000+301
c                 1.0     .2193839e+00
c                 2.0     .4890051e-01
c                 3.0     .1304838e-01
c                 4.0     .3779352e-02
c                 5.0     .1148296e-02
c       =========================================================
c
        double precision e1,x
        write(*,*)'please enter x '
        read(*,*) x
        write(*,*)'   x        e1(x)'
        write(*,*)' ----------------------'
        call e1xa(x,e1)
        write(*,10)x,e1
10      format(1x,f5.1,e17.7)
        end


        subroutine e1xa(x,e1)
c
c       ============================================
c       purpose: compute exponential integral e1(x)
c       input :  x  --- argument of e1(x) 
c       output:  e1 --- e1(x) ( x > 0 )
c       ============================================
c
        implicit double precision (a-h,o-z)
        if (x.eq.0.0) then
           e1=1.0d+300
        else if (x.le.1.0) then
           e1=-dlog(x)+((((1.07857d-3*x-9.76004d-3)*x+5.519968d-2)*x
     &        -0.24991055d0)*x+0.99999193d0)*x-0.57721566d0
        else
           es1=(((x+8.5733287401d0)*x+18.059016973d0)*x
     &         +8.6347608925d0)*x+0.2677737343d0
           es2=(((x+9.5733223454d0)*x+25.6329561486d0)*x
     &         +21.0996530827d0)*x+3.9584969228d0
           e1=dexp(-x)/x*es1/es2
        endif
        return
        end
