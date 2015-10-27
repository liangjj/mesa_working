        program mittika
c
c       ============================================================
c       purpose: this program computes the integral of [i0(t)-1]/t
c                with respect to t from 0 to x and k0(t)/t with 
c                respect to t from x to ì using subroutine ittika
c       input :  x   --- variable in the limits  ( x ò 0 )
c       output:  tti --- integration of [i0(t)-1]/t from 0 to x
c                ttk --- integration of k0(t)/t from x to ì
c       example:
c                   x     [1-i0(t)]/tdt     k0(t)/tdt
c                ---------------------------------------
c                  5.0   .71047763d+01   .58635626d-03
c                 10.0   .34081537d+03   .15629282d-05
c                 15.0   .25437619d+05   .59837472d-08
c                 20.0   .23673661d+07   .26790545d-10
c                 25.0   .24652751d+09   .13100706d-12
c       ============================================================
c
        double precision x,tti,ttk
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x     [1-i0(t)]/tdt     k0(t)/tdt'
        write(*,*)'---------------------------------------'
        call ittika(x,tti,ttk)
        write(*,10)x,tti,ttk
10      format(1x,f5.1,2d16.8)
        end


        subroutine ittika(x,tti,ttk)
c
c       =========================================================
c       purpose: integrate [i0(t)-1]/t with respect to t from 0
c                to x, and k0(t)/t with respect to t from x to ì
c       input :  x   --- variable in the limits  ( x ò 0 )
c       output:  tti --- integration of [i0(t)-1]/t from 0 to x
c                ttk --- integration of k0(t)/t from x to ì
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension c(8)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        data c/1.625d0,4.1328125d0,
     &       1.45380859375d+1,6.553353881835d+1,
     &       3.6066157150269d+2,2.3448727161884d+3,
     &       1.7588273098916d+4,1.4950639538279d+5/
        if (x.eq.0.0d0) then
           tti=0.0d0
           ttk=1.0d+300
           return
        endif
        if (x.lt.40.0d0) then
           tti=1.0d0
           r=1.0d0
           do 10 k=2,50
              r=.25d0*r*(k-1.0d0)/(k*k*k)*x*x
              tti=tti+r
              if (dabs(r/tti).lt.1.0d-12) go to 15
10         continue
15         tti=tti*.125d0*x*x
        else
           tti=1.0d0
           r=1.0d0
           do 20 k=1,8
              r=r/x
20            tti=tti+c(k)*r
           rc=x*dsqrt(2.0d0*pi*x)
           tti=tti*dexp(x)/rc
        endif
        if (x.le.12.0d0) then
           e0=(.5d0*dlog(x/2.0d0)+el)*dlog(x/2.0d0)
     &        +pi*pi/24.0d0+.5d0*el*el
           b1=1.5d0-(el+dlog(x/2.0d0))
           rs=1.0d0
           r=1.0d0
           do 25 k=2,50
              r=.25d0*r*(k-1.0d0)/(k*k*k)*x*x
              rs=rs+1.0d0/k
              r2=r*(rs+1.0d0/(2.0d0*k)-(el+dlog(x/2.0d0)))
              b1=b1+r2
              if (dabs(r2/b1).lt.1.0d-12) go to 30
25         continue
30         ttk=e0-.125d0*x*x*b1
        else
           ttk=1.0d0
           r=1.0d0
           do 35 k=1,8
              r=-r/x
35            ttk=ttk+c(k)*r
           rc=x*dsqrt(2.0d0/pi*x)
           ttk=ttk*dexp(-x)/rc
        endif
        return
        end
