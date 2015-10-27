        program mstvh0
c
c       ====================================================
c       purpose: this program computes struve function 
c                h0(x) using subroutine stvh0
c       input :  x   --- argument of h0(x) ( x ò 0 )
c       output:  sh0 --- h0(x)
c       example:
c                   x          h0(x)
c                ----------------------
c                  0.0       .00000000
c                  5.0      -.18521682
c                 10.0       .11874368
c                 15.0       .24772383
c                 20.0       .09439370
c                 25.0      -.10182519
c       ====================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x          h0(x)'
        write(*,*)'----------------------'
        call stvh0(x,sh0)
        write(*,10)x,sh0
10      format(1x,f5.1,e16.8)
        end


        subroutine stvh0(x,sh0)
c
c       =============================================
c       purpose: compute struve function h0(x)
c       input :  x   --- argument of h0(x) ( x ò 0 )
c       output:  sh0 --- h0(x)
c       =============================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        s=1.0d0
        r=1.0d0
        if (x.le.20.0d0) then
           a0=2.0*x/pi
           do 10 k=1,60
              r=-r*x/(2.0d0*k+1.0d0)*x/(2.0d0*k+1.0d0)
              s=s+r
              if (dabs(r).lt.dabs(s)*1.0d-12) go to 15
10         continue
15         sh0=a0*s
        else
           km=int(.5*(x+1.0))
           if (x.ge.50.0) km=25
           do 20 k=1,km
              r=-r*((2.0d0*k-1.0d0)/x)**2
              s=s+r
              if (dabs(r).lt.dabs(s)*1.0d-12) go to 25
20         continue
25         t=4.0d0/x
           t2=t*t
           p0=((((-.37043d-5*t2+.173565d-4)*t2-.487613d-4)
     &        *t2+.17343d-3)*t2-.1753062d-2)*t2+.3989422793d0
           q0=t*(((((.32312d-5*t2-.142078d-4)*t2+.342468d-4)*
     &        t2-.869791d-4)*t2+.4564324d-3)*t2-.0124669441d0)
           ta0=x-.25d0*pi
           by0=2.0d0/dsqrt(x)*(p0*dsin(ta0)+q0*dcos(ta0))
           sh0=2.0d0/(pi*x)*s+by0
        endif
        return
        end
