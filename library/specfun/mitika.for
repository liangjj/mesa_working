        program mitika
c
c       ============================================================
c       purpose: this program evaluates the integral of modified 
c                bessel functions i0(t) and k0(t) with respect to t
c                from 0 to x using subroutine itika
c       input :  x  --- upper limit of the integral  ( x ò 0 )
c       output:  ti --- integration of i0(t) from 0 to x
c                tk --- integration of k0(t) from 0 to x
c       example:
c                    x         i0(t)dt         k0(t)dt
c                 --------------------------------------
c                   5.0    .31848668d+02     1.56738739
c                  10.0    .29930445d+04     1.57077931
c                  15.0    .35262048d+06     1.57079623
c                  20.0    .44758593d+08     1.57079633
c                  25.0    .58991731d+10     1.57079633
c       ============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x         i0(t)dt         k0(t)dt'
        write(*,*)' --------------------------------------'
        call itika(x,ti,tk)
        write(*,10)x,ti,tk
10      format(1x,f5.1,d17.8,f15.8)
        end


        subroutine itika(x,ti,tk)
c
c       =======================================================
c       purpose: integrate modified bessel functions i0(t) and
c                k0(t) with respect to t from 0 to x
c       input :  x  --- upper limit of the integral  ( x ò 0 )
c       output:  ti --- integration of i0(t) from 0 to x
c                tk --- integration of k0(t) from 0 to x
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension a(10)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        data a/.625d0,1.0078125d0,
     &       2.5927734375d0,9.1868591308594d0,
     &       4.1567974090576d+1,2.2919635891914d+2,
     &       1.491504060477d+3,1.1192354495579d+4,
     &       9.515939374212d+4,9.0412425769041d+5/
        if (x.eq.0.0d0) then
           ti=0.0d0
           tk=0.0d0
           return
        else if (x.lt.20.0d0) then
           x2=x*x
           ti=1.0d0
           r=1.0d0
           do 10 k=1,50
              r=.25d0*r*(2*k-1.0d0)/(2*k+1.0d0)/(k*k)*x2
              ti=ti+r
              if (dabs(r/ti).lt.1.0d-12) go to 15
10         continue
15         ti=ti*x
        else
           ti=1.0d0
           r=1.0d0
           do 20 k=1,10
              r=r/x
20            ti=ti+a(k)*r
           rc1=1.0d0/dsqrt(2.0d0*pi*x)
           ti=rc1*dexp(x)*ti
        endif
        if (x.lt.12.0d0) then
           e0=el+dlog(x/2.0d0)
           b1=1.0d0-e0
           b2=0.0d0
           rs=0.0d0
           r=1.0d0
           do 25 k=1,50
              r=.25d0*r*(2*k-1.0d0)/(2*k+1.0d0)/(k*k)*x2
              b1=b1+r*(1.0d0/(2*k+1)-e0)
              rs=rs+1.0d0/k
              b2=b2+r*rs
              tk=b1+b2
              if (dabs((tk-tw)/tk).lt.1.0d-12) go to 30
25            tw=tk
30         tk=tk*x
        else
           tk=1.0d0
           r=1.0d0
           do 35 k=1,10
              r=-r/x
35            tk=tk+a(k)*r
           rc2=dsqrt(pi/(2.0d0*x))
           tk=pi/2.0d0-rc2*tk*dexp(-x)
        endif
        return
        end
