        program mitikb
c
c       ============================================================
c       purpose: this program evaluates the integral of modified 
c                bessel functions i0(t) and k0(t) with respect to t
c                from 0 to x using subroutine itikb
c       input :  x  --- upper limit of the integral  ( x ò 0 )
c       output:  ti --- integration of i0(t) from 0 to x
c                tk --- integration of k0(t) from 0 to x
c       example:
c                    x         i0(t)dt         k0(t)dt
c                 -------------------------------------
c                   5.0     .318487d+02       1.567387
c                  10.0     .299305d+04       1.570779
c                  15.0     .352619d+06       1.570796
c                  20.0     .447586d+08       1.570796
c                  25.0     .589919d+10       1.570796
c       ============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x         i0(t)dt         k0(t)dt'
        write(*,*)'--------------------------------------'
        call itikb(x,ti,tk)
        write(*,10)x,ti,tk
10      format(1x,f5.1,d16.6,f15.6)
        end


        subroutine itikb(x,ti,tk)
c
c       =======================================================
c       purpose: integrate bessel functions i0(t) and k0(t)
c                with respect to t from 0 to x
c       input :  x  --- upper limit of the integral ( x ò 0 )
c       output:  ti --- integration of i0(t) from 0 to x
c                tk --- integration of k0(t) from 0 to x
c       =======================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        if (x.eq.0.0d0) then
           ti=0.0d0
        else if (x.lt.5.0d0) then
           t1=x/5.0d0
           t=t1*t1
           ti=((((((((.59434d-3*t+.4500642d-2)*t
     &        +.044686921d0)*t+.300704878d0)*t+1.471860153d0)
     &        *t+4.844024624d0)*t+9.765629849d0)*t
     &        +10.416666367d0)*t+5.0d0)*t1
        else if (x.ge.5.0.and.x.le.8.0d0) then
           t=5.0d0/x
           ti=(((-.015166d0*t-.0202292d0)*t+.1294122d0)*t
     &        -.0302912d0)*t+.4161224d0
           ti=ti*dexp(x)/dsqrt(x)
        else
           t=8.0d0/x
           ti=(((((-.0073995d0*t+.017744d0)*t-.0114858d0)*t
     &        +.55956d-2)*t+.59191d-2)*t+.0311734d0)*t
     &        +.3989423d0
           ti=ti*dexp(x)/dsqrt(x)
        endif
        if (x.eq.0.0d0) then
           tk=0.0d0
        else if (x.le.2.0d0) then
           t1=x/2.0d0
           t=t1*t1
           tk=((((((.116d-5*t+.2069d-4)*t+.62664d-3)*t
     &        +.01110118d0)*t+.11227902d0)*t+.50407836d0)*t
     &        +.84556868d0)*t1
              tk=tk-dlog(x/2.0d0)*ti
        else if (x.gt.2.0.and.x.le.4.0d0) then
           t=2.0d0/x
           tk=(((.0160395d0*t-.0781715d0)*t+.185984d0)*t
     &        -.3584641d0)*t+1.2494934d0
           tk=pi/2.0d0-tk*dexp(-x)/dsqrt(x)
        else if (x.gt.4.0.and.x.le.7.0d0) then
           t=4.0d0/x
           tk=(((((.37128d-2*t-.0158449d0)*t+.0320504d0)*t
     &        -.0481455d0)*t+.0787284d0)*t-.1958273d0)*t
     &        +1.2533141d0
           tk=pi/2.0d0-tk*dexp(-x)/dsqrt(x)
        else
           t=7.0d0/x
           tk=(((((.33934d-3*t-.163271d-2)*t+.417454d-2)*t
     &        -.933944d-2)*t+.02576646d0)*t-.11190289d0)*t
     &        +1.25331414d0
           tk=pi/2.0d0-tk*dexp(-x)/dsqrt(x)
        endif
        return
        end
