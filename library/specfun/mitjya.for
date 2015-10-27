        program mitjya
c
c       ===========================================================
c       purpose: this program evaluates the integral of bessel
c                functions j0(t) and y0(t) with respect to t 
c                from 0 to x using subroutine itjya
c       input :  x  --- upper limit of the integral ( x ò 0 )
c       output:  tj --- integration of j0(t) from 0 to x
c                ty --- integration of y0(t) from 0 to x
c       example:
c                   x         j0(t)dt          y0(t)dt
c                ---------------------------------------
c                  5.0       .71531192       .19971938
c                 10.0      1.06701130       .24129032
c                 15.0      1.20516194       .00745772
c                 20.0      1.05837882      -.16821598
c                 25.0       .87101492      -.09360793
c                 30.0       .88424909       .08822971
c       ===========================================================
c
        double precision x,tj,ty
        write(*,*)'please enter x '
        read(*,*)x
        write(*,*)'   x         j0(t)dt          y0(t)dt'
        write(*,*)'---------------------------------------'
        call itjya(x,tj,ty)
        write(*,10)x,tj,ty
10      format(1x,f5.1,2f16.8)
        end


        subroutine itjya(x,tj,ty)
c
c       ==========================================================
c       purpose: integrate bessel functions j0(t) & y0(t) with
c                respect to t from 0 to x
c       input :  x  --- upper limit of the integral ( x ò 0 )
c       output:  tj --- integration of j0(t) from 0 to x
c                ty --- integration of y0(t) from 0 to x
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension a(18)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        eps=1.0d-12
        if (x.eq.0.0d0) then
           tj=0.0d0
           ty=0.0d0
        else if (x.le.20.0d0) then
           x2=x*x
           tj=x
           r=x
           do 10 k=1,60
              r=-.25d0*r*(2*k-1.0d0)/(2*k+1.0d0)/(k*k)*x2
              tj=tj+r
              if (dabs(r).lt.dabs(tj)*eps) go to 15
10         continue
15         ty1=(el+dlog(x/2.0d0))*tj
           rs=0.0d0
           ty2=1.0d0
           r=1.0d0
           do 20 k=1,60
              r=-.25d0*r*(2*k-1.0d0)/(2*k+1.0d0)/(k*k)*x2
              rs=rs+1.0d0/k
              r2=r*(rs+1.0d0/(2.0d0*k+1.0d0))
              ty2=ty2+r2
              if (dabs(r2).lt.dabs(ty2)*eps) go to 25
20         continue
25         ty=(ty1-x*ty2)*2.0d0/pi
        else
           a0=1.0d0
           a1=5.0d0/8.0d0
           a(1)=a1
           do 30 k=1,16
              af=((1.5d0*(k+.5d0)*(k+5.0d0/6.0d0)*a1-.5d0
     &           *(k+.5d0)*(k+.5d0)*(k-.5d0)*a0))/(k+1.0d0)
              a(k+1)=af
              a0=a1
30            a1=af
           bf=1.0d0
           r=1.0d0
           do 35 k=1,8
              r=-r/(x*x)
35            bf=bf+a(2*k)*r
           bg=a(1)/x
           r=1.0d0/x
           do 40 k=1,8
              r=-r/(x*x)
40            bg=bg+a(2*k+1)*r
           xp=x+.25d0*pi
           rc=dsqrt(2.0d0/(pi*x))
           tj=1.0d0-rc*(bf*dcos(xp)+bg*dsin(xp))
           ty=rc*(bg*dcos(xp)-bf*dsin(xp))
        endif
        return
        end
