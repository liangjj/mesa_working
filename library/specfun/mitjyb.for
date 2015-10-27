        program mitjyb
c
c       ===========================================================
c       purpose: this program evaluates the integral of bessel
c                functions j0(t) and y0(t) with respect to t 
c                from 0 to x using subroutine itjyb
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
        write(*,*)'pleas enter x '
        read(*,*)x
        write(*,*)'   x         j0(t)dt          y0(t)dt'
        write(*,*)'---------------------------------------'
        call itjyb(x,tj,ty)
        write(*,10)x,tj,ty
10      format(1x,f5.1,2f16.8)
        end


        subroutine itjyb(x,tj,ty)
c
c       =======================================================
c       purpose: integrate bessel functions j0(t) and y0(t)
c                with respect to t from 0 to x ( x ò 0 )
c       input :  x  --- upper limit of the integral
c       output:  tj --- integration of j0(t) from 0 to x
c                ty --- integration of y0(t) from 0 to x
c       =======================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        if (x.eq.0.0d0) then
           tj=0.0d0
           ty=0.0d0
        else if (x.le.4.0d0) then
           x1=x/4.0d0
           t=x1*x1
           tj=(((((((-.133718d-3*t+.2362211d-2)*t
     &        -.025791036d0)*t+.197492634d0)*t-1.015860606d0)
     &        *t+3.199997842d0)*t-5.333333161d0)*t+4.0d0)*x1
           ty=((((((((.13351d-4*t-.235002d-3)*t+.3034322d-2)*
     &        t-.029600855d0)*t+.203380298d0)*t-.904755062d0)
     &        *t+2.287317974d0)*t-2.567250468d0)*t
     &        +1.076611469d0)*x1
           ty=2.0d0/pi*dlog(x/2.0d0)*tj-ty
        else if (x.le.8.0d0) then
           xt=x-.25d0*pi
           t=16.0d0/(x*x)
           f0=((((((.1496119d-2*t-.739083d-2)*t+.016236617d0)
     &        *t-.022007499d0)*t+.023644978d0)
     &        *t-.031280848d0)*t+.124611058d0)*4.0d0/x
           g0=(((((.1076103d-2*t-.5434851d-2)*t+.01242264d0)
     &        *t-.018255209)*t+.023664841d0)*t-.049635633d0)
     &        *t+.79784879d0
           tj=1.0d0-(f0*dcos(xt)-g0*dsin(xt))/dsqrt(x)
           ty=-(f0*dsin(xt)+g0*dcos(xt))/dsqrt(x)
        else
           t=64.0d0/(x*x)
           xt=x-.25d0*pi
           f0=(((((((-.268482d-4*t+.1270039d-3)*t
     &        -.2755037d-3)*t+.3992825d-3)*t-.5366169d-3)*t
     &        +.10089872d-2)*t-.40403539d-2)*t+.0623347304d0)
     &        *8.0d0/x
           g0=((((((-.226238d-4*t+.1107299d-3)*t-.2543955d-3)
     &        *t+.4100676d-3)*t-.6740148d-3)*t+.17870944d-2)
     &        *t-.01256424405d0)*t+.79788456d0
           tj=1.0d0-(f0*dcos(xt)-g0*dsin(xt))/dsqrt(x)
           ty=-(f0*dsin(xt)+g0*dcos(xt))/dsqrt(x)
        endif
        return
        end

