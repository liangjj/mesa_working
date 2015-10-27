        program mfcs
c
c       =======================================================
c       purpose: this program computes the fresnel integrals 
c                c(x) and s(x) using subroutine fcs
c       input :  x --- argument of c(x) and s(x)
c       output:  c --- c(x)
c                s --- s(x)
c       example:
c                  x          c(x)          s(x)
c                -----------------------------------
c                 0.0      .00000000      .00000000
c                 0.5      .49234423      .06473243
c                 1.0      .77989340      .43825915
c                 1.5      .44526118      .69750496
c                 2.0      .48825341      .34341568
c                 2.5      .45741301      .61918176
c       =======================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x '
        read(*,*) x
        write(*,*)'   x          c(x)          s(x)'
        write(*,*)' -----------------------------------'
        call fcs(x,c,s)
        write(*,10)x,c,s
10      format(1x,f5.1,2f15.8)
        end


        subroutine fcs(x,c,s)
c
c       =================================================
c       purpose: compute fresnel integrals c(x) and s(x)
c       input :  x --- argument of c(x) and s(x)
c       output:  c --- c(x)
c                s --- s(x)
c       =================================================
c
        implicit double precision (a-h,o-z)
        eps=1.0d-15
        pi=3.141592653589793d0
        xa=dabs(x)
        px=pi*xa
        t=.5d0*px*xa
        t2=t*t
        if (xa.eq.0.0) then
           c=0.0d0
           s=0.0d0
        else if (xa.lt.2.5d0) then
           r=xa
           c=r
           do 10 k=1,50
              r=-.5d0*r*(4.0d0*k-3.0d0)/k/(2.0d0*k-1.0d0)
     &          /(4.0d0*k+1.0d0)*t2
              c=c+r
              if (dabs(r).lt.dabs(c)*eps) go to 15
10         continue
15         s=xa*t/3.0d0
           r=s
           do 20 k=1,50
              r=-.5d0*r*(4.0d0*k-1.0d0)/k/(2.0d0*k+1.0d0)
     &          /(4.0d0*k+3.0d0)*t2
              s=s+r
              if (dabs(r).lt.dabs(s)*eps) go to 40
20         continue
        else if (xa.lt.4.5d0) then
           m=int(42.0+1.75*t)
           su=0.0d0
           c=0.0d0
           s=0.0d0
           f1=0.0d0
           f0=1.0d-100
           do 25 k=m,0,-1
              f=(2.0d0*k+3.0d0)*f0/t-f1
              if (k.eq.int(k/2)*2) then
                 c=c+f
              else
                 s=s+f
              endif
              su=su+(2.0d0*k+1.0d0)*f*f
              f1=f0
25            f0=f
           q=dsqrt(su)
           c=c*xa/q
           s=s*xa/q
        else
           r=1.0d0
           f=1.0d0
           do 30 k=1,20
              r=-.25d0*r*(4.0d0*k-1.0d0)*(4.0d0*k-3.0d0)/t2
30            f=f+r
           r=1.0d0/(px*xa)
           g=r
           do 35 k=1,12
              r=-.25d0*r*(4.0d0*k+1.0d0)*(4.0d0*k-1.0d0)/t2
35            g=g+r
           t0=t-int(t/(2.0d0*pi))*2.0d0*pi
           c=.5d0+(f*dsin(t0)-g*dcos(t0))/px
           s=.5d0-(f*dcos(t0)+g*dsin(t0))/px
        endif
40      if (x.lt.0.0d0) then
           c=-c
           s=-s
        endif
        return
        end
