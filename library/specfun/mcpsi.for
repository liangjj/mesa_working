        program mcpsi
c
c       =========================================================
c       purpose: this program computes the psi function psi(z)
c                for a complex argument using subroutine cpsi
c       input :  x   --- real part of z
c                y   --- imaginary part of z
c       output:  psr --- real part of psi(z)
c                psi --- imaginary part of psi(z)
c       examples:
c                   x       y      re[psi(z)]     im[psi(z)]
c                 -------------------------------------------
c                  3.0     2.0     1.16459152      .67080728
c                  3.0    -2.0     1.16459152     -.67080728
c                 -3.0     2.0     1.39536075     2.62465344
c                 -3.0    -2.0     1.39536075    -2.62465344
c       =========================================================
c
        double precision x,y,psr,psi
        write(*,*)'please enter x and y ( z=x+iy)'
        read(*,*)x,y
        write(*,20)x,y
        write(*,*)
        write(*,*)'   x       y      re[psi(z)]      im[psi(z)]'
        write(*,*)' ----------------------------------------------'
        call cpsi(x,y,psr,psi)
        write(*,10)x,y,psr,psi
10      format(1x,f5.1,3x,f5.1,2x,2e16.8)
20      format(1x,2hx=,f6.2,6x,2hy=,f6.2)
        end


        subroutine cpsi(x,y,psr,psi)
c
c       =============================================
c       purpose: compute the psi function for a
c                complex argument
c       input :  x   --- real part of z
c                y   --- imaginary part of z
c       output:  psr --- real part of psi(z)
c                psi --- imaginary part of psi(z)
c       =============================================
c
        implicit double precision (a-h,o-z)
        dimension a(8)
        data a/-.8333333333333d-01,.83333333333333333d-02,
     &       -.39682539682539683d-02,.41666666666666667d-02,
     &       -.75757575757575758d-02,.21092796092796093d-01,
     &       -.83333333333333333d-01,.4432598039215686d0/
        pi=3.141592653589793d0
        if (y.eq.0.0d0.and.x.eq.int(x).and.x.le.0.0d0) then
           psr=1.0d+300
           psi=0.0d0
        else
           if (x.lt.0.0d0) then
              x1=x
              y1=y
              x=-x
              y=-y
           endif
           x0=x
           if (x.lt.8.0d0) then
              n=8-int(x)
              x0=x+n
           endif
           if (x0.eq.0.0d0.and.y.ne.0.0d0) th=0.5d0*pi
           if (x0.ne.0.0d0) th=datan(y/x0)
           z2=x0*x0+y*y
           z0=dsqrt(z2)
           psr=dlog(z0)-0.5d0*x0/z2
           psi=th+0.5d0*y/z2
           do 10 k=1,8
              psr=psr+a(k)*z2**(-k)*dcos(2.0d0*k*th)
10            psi=psi-a(k)*z2**(-k)*dsin(2.0d0*k*th)
           if (x.lt.8.0d0) then
              rr=0.0d0
              ri=0.0d0
              do 20 k=1,n
                 rr=rr+(x0-k)/((x0-k)**2.0d0+y*y)
20               ri=ri+y/((x0-k)**2.0d0+y*y)
              psr=psr-rr
              psi=psi+ri
           endif
           if (x1.lt.0.0d0) then
              tn=dtan(pi*x)
              tm=dtanh(pi*y)
              ct2=tn*tn+tm*tm
              psr=psr+x/(x*x+y*y)+pi*(tn-tn*tm*tm)/ct2
              psi=psi-y/(x*x+y*y)-pi*tm*(1.0d0+tn*tn)/ct2
              x=x1
              y=y1
           endif
        endif
        return
        end
