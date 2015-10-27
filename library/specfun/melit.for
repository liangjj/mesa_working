        program melit
c
c       ==========================================================
c       purpose: this program computes complete and incomplete 
c                elliptic integrals f(k,phi) and e(k,phi) using
c                subroutine elit
c       input  : hk  --- modulus k ( 0 ó k ó 1 )
c                phi --- argument ( in degrees )
c       output : fe  --- f(k,phi)
c                ee  --- e(k,phi)
c       example:
c                k = .5
c
c                 phi     f(k,phi)       e(k,phi)
c                -----------------------------------
c                   0      .00000000      .00000000
c                  15      .26254249      .26106005
c                  30      .52942863      .51788193
c                  45      .80436610      .76719599
c                  60     1.08955067     1.00755556
c                  75     1.38457455     1.23988858
c                  90     1.68575035     1.46746221
c       ==========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter k and phi (in degs.) '
        read(*,*)hk,phi
        write(*,*)
        write(*,*)'  phi     f(k,phi)       e(k,phi)'
        write(*,*)' -----------------------------------'
        call elit(hk,phi,fe,ee)
        write(*,10)phi,fe,ee
10      format(1x,i4,2f15.8)
        end


        subroutine elit(hk,phi,fe,ee)
c
c       ==================================================
c       purpose: compute complete and incomplete elliptic
c                integrals f(k,phi) and e(k,phi)
c       input  : hk  --- modulus k ( 0 ó k ó 1 )
c                phi --- argument ( in degrees )
c       output : fe  --- f(k,phi)
c                ee  --- e(k,phi)
c       ==================================================
c
        implicit double precision (a-h,o-z)
        g=0.0d0
        pi=3.14159265358979d0
        a0=1.0d0
        b0=dsqrt(1.0d0-hk*hk)
        d0=(pi/180.0d0)*phi
        r=hk*hk
        if (hk.eq.1.0d0.and.phi.eq.90.0d0) then
           fe=1.0d+300
           ee=1.0d0
        else if (hk.eq.1.0d0) then
           fe=dlog((1.0d0+dsin(d0))/dcos(d0))
           ee=dsin(d0)
        else
           fac=1.0d0
           do 10 n=1,40
              a=(a0+b0)/2.0d0
              b=dsqrt(a0*b0)
              c=(a0-b0)/2.0d0
              fac=2.0d0*fac
              r=r+fac*c*c
              if (phi.ne.90.0d0) then
                 d=d0+datan((b0/a0)*dtan(d0))
                 g=g+c*dsin(d)
                 d0=d+pi*int(d/pi+.5d0)
              endif
              a0=a
              b0=b
              if (c.lt.1.0d-7) go to 15
10         continue
15         ck=pi/(2.0d0*a)
           ce=pi*(2.0d0-r)/(4.0d0*a)
           if (phi.eq.90.0d0) then
              fe=ck
              ee=ce
           else
              fe=d/(fac*a)
              ee=fe*ce/ck+g
           endif
        endif
        return
        end
