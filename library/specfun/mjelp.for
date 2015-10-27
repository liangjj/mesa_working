        program mjelp
c
c       ============================================================
c       purpose: this program computes jacobian elliptic functions 
c                sn u, cn u and dn u using subroutine jelp
c       input  : u   --- argument of jacobian elliptic fuctions
c                hk  --- modulus k ( 0 ó k ó 1 )
c       output : esn --- sn u
c                ecn --- cn u
c                edn --- dn u
c                eph --- phi ( in degrees )
c       example:
c                k = .5, ( k(k) = 1.68575035 ), and u = u0*k
c
c                u0       phi       sn u        cn u        dn u
c              ----------------------------------------------------
c               0.0      .0000    .0000000   1.0000000   1.0000000
c               0.5    47.0586    .7320508    .6812500    .9306049
c               1.0    90.0000   1.0000000    .0000000    .8660254
c               1.5   132.9414    .7320508   -.6812500    .9306049
c               2.0   180.0000    .0000000  -1.0000000   1.0000000
c       ============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter k and u '
        read(*,*)hk,u
        write(*,*)
        write(*,*)'   k        u          phi        sn u',
     &            '        cn u        dn u'
        write(*,*)' -------------------------------------',
     &            '---------------------------'
        call jelp(u,hk,esn,ecn,edn,eph)
        write(*,10)hk,u,eph,esn,ecn,edn
10      format(1x,f5.3,f12.7,2x,f9.5,3f12.7)
        end


        subroutine jelp(u,hk,esn,ecn,edn,eph)
c
c       ========================================================
c       purpose: compute jacobian elliptic functions sn u, cn u
c                and dn u
c       input  : u   --- argument of jacobian elliptic fuctions
c                hk  --- modulus k ( 0 ó k ó 1 )
c       output : esn --- sn u
c                ecn --- cn u
c                edn --- dn u
c                eph --- phi ( in degrees )
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension r(40)
        pi=3.14159265358979d0
        a0=1.0d0
        b0=dsqrt(1.0d0-hk*hk)
        do 10 n=1,40
           a=(a0+b0)/2.0d0
           b=dsqrt(a0*b0)
           c=(a0-b0)/2.0d0
           r(n)=c/a
           if (c.lt.1.0d-7) go to 15
           a0=a
10         b0=b
15      dn=2.0d0**n*a*u
        do 20 j=n,1,-1
           t=r(j)*dsin(dn)
           sa=datan(t/dsqrt(dabs(1.0d0-t*t)))
           d=.5d0*(dn+sa)
20         dn=d
        eph=d*180.0d0/pi
        esn=dsin(d)
        ecn=dcos(d)
        edn=dsqrt(1.0d0-hk*hk*esn*esn)
        return
        end
