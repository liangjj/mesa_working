        program mlpni
c
c       ========================================================
c       purpose: this program computes the legendre polynomials 
c                pn(x), pn'(x) and the integral of pn(t) from 0 
c                to x using subroutine lpni
c       input :  x --- argument of pn(x)
c                n --- degree of pn(x) ( n = 0,1,... )
c       output:  pn(n) --- pn(x)
c                pd(n) --- pn'(x)
c                pl(n) --- integral of pn(t) from 0 to x
c       example: x = 0.50
c                n       pn(x)         pn'(x)        pn(t)dt
c               ---------------------------------------------
c                0    1.00000000     .00000000     .50000000
c                1     .50000000    1.00000000     .12500000
c                2    -.12500000    1.50000000    -.18750000
c                3    -.43750000     .37500000    -.14843750
c                4    -.28906250   -1.56250000     .05859375
c                5     .08984375   -2.22656250     .11816406
c       ========================================================
c
        double precision pn,pd,pl,x
        dimension pn(0:100),pd(0:100),pl(0:100)
        write(*,*)'  please enter nmax and x'
        read(*,*)n,x
        write(*,30)x
        write(*,*)
        write(*,*)'  n        pn(x)          pn''(x)         pn(t)dt'
        write(*,*)' ---------------------------------------------------'
        call lpni(n,x,pn,pd,pl)
        do 10 k=0,n
10         write(*,20)k,pn(k),pd(k),pl(k)
20      format(1x,i3,3e16.8)
30      format(3x,'x =',f5.2)
        end


        subroutine lpni(n,x,pn,pd,pl)
c
c       =====================================================
c       purpose: compute legendre polynomials pn(x), pn'(x)
c                and the integral of pn(t) from 0 to x
c       input :  x --- argument of pn(x)
c                n --- degree of pn(x) ( n = 0,1,... )
c       output:  pn(n) --- pn(x)
c                pd(n) --- pn'(x)
c                pl(n) --- integral of pn(t) from 0 to x
c       =====================================================
c
        implicit double precision (p,r,x)
        dimension pn(0:n),pd(0:n),pl(0:n)
        pn(0)=1.0d0
        pn(1)=x
        pd(0)=0.0d0
        pd(1)=1.0d0
        pl(0)=x
        pl(1)=0.5d0*x*x
        p0=1.0d0
        p1=x
        do 15 k=2,n
           pf=(2.0d0*k-1.0d0)/k*x*p1-(k-1.0d0)/k*p0
           pn(k)=pf
           if (dabs(x).eq.1.0d0) then
              pd(k)=0.5d0*x**(k+1)*k*(k+1.0d0)
           else
              pd(k)=k*(p1-x*pf)/(1.0d0-x*x)
           endif
           pl(k)=(x*pn(k)-pn(k-1))/(k+1.0d0)
           p0=p1
           p1=pf
           if (k.eq.2*int(k/2)) go to 15
           r=1.0d0/(k+1.0d0)
           n1=(k-1)/2
           do 10 j=1,n1
10            r=(0.5d0/j-1.0d0)*r
           pl(k)=pl(k)+r
15      continue
        return
        end
