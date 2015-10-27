        program mlpn
c
c       ========================================================
c       purpose: this program computes the legendre polynomials 
c                pn(x) and their derivatives pn'(x) using
c                subroutine lpn
c       input :  x --- argument of pn(x)
c                n --- degree of pn(x) ( n = 0,1,...)
c       output:  pn(n) --- pn(x)
c                pd(n) --- pn'(x)
c       example:    x = 0.5
c                  n          pn(x)            pn'(x)
c                ---------------------------------------
c                  0       1.00000000        .00000000
c                  1        .50000000       1.00000000
c                  2       -.12500000       1.50000000
c                  3       -.43750000        .37500000
c                  4       -.28906250      -1.56250000
c                  5        .08984375      -2.22656250
c       ========================================================
c
        double precision pn,pd,x
        dimension pn(0:100),pd(0:100)
        write(*,*)'  please enter nmax and x '
        read(*,*)n,x
        write(*,30)x
        write(*,*)
        call lpn(n,x,pn,pd)
        write(*,*)'  n         pn(x)           pn''(x)'
        write(*,*)'---------------------------------------'
        do 10 k=0,n
10         write(*,20)k,pn(k),pd(k)
20      format(1x,i3,2e17.8)
30      format(3x,'x =',f5.1)
        end


        subroutine lpn(n,x,pn,pd)
c
c       ===============================================
c       purpose: compute legendre polynomials pn(x)
c                and their derivatives pn'(x)
c       input :  x --- argument of pn(x)
c                n --- degree of pn(x) ( n = 0,1,...)
c       output:  pn(n) --- pn(x)
c                pd(n) --- pn'(x)
c       ===============================================
c
        implicit double precision (p,x)
        dimension pn(0:n),pd(0:n)
        pn(0)=1.0d0
        pn(1)=x
        pd(0)=0.0d0
        pd(1)=1.0d0
        p0=1.0d0
        p1=x
        do 10 k=2,n
           pf=(2.0d0*k-1.0d0)/k*x*p1-(k-1.0d0)/k*p0
           pn(k)=pf
           if (dabs(x).eq.1.0d0) then
              pd(k)=0.5d0*x**(k+1)*k*(k+1.0d0)
           else
              pd(k)=k*(p1-x*pf)/(1.0d0-x*x)
           endif
           p0=p1
10         p1=pf
        return
        end
