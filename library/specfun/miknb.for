        program miknb
c
c       =============================================================
c       purpose: this program computes modified bessel functions 
c                in(x) and kn(x), and their derivatives using
c                subroutine iknb
c       input:   x --- argument of in(x) and kn(x) ( 0 ó x ó 700 )
c                n --- order of in(x) and kn(x)
c                      ( n = 0,1,..., n ó 250 )
c       output:  bi(n) --- in(x)
c                di(n) --- in'(x)
c                bk(n) --- kn(x)
c                dk(n) --- kn'(x)
c       example: nmax = 5,    x = 10.0
c
c     n      in(x)          in'(x)         kn(x)         kn'(x)
c    ---------------------------------------------------------------
c     0   .2815717d+04   .2670988d+04   .1778006d-04  -.1864877d-04
c     1   .2670988d+04   .2548618d+04   .1864877d-04  -.1964494d-04
c     2   .2281519d+04   .2214685d+04   .2150982d-04  -.2295074d-04
c     3   .1758381d+04   .1754005d+04   .2725270d-04  -.2968563d-04
c     4   .1226491d+04   .1267785d+04   .3786144d-04  -.4239728d-04
c     5   .7771883d+03   .8378964d+03   .5754185d-04  -.6663236d-04
c       =============================================================
c
        implicit double precision (a-h,o-z)
        dimension bi(0:250),di(0:250),bk(0:250),dk(0:250)
        write(*,*)'  please enter n, x '
        read(*,*)n,x
        write(*,25)n,x
        write(*,*)
        if (n.le.10) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call iknb(n,x,nm,bi,di,bk,dk)
        write(*,*)'  n      in(x)          in''(x) ',
     &            '        kn(x)         kn''(x) '
        write(*,*)' -------------------------------',
     &            '--------------------------------'
        do 10 k=0,nm,ns
           write(*,20)k,bi(k),di(k),bk(k),dk(k)
10      continue
20      format(1x,i3,4d15.7)
25      format(3x,6hnmax =,i3,',    ',3hx =,f5.1)
        end


        subroutine iknb(n,x,nm,bi,di,bk,dk)
c
c       ============================================================
c       purpose: compute modified bessel functions in(x) and kn(x),
c                and their derivatives
c       input:   x --- argument of in(x) and kn(x) ( 0 ó x ó 700 )
c                n --- order of in(x) and kn(x)
c       output:  bi(n) --- in(x)
c                di(n) --- in'(x)
c                bk(n) --- kn(x)
c                dk(n) --- kn'(x)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the starting point 
c                for backward recurrence
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        dimension bi(0:n),di(0:n),bk(0:n),dk(0:n)
        pi=3.141592653589793d0
        el=0.5772156649015329d0
        nm=n
        if (x.le.1.0d-100) then
           do 10 k=0,n
              bi(k)=0.0d0
              di(k)=0.0d0
              bk(k)=1.0d+300
10            dk(k)=-1.0d+300
           bi(0)=1.0d0
           di(1)=0.5d0
           return
        endif
        if (n.eq.0) nm=1
        m=msta1(x,200)
        if (m.lt.nm) then
           nm=m
        else
           m=msta2(x,nm,15)
        endif
        bs=0.0d0
        sk0=0.0d0
        f0=0.0d0
        f1=1.0d-100
        do 15 k=m,0,-1
           f=2.0d0*(k+1.0d0)/x*f1+f0
           if (k.le.nm) bi(k)=f
           if (k.ne.0.and.k.eq.2*int(k/2)) sk0=sk0+4.0d0*f/k
           bs=bs+2.0d0*f
           f0=f1
15         f1=f
        s0=dexp(x)/(bs-f)
        do 20 k=0,nm
20         bi(k)=s0*bi(k)
        if (x.le.8.0d0) then
           bk(0)=-(dlog(0.5d0*x)+el)*bi(0)+s0*sk0
           bk(1)=(1.0d0/x-bi(1)*bk(0))/bi(0)
        else
           a0=dsqrt(pi/(2.0d0*x))*dexp(-x)
           k0=16
           if (x.ge.25.0) k0=10
           if (x.ge.80.0) k0=8
           if (x.ge.200.0) k0=6
           do 30 l=0,1
              bkl=1.0d0
              vt=4.0d0*l
              r=1.0d0
              do 25 k=1,k0
                 r=0.125d0*r*(vt-(2.0*k-1.0)**2)/(k*x)
25               bkl=bkl+r
              bk(l)=a0*bkl
30         continue
        endif
        g0=bk(0)
        g1=bk(1)
        do 35 k=2,nm
           g=2.0d0*(k-1.0d0)/x*g1+g0
           bk(k)=g
           g0=g1
35         g1=g
        di(0)=bi(1)
        dk(0)=-bk(1)
        do 40 k=1,nm
           di(k)=bi(k-1)-k/x*bi(k)
40         dk(k)=-bk(k-1)-k/x*bk(k)
        return
        end


        integer function msta1(x,mp)
c
c       ===================================================
c       purpose: determine the starting point for backward  
c                recurrence such that the magnitude of    
c                jn(x) at that point is about 10^(-mp)
c       input :  x     --- argument of jn(x)
c                mp    --- value of magnitude
c       output:  msta1 --- starting point   
c       ===================================================
c
        implicit double precision (a-h,o-z)
        a0=dabs(x)
        n0=int(1.1*a0)+1
        f0=envj(n0,a0)-mp
        n1=n0+5
        f1=envj(n1,a0)-mp
        do 10 it=1,20             
           nn=n1-(n1-n0)/(1.0d0-f0/f1)                  
           f=envj(nn,a0)-mp
           if(abs(nn-n1).lt.1) go to 20
           n0=n1
           f0=f1
           n1=nn
 10        f1=f
 20     msta1=nn
        return
        end


        integer function msta2(x,n,mp)
c
c       ===================================================
c       purpose: determine the starting point for backward
c                recurrence such that all jn(x) has mp
c                significant digits
c       input :  x     --- argument of jn(x)
c                n     --- order of jn(x)
c                mp    --- significant digit
c       output:  msta2 --- starting point
c       ===================================================
c
        implicit double precision (a-h,o-z)
        a0=dabs(x)
        hmp=0.5d0*mp
        ejn=envj(n,a0)
        if (ejn.le.hmp) then
           obj=mp
           n0=int(1.1*a0)
        else
           obj=hmp+ejn
           n0=n
        endif
        f0=envj(n0,a0)-obj
        n1=n0+5
        f1=envj(n1,a0)-obj
        do 10 it=1,20
           nn=n1-(n1-n0)/(1.0d0-f0/f1)
           f=envj(nn,a0)-obj
           if (abs(nn-n1).lt.1) go to 20
           n0=n1
           f0=f1
           n1=nn
10         f1=f
20      msta2=nn+10
        return
        end

        real*8 function envj(n,x)
        double precision x
        envj=0.5d0*dlog10(6.28d0*n)-n*dlog10(1.36d0*x/n)
        return
        end
