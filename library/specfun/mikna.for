        program mikna
c
c       =============================================================
c       purpose: this program computes modified bessel functions 
c                in(x) and kn(x), and their derivatives using
c                subroutine ikna
c       input:   x --- argument of in(x) and kn(x) ( x ע 0 )
c                n --- order of in(x) and kn(x)
c                      ( n = 0,1,תתת, n ף 250 )
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
        call ikna(n,x,nm,bi,di,bk,dk)
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


        subroutine ikna(n,x,nm,bi,di,bk,dk)
c
c       ========================================================
c       purpose: compute modified bessel functions in(x) and
c                kn(x), and their derivatives
c       input:   x --- argument of in(x) and kn(x) ( x ע 0 )
c                n --- order of in(x) and kn(x)
c       output:  bi(n) --- in(x)
c                di(n) --- in'(x)
c                bk(n) --- kn(x)
c                dk(n) --- kn'(x)
c                nm --- highest order computed
c       routines called:
c            (1) ik01a for computing i0(x),i1(x),k0(x) & k1(x)
c            (2) msta1 and msta2 for computing the starting 
c                point for backward recurrence
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension bi(0:n),di(0:n),bk(0:n),dk(0:n)
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
        call ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
        bi(0)=bi0
        bi(1)=bi1
        bk(0)=bk0
        bk(1)=bk1
        di(0)=di0
        di(1)=di1
        dk(0)=dk0
        dk(1)=dk1
        if (n.le.1) return
        if (x.gt.40.0.and.n.lt.int(0.25*x)) then
           h0=bi0
           h1=bi1
           do 15 k=2,n
           h=-2.0d0*(k-1.0d0)/x*h1+h0
           bi(k)=h
           h0=h1
15         h1=h
        else
           m=msta1(x,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(x,n,15)
           endif
           f0=0.0d0
           f1=1.0d-100
           do 20 k=m,0,-1
              f=2.0d0*(k+1.0d0)*f1/x+f0
              if (k.le.nm) bi(k)=f
              f0=f1
20            f1=f
           s0=bi0/f
           do 25 k=0,nm
25            bi(k)=s0*bi(k)
        endif
        g0=bk0
        g1=bk1
        do 30 k=2,nm
           g=2.0d0*(k-1.0d0)/x*g1+g0
           bk(k)=g
           g0=g1
30         g1=g
        do 40 k=2,nm
           di(k)=bi(k-1)-k/x*bi(k)
40         dk(k)=-bk(k-1)-k/x*bk(k)
        return
        end


        subroutine ik01a(x,bi0,di0,bi1,di1,bk0,dk0,bk1,dk1)
c
c       =========================================================
c       purpose: compute modified bessel functions i0(x), i1(1),
c                k0(x) and k1(x), and their derivatives
c       input :  x   --- argument ( x ע 0 )
c       output:  bi0 --- i0(x)
c                di0 --- i0'(x)
c                bi1 --- i1(x)
c                di1 --- i1'(x)
c                bk0 --- k0(x)
c                dk0 --- k0'(x)
c                bk1 --- k1(x)
c                dk1 --- k1'(x)
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension a(12),b(12),a1(8)
        pi=3.141592653589793d0
        el=0.5772156649015329d0
        x2=x*x
        if (x.eq.0.0d0) then
           bi0=1.0d0
           bi1=0.0d0
           bk0=1.0d+300
           bk1=1.0d+300
           di0=0.0d0
           di1=0.5d0
           dk0=-1.0d+300
           dk1=-1.0d+300
           return
        else if (x.le.18.0) then
           bi0=1.0d0
           r=1.0d0
           do 15 k=1,50
              r=0.25d0*r*x2/(k*k)
              bi0=bi0+r
              if (dabs(r/bi0).lt.1.0d-15) go to 20
15         continue
20         bi1=1.0d0
           r=1.0d0
           do 25 k=1,50
              r=0.25d0*r*x2/(k*(k+1))
              bi1=bi1+r
              if (dabs(r/bi1).lt.1.0d-15) go to 30
25         continue
30         bi1=0.5d0*x*bi1
        else
           data a/0.125d0,7.03125d-2,
     &            7.32421875d-2,1.1215209960938d-1,
     &            2.2710800170898d-1,5.7250142097473d-1,
     &            1.7277275025845d0,6.0740420012735d0,
     &            2.4380529699556d01,1.1001714026925d02,
     &            5.5133589612202d02,3.0380905109224d03/
           data b/-0.375d0,-1.171875d-1,
     &            -1.025390625d-1,-1.4419555664063d-1,
     &            -2.7757644653320d-1,-6.7659258842468d-1,
     &            -1.9935317337513d0,-6.8839142681099d0,
     &            -2.7248827311269d01,-1.2159789187654d02,
     &            -6.0384407670507d02,-3.3022722944809d03/
           k0=12
           if (x.ge.35.0) k0=9
           if (x.ge.50.0) k0=7
           ca=dexp(x)/dsqrt(2.0d0*pi*x)
           bi0=1.0d0
           xr=1.0d0/x
           do 35 k=1,k0
35            bi0=bi0+a(k)*xr**k
           bi0=ca*bi0
           bi1=1.0d0
           do 40 k=1,k0
40            bi1=bi1+b(k)*xr**k
           bi1=ca*bi1
        endif
        if (x.le.9.0d0) then
           ct=-(dlog(x/2.0d0)+el)
           bk0=0.0d0
           w0=0.0d0
           r=1.0d0
           do 65 k=1,50
              w0=w0+1.0d0/k
              r=0.25d0*r/(k*k)*x2
              bk0=bk0+r*(w0+ct)
              if (dabs((bk0-ww)/bk0).lt.1.0d-15) go to 70
65            ww=bk0
70         bk0=bk0+ct
        else
           data a1/0.125d0,0.2109375d0,
     &             1.0986328125d0,1.1775970458984d01,
     &             2.1461706161499d02,5.9511522710323d03,
     &             2.3347645606175d05,1.2312234987631d07/
           cb=0.5d0/x
           xr2=1.0d0/x2
           bk0=1.0d0
           do 75 k=1,8
75            bk0=bk0+a1(k)*xr2**k
           bk0=cb*bk0/bi0
        endif
        bk1=(1.0d0/x-bi1*bk0)/bi0
        di0=bi1
        di1=bi0-bi1/x
        dk0=-bk1
        dk1=-bk0-bk1/x
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
c       input :  x  --- argument of jn(x)
c                n  --- order of jn(x)
c                mp --- significant digit
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
