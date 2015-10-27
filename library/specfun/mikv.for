        program mikv
c
c       ============================================================
c       purpose: this program computes modified bessel functions 
c                iv(x) and kv(x) with an arbitrary order, and
c                their derivatives using subroutine ikv
c       input :  x --- argument ( x ע 0 )
c                v --- order of iv(x) and kv(x)
c                      ( v = n+v0, 0 ף n ף 250 , 0 ף v0 < 1 )
c       output:  bi(n) --- in+v0(x)
c                di(n) --- in+v0'(x)
c                bk(n) --- kn+v0(x)
c                dk(n) --- kn+v0'(x)
c       example: v = n+v0,   v0 = .25,   x = 10.0
c
c     n       iv(x)          iv'(x)         kv(x)          kv'(x)
c    ----------------------------------------------------------------
c     0   .28064359d+04  .26631677d+04  .17833184d-04 -.18709581d-04
c     1   .25930068d+04  .24823101d+04  .19155411d-04 -.20227611d-04
c     2   .21581842d+04  .21074153d+04  .22622037d-04 -.24245369d-04
c     3   .16218239d+04  .16310915d+04  .29335327d-04 -.32156018d-04
c     4   .11039987d+04  .11526244d+04  .41690000d-04 -.47053577d-04
c     5   .68342498d+03  .74520058d+03  .64771827d-04 -.75695209d-04
c       =============================================================
c
        implicit double precision (a-h,o-z)
        common bi(0:250),di(0:250),bk(0:250),dk(0:250)
        write(*,*)'  please enter v, x '
        read(*,*)v,x
        n=int(v)
        v0=v-n
        write(*,30)v0,x
        if (n.le.10) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call ikv(v,x,vm,bi,di,bk,dk)
        nm=int(vm)
        write(*,*)
        write(*,*)'    n       iv(x)           iv''(x)          ',
     &            'kv(x)           kv''(x)'
        write(*,*)'  -------------------------------------',
     &            '--------------------------------'
        do 10 k=0,nm,ns
10         write(*,20) k,bi(k),di(k),bk(k),dk(k)
20      format(2x,i4,1x,4d16.8)
30      format(8x,'v = n+v0',',   ','v0 =',f7.5,',   ','x =',f5.1)
        end


        subroutine ikv(v,x,vm,bi,di,bk,dk)
c
c       =======================================================
c       purpose: compute modified bessel functions iv(x) and
c                kv(x), and their derivatives
c       input :  x --- argument ( x ע 0 )
c                v --- order of iv(x) and kv(x)
c                      ( v = n+v0, n = 0,1,2,..., 0 ף v0 < 1 )
c       output:  bi(n) --- in+v0(x)
c                di(n) --- in+v0'(x)
c                bk(n) --- kn+v0(x)
c                dk(n) --- kn+v0'(x)
c                vm --- highest order computed
c       routines called:
c            (1) gamma for computing the gamma function
c            (2) msta1 and msta2 to compute the starting
c                point for backward recurrence
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension bi(0:*),di(0:*),bk(0:*),dk(0:*)
        pi=3.141592653589793d0
        x2=x*x
        n=int(v)
        v0=v-n
        if (n.eq.0) n=1
        if (x.lt.1.0d-100) then
           do 10 k=0,n
              bi(k)=0.0d0
              di(k)=0.0d0
              bk(k)=-1.0d+300
10            dk(k)=1.0d+300
           if (v.eq.0.0) then
              bi(0)=1.0d0
              di(1)=0.5d0
           endif
           vm=v
           return
        endif
        piv=pi*v0
        vt=4.0d0*v0*v0
        if (v0.eq.0.0d0) then
           a1=1.0d0
        else
           v0p=1.0d0+v0
           call gamma(v0p,gap)
           a1=(0.5d0*x)**v0/gap
        endif
        k0=14
        if (x.ge.35.0) k0=10
        if (x.ge.50.0) k0=8
        if (x.le.18.0) then
           bi0=1.0d0
           r=1.0d0
           do 15 k=1,30
              r=0.25d0*r*x2/(k*(k+v0))
              bi0=bi0+r
              if (dabs(r/bi0).lt.1.0d-15) go to 20
15         continue
20         bi0=bi0*a1
        else
           ca=dexp(x)/dsqrt(2.0d0*pi*x)
           sum=1.0d0
           r=1.0d0
           do 25 k=1,k0
              r=-0.125d0*r*(vt-(2.0d0*k-1.0d0)**2.0)/(k*x)
25            sum=sum+r
           bi0=ca*sum
        endif
        m=msta1(x,200)
        if (m.lt.n) then
           n=m
        else
           m=msta2(x,n,15)
        endif
        f2=0.0d0
        f1=1.0d-100
        do 30 k=m,0,-1
           f=2.0d0*(v0+k+1.0d0)/x*f1+f2
           if (k.le.n) bi(k)=f
           f2=f1
30         f1=f
        cs=bi0/f
        do 35 k=0,n
35         bi(k)=cs*bi(k)
        di(0)=v0/x*bi(0)+bi(1)
        do 40 k=1,n
40         di(k)=-(k+v0)/x*bi(k)+bi(k-1)
        if (x.le.9.0d0) then
           if (v0.eq.0.0d0) then
              ct=-dlog(0.5d0*x)-0.5772156649015329d0
              cs=0.0d0
              w0=0.0d0
              r=1.0d0
              do 45 k=1,50
                 w0=w0+1.0d0/k
                 r=0.25d0*r/(k*k)*x2
                 cs=cs+r*(w0+ct)
                 wa=dabs(cs)
                 if (dabs((wa-ww)/wa).lt.1.0d-15) go to 50
45               ww=wa
50            bk0=ct+cs
           else
              v0n=1.0d0-v0
              call gamma(v0n,gan)
              a2=1.0d0/(gan*(0.5d0*x)**v0)
              a1=(0.5d0*x)**v0/gap
              sum=a2-a1
              r1=1.0d0
              r2=1.0d0
              do 55 k=1,120
                 r1=0.25d0*r1*x2/(k*(k-v0))
                 r2=0.25d0*r2*x2/(k*(k+v0))
                 sum=sum+a2*r1-a1*r2
                 wa=dabs(sum)
                 if (dabs((wa-ww)/wa).lt.1.0d-15) go to 60
55               ww=wa
60            bk0=0.5d0*pi*sum/dsin(piv)
           endif
        else
           cb=dexp(-x)*dsqrt(0.5d0*pi/x)
           sum=1.0d0
           r=1.0d0
           do 65 k=1,k0
              r=0.125d0*r*(vt-(2.0*k-1.0)**2.0)/(k*x)
65            sum=sum+r
           bk0=cb*sum
        endif
        bk1=(1.0d0/x-bi(1)*bk0)/bi(0)
        bk(0)=bk0
        bk(1)=bk1
        do 70 k=2,n
           bk2=2.0d0*(v0+k-1.0d0)/x*bk1+bk0
           bk(k)=bk2
           bk0=bk1
70         bk1=bk2
        dk(0)=v0/x*bk(0)-bk(1)
        do 80 k=1,n
80         dk(k)=-(k+v0)/x*bk(k)-bk(k-1)
        vm=n+v0
        return
        end


        subroutine gamma(x,ga)
c
c       ==================================================
c       purpose: compute gamma function ג(x)
c       input :  x  --- argument of ג(x)
c                       ( x is not equal to 0,-1,-2,תתת)
c       output:  ga --- ג(x)
c       ==================================================
c
        implicit double precision (a-h,o-z)
        dimension g(26)
        pi=3.141592653589793d0
        if (x.eq.int(x)) then
           if (x.gt.0.0d0) then
              ga=1.0d0
              m1=x-1
              do 10 k=2,m1
10               ga=ga*k
           else
              ga=1.0d+300
           endif
        else
           if (dabs(x).gt.1.0d0) then
              z=dabs(x)
              m=int(z)
              r=1.0d0
              do 15 k=1,m
15               r=r*(z-k)
              z=z-m
           else
              z=x
           endif
           data g/1.0d0,0.5772156649015329d0,
     &          -0.6558780715202538d0, -0.420026350340952d-1,
     &          0.1665386113822915d0,-.421977345555443d-1,
     &          -.96219715278770d-2, .72189432466630d-2,
     &          -.11651675918591d-2, -.2152416741149d-3,
     &          .1280502823882d-3, -.201348547807d-4,
     &          -.12504934821d-5, .11330272320d-5,
     &          -.2056338417d-6, .61160950d-8,
     &          .50020075d-8, -.11812746d-8,
     &          .1043427d-9, .77823d-11,
     &          -.36968d-11, .51d-12,
     &          -.206d-13, -.54d-14, .14d-14, .1d-15/
           gr=g(26)
           do 20 k=25,1,-1
20            gr=gr*z+g(k)
           ga=1.0d0/(gr*z)
           if (dabs(x).gt.1.0d0) then
              ga=ga*r
              if (x.lt.0.0d0) ga=-pi/(x*ga*dsin(pi*x))
           endif
        endif
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
