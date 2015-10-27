        program mjyzo
c
c       ==========================================================
c       purpose: this program computes the zeros of bessel 
c                functions jn(x), yn(x), and their derivatives 
c                using subroutine jyzo
c       input :  n --- order of bessel functions ( n ó 100 )
c                nt --- number of zeros
c       output:  rj0(m) --- m-th zero of jn(x),  m=1,2,...,nt
c                rj1(m) --- m-th zero of jn'(x), m=1,2,...,nt
c                ry0(m) --- m-th zero of yn(x),  m=1,2,...,nt
c                ry1(m) --- m-th zero of yn'(x), m=1,2,...,nt
c       example: n = 1, nt =5
c
c      zeros of bessel funcions jn(x), yn(x) and their derivatives
c                                 ( n = 1 )
c       m       jnm           j'nm          ynm           y'nm
c      -----------------------------------------------------------
c       1     3.8317060     1.8411838     2.1971413     3.6830229
c       2     7.0155867     5.3314428     5.4296810     6.9415000
c       3    10.1734681     8.5363164     8.5960059    10.1234047
c       4    13.3236919    11.7060049    11.7491548    13.2857582
c       5    16.4706301    14.8635886    14.8974421    16.4400580
c       ==========================================================
c
        implicit double precision (a-h,o-z)
        dimension rj0(101),rj1(101),ry0(101),ry1(101)
        write(*,*)'please enter n and nt '
        read(*,*)n,nt
        write(*,*)
        call jyzo(n,nt,rj0,rj1,ry0,ry1)
        write(*,30)
        write(*,40)n
        write(*,*)'  m       jnm           j''nm          ynm',
     &            '           y''nm'
        write(*,*)' ----------------------------------------',
     &            '-------------------'
        do 10 m=1,nt
10         write(*,50)m,rj0(m),rj1(m),ry0(m),ry1(m)
30      format(2x,'zeros of bessel funcions jn(x), yn(x)',
     &         ' and their derivatives')
40      format(30x,'( n =',i2,' )')
50      format(1x,i3,4f14.7)
        end


        subroutine jyzo(n,nt,rj0,rj1,ry0,ry1)
c
c       ======================================================
c       purpose: compute the zeros of bessel functions jn(x),
c                yn(x), and their derivatives
c       input :  n  --- order of bessel functions ( n ó 101 )
c                nt --- number of zeros (roots)
c       output:  rj0(l) --- l-th zero of jn(x),  l=1,2,...,nt
c                rj1(l) --- l-th zero of jn'(x), l=1,2,...,nt
c                ry0(l) --- l-th zero of yn(x),  l=1,2,...,nt
c                ry1(l) --- l-th zero of yn'(x), l=1,2,...,nt
c       routine called: jyndd for computing jn(x), yn(x), and
c                       their first and second derivatives
c       ======================================================
c
        implicit double precision (a-h,o-z)
        dimension rj0(nt),rj1(nt),ry0(nt),ry1(nt)
        if (n.le.20) then
           x=2.82141+1.15859*n
        else
           x=n+1.85576*n**0.33333+1.03315/n**0.33333
        endif
        l=0
10      x0=x
        call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
        x=x-bjn/djn
        if (dabs(x-x0).gt.1.0d-9) go to 10
        l=l+1
        rj0(l)=x
        x=x+3.1416+(0.0972+0.0679*n-0.000354*n**2)/l
        if (l.lt.nt) go to 10
        if (n.le.20) then
           x=0.961587+1.07703*n
        else
           x=n+0.80861*n**0.33333+0.07249/n**0.33333
        endif
        if (n.eq.0) x=3.8317
        l=0
15      x0=x
        call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
        x=x-djn/fjn
        if (dabs(x-x0).gt.1.0d-9) go to 15
        l=l+1
        rj1(l)=x
        x=x+3.1416+(0.4955+0.0915*n-0.000435*n**2)/l
        if (l.lt.nt) go to 15
        if (n.le.20) then
           x=1.19477+1.08933*n
        else
           x=n+0.93158*n**0.33333+0.26035/n**0.33333
        endif           
        l=0
20      x0=x
        call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
        x=x-byn/dyn
        if (dabs(x-x0).gt.1.0d-9) go to 20
        l=l+1
        ry0(l)=x
        x=x+3.1416+(0.312+0.0852*n-0.000403*n**2)/l
        if (l.lt.nt) go to 20
        if (n.le.20) then
           x=2.67257+1.16099*n
        else
           x=n+1.8211*n**0.33333+0.94001/n**0.33333
        endif  
        l=0
25      x0=x
        call jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
        x=x-dyn/fyn
        if (dabs(x-x0).gt.1.0d-9) go to 25
        l=l+1
        ry1(l)=x
        x=x+3.1416+(0.197+0.0643*n-0.000286*n**2)/l 
        if (l.lt.nt) go to 25
        return
        end


        subroutine jyndd(n,x,bjn,djn,fjn,byn,dyn,fyn)
c
c       ===========================================================
c       purpose: compute bessel functions jn(x) and yn(x), and
c                their first and second derivatives 
c       input:   x   ---  argument of jn(x) and yn(x) ( x > 0 )
c                n   ---  order of jn(x) and yn(x)
c       output:  bjn ---  jn(x)
c                djn ---  jn'(x)
c                fjn ---  jn"(x)
c                byn ---  yn(x)
c                dyn ---  yn'(x)
c                fyn ---  yn"(x)
c       ===========================================================
c
        implicit double precision (a-h,o-z)
        dimension bj(102),by(102)
        do 10 nt=1,900
           mt=int(0.5*log10(6.28*nt)-nt*log10(1.36*dabs(x)/nt))
           if (mt.gt.20) go to 15
10      continue
15      m=nt
        bs=0.0d0
        f0=0.0d0
        f1=1.0d-35
        su=0.0d0
        do 20 k=m,0,-1
           f=2.0d0*(k+1.0d0)*f1/x-f0
           if (k.le.n+1) bj(k+1)=f
           if (k.eq.2*int(k/2)) then
              bs=bs+2.0d0*f
              if (k.ne.0) su=su+(-1)**(k/2)*f/k
           endif
           f0=f1
20         f1=f
        do 25 k=0,n+1
25         bj(k+1)=bj(k+1)/(bs-f)
        bjn=bj(n+1)
        ec=0.5772156649015329d0
        e0=0.3183098861837907d0
        s1=2.0d0*e0*(dlog(x/2.0d0)+ec)*bj(1)
        f0=s1-8.0d0*e0*su/(bs-f)
        f1=(bj(2)*f0-2.0d0*e0/x)/bj(1)
        by(1)=f0
        by(2)=f1
        do 30 k=2,n+1
           f=2.0d0*(k-1.0d0)*f1/x-f0
           by(k+1)=f
           f0=f1
30         f1=f
        byn=by(n+1)
        djn=-bj(n+2)+n*bj(n+1)/x
        dyn=-by(n+2)+n*by(n+1)/x
        fjn=(n*n/(x*x)-1.0d0)*bjn-djn/x
        fyn=(n*n/(x*x)-1.0d0)*byn-dyn/x
        return
        end


 
