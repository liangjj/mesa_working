        program mchgm
c
c       =======================================================
c       purpose: this program computes the confluent
c                hypergeometric function m(a,b,x) using
c                subroutine chgm
c       input  : a  --- parameter
c                b  --- parameter ( b <> 0,-1,-2,... )
c                x  --- argument
c       output:  hg --- m(a,b,x)
c       example:
c                   a       b       x          m(a,b,x)
c                 -----------------------------------------
c                  1.5     2.0    20.0     .1208527185d+09
c                  4.5     2.0    20.0     .1103561117d+12
c                 -1.5     2.0    20.0     .1004836854d+05
c                 -4.5     2.0    20.0    -.3936045244d+03
c                  1.5     2.0    50.0     .8231906643d+21
c                  4.5     2.0    50.0     .9310512715d+25
c                 -1.5     2.0    50.0     .2998660728d+16
c                 -4.5     2.0    50.0    -.1806547113d+13
c       =======================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter a, b and x '
        read(*,*)a,b,x
        write(*,*)'   a       b       x          m(a,b,x)'
        write(*,*)' -----------------------------------------'
        call chgm(a,b,x,hg)
        write(*,10)a,b,x,hg
10      format(1x,f5.1,3x,f5.1,3x,f5.1,d20.10)
        end


        subroutine chgm(a,b,x,hg)
c
c       ===================================================
c       purpose: compute confluent hypergeometric function
c                m(a,b,x)
c       input  : a  --- parameter
c                b  --- parameter ( b <> 0,-1,-2,... )
c                x  --- argument
c       output:  hg --- m(a,b,x)
c       routine called: gamma for computing ג(x)
c       ===================================================
c
        implicit double precision (a-h,o-z)
        pi=3.141592653589793d0
        a0=a
        a1=a
        x0=x
        hg=0.0d0
        if (b.eq.0.0d0.or.b.eq.-abs(int(b))) then
           hg=1.0d+300
        else if (a.eq.0.0d0.or.x.eq.0.0d0) then
           hg=1.0d0
        else if (a.eq.-1.0d0) then
           hg=1.0d0-x/b
        else if (a.eq.b) then
           hg=dexp(x)
        else if (a-b.eq.1.0d0) then
           hg=(1.0d0+x/b)*dexp(x)
        else if (a.eq.1.0d0.and.b.eq.2.0d0) then
           hg=(dexp(x)-1.0d0)/x
        else if (a.eq.int(a).and.a.lt.0.0d0) then
           m=int(-a)
           r=1.0d0
           hg=1.0d0
           do 10 k=1,m
              r=r*(a+k-1.0d0)/k/(b+k-1.0d0)*x
10            hg=hg+r
        endif
        if (hg.ne.0.0d0) return
        if (x.lt.0.0d0) then
           a=b-a
           a0=a
           x=dabs(x)
        endif
        if (a.lt.2.0d0) nl=0
        if (a.ge.2.0d0) then
           nl=1
           la=int(a)
           a=a-la-1.0d0
        endif
        do 30 n=0,nl
           if (a0.ge.2.0d0) a=a+1.0d0
           if (x.le.30.0d0+dabs(b).or.a.lt.0.0d0) then
              hg=1.0d0
              rg=1.0d0
              do 15 j=1,500
                 rg=rg*(a+j-1.0d0)/(j*(b+j-1.0d0))*x
                 hg=hg+rg
                 if (dabs(rg/hg).lt.1.0d-15) go to 25
15            continue
           else
              call gamma(a,ta)
              call gamma(b,tb)
              xg=b-a
              call gamma(xg,tba)
              sum1=1.0d0
              sum2=1.0d0
              r1=1.0d0
              r2=1.0d0
              do 20 i=1,8
                 r1=-r1*(a+i-1.0d0)*(a-b+i)/(x*i)
                 r2=-r2*(b-a+i-1.0d0)*(a-i)/(x*i)
                 sum1=sum1+r1
20               sum2=sum2+r2
              hg1=tb/tba*x**(-a)*dcos(pi*a)*sum1
              hg2=tb/ta*dexp(x)*x**(a-b)*sum2
              hg=hg1+hg2
           endif
25         if (n.eq.0) y0=hg
           if (n.eq.1) y1=hg
30      continue
        if (a0.ge.2.0d0) then
           do 35 i=1,la-1
              hg=((2.0d0*a-b+x)*y1+(b-a)*y0)/a
              y0=y1
              y1=hg
35            a=a+1.0d0
        endif
        if (x0.lt.0.0d0) hg=hg*dexp(x0)
        a=a1
        x=x0
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
