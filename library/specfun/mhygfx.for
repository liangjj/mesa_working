        program mhygfx
c
c       ============================================================
c       purpose: this program computes the hypergeometric function 
c                f(a,b,c,x) using subroutine hygfx
c       input :  a --- parameter
c                b --- parameter
c                c --- parameter, c <> 0,-1,-2,...
c                x --- argument ( x ף 1 )
c       output:  hf --- f(a,b,c,x)
c       example:
c              b = 3.30,  c = 6.70
c              a     f(a,b,c,.25)     f(a,b,c,.55)    f(a,b,c,.85)
c            ------------------------------------------------------
c            -2.5   .72356129d+00    .46961432d+00   .29106096d+00
c            -0.5   .93610145d+00    .85187390d+00   .75543187d+00
c             0.5   .10689695d+01    .11795358d+01   .13510497d+01
c             2.5   .14051563d+01    .23999063d+01   .57381566d+01
c
c              a = 3.30,  b = 6.70
c              c     f(a,b,c,.25)     f(a,b,c,.55)    f(a,b,c,.85)
c            ------------------------------------------------------
c            -5.5   .15090670d+05    .10170778d+11   .58682088d+19
c            -0.5  -.21631479d+04   -.30854772d+07  -.10217370d+13
c             0.5   .26451677d+03    .11967860d+06   .92370648d+10
c             4.5   .41946916d+01    .58092729d+02   .20396914d+05
c       ============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter a,b,c and x '
        read(*,*)a,b,c,x
        call hygfx(a,b,c,x,hf)
        write(*,10)a,b,c,x
        write(*,20)hf
10      format(1x,'a =',f5.2,',    ','b =',f5.2,',   ','c =',f5.2,
     &         ',    ','x =',f5.2)
20      format(1x,'f(a,b,c,x)=',d16.8)
        end


        subroutine hygfx(a,b,c,x,hf)
c
c       ====================================================
c       purpose: compute hypergeometric function f(a,b,c,x)
c       input :  a --- parameter
c                b --- parameter
c                c --- parameter, c <> 0,-1,-2,...
c                x --- argument   ( x < 1 )
c       output:  hf --- f(a,b,c,x)
c       routines called:
c            (1) gamma for computing gamma function
c            (2) psi for computing psi function
c       ====================================================
c
        implicit double precision (a-h,o-z)
        logical l0,l1,l2,l3,l4,l5
        pi=3.141592653589793d0
        el=.5772156649015329d0
        l0=c.eq.int(c).and.c.lt.0.0
        l1=1.0d0-x.lt.1.0d-15.and.c-a-b.le.0.0
        l2=a.eq.int(a).and.a.lt.0.0
        l3=b.eq.int(b).and.b.lt.0.0
        l4=c-a.eq.int(c-a).and.c-a.le.0.0
        l5=c-b.eq.int(c-b).and.c-b.le.0.0
        if (l0.or.l1) then
           write(*,*)'the hypergeometric series is divergent'
           return
        endif
        eps=1.0d-15
        if (x.gt.0.95) eps=1.0d-8
        if (x.eq.0.0.or.a.eq.0.0.or.b.eq.0.0) then
           hf=1.0d0
           return
        else if (1.0d0-x.eq.eps.and.c-a-b.gt.0.0) then
           call gamma(c,gc)
           call gamma(c-a-b,gcab)
           call gamma(c-a,gca)
           call gamma(c-b,gcb)
           hf=gc*gcab/(gca*gcb)
           return
        else if (1.0d0+x.le.eps.and.dabs(c-a+b-1.0).le.eps) then
           g0=dsqrt(pi)*2.0d0**(-a)
           call gamma(c,g1)
           call gamma(1.0d0+a/2.0-b,g2)
           call gamma(0.5d0+0.5*a,g3)
           hf=g0*g1/(g2*g3)
           return
        else if (l2.or.l3) then
           if (l2) nm=int(abs(a))
           if (l3) nm=int(abs(b))
           hf=1.0d0
           r=1.0d0
           do 10 k=1,nm
              r=r*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*x
10            hf=hf+r
           return
        else if (l4.or.l5) then
           if (l4) nm=int(abs(c-a))
           if (l5) nm=int(abs(c-b))
           hf=1.0d0
           r=1.0d0
           do 15 k=1,nm
              r=r*(c-a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c+k-1.0d0))*x
15            hf=hf+r
           hf=(1.0d0-x)**(c-a-b)*hf
           return
        endif
        aa=a
        bb=b
        x1=x
        if (x.lt.0.0d0) then
           x=x/(x-1.0d0)
           if (c.gt.a.and.b.lt.a.and.b.gt.0.0) then
              a=bb
              b=aa
           endif
           b=c-b
        endif
        if (x.ge.0.75d0) then
           gm=0.0d0
           if (dabs(c-a-b-int(c-a-b)).lt.1.0d-15) then
              m=int(c-a-b)
              call gamma(a,ga)
              call gamma(b,gb)
              call gamma(c,gc)
              call gamma(a+m,gam)
              call gamma(b+m,gbm)
              call psi(a,pa)
              call psi(b,pb)
              if (m.ne.0) gm=1.0d0
              do 30 j=1,abs(m)-1
30               gm=gm*j
              rm=1.0d0
              do 35 j=1,abs(m)
35               rm=rm*j
              f0=1.0d0
              r0=1.0d0
              r1=1.0d0
              sp0=0.d0
              sp=0.0d0
              if (m.ge.0) then
                 c0=gm*gc/(gam*gbm)
                 c1=-gc*(x-1.0d0)**m/(ga*gb*rm)
                 do 40 k=1,m-1
                    r0=r0*(a+k-1.0d0)*(b+k-1.0)/(k*(k-m))*(1.0-x)
40                  f0=f0+r0
                 do 45 k=1,m
45                  sp0=sp0+1.0d0/(a+k-1.0)+1.0/(b+k-1.0)-1.0/k
                 f1=pa+pb+sp0+2.0d0*el+dlog(1.0d0-x)
                 do 55 k=1,250
                    sp=sp+(1.0d0-a)/(k*(a+k-1.0))+(1.0-b)/(k*(b+k-1.0))
                    sm=0.0d0
                    do 50 j=1,m
50                     sm=sm+(1.0d0-a)/((j+k)*(a+j+k-1.0))+1.0/
     &                    (b+j+k-1.0)
                    rp=pa+pb+2.0d0*el+sp+sm+dlog(1.0d0-x)
                    r1=r1*(a+m+k-1.0d0)*(b+m+k-1.0)/(k*(m+k))*(1.0-x)
                    f1=f1+r1*rp
                    if (dabs(f1-hw).lt.dabs(f1)*eps) go to 60
55                  hw=f1
60               hf=f0*c0+f1*c1
              else if (m.lt.0) then
                 m=-m
                 c0=gm*gc/(ga*gb*(1.0d0-x)**m)
                 c1=-(-1)**m*gc/(gam*gbm*rm)
                 do 65 k=1,m-1
                    r0=r0*(a-m+k-1.0d0)*(b-m+k-1.0)/(k*(k-m))*(1.0-x)
65                  f0=f0+r0
                 do 70 k=1,m
70                  sp0=sp0+1.0d0/k
                 f1=pa+pb-sp0+2.0d0*el+dlog(1.0d0-x)
                 do 80 k=1,250
                    sp=sp+(1.0d0-a)/(k*(a+k-1.0))+(1.0-b)/(k*(b+k-1.0))
                    sm=0.0d0
                    do 75 j=1,m
75                     sm=sm+1.0d0/(j+k)
                    rp=pa+pb+2.0d0*el+sp-sm+dlog(1.0d0-x)
                    r1=r1*(a+k-1.0d0)*(b+k-1.0)/(k*(m+k))*(1.0-x)
                    f1=f1+r1*rp
                    if (dabs(f1-hw).lt.dabs(f1)*eps) go to 85
80                  hw=f1
85               hf=f0*c0+f1*c1
              endif
           else
              call gamma(a,ga)
              call gamma(b,gb)
              call gamma(c,gc)
              call gamma(c-a,gca)
              call gamma(c-b,gcb)
              call gamma(c-a-b,gcab)
              call gamma(a+b-c,gabc)
              c0=gc*gcab/(gca*gcb)
              c1=gc*gabc/(ga*gb)*(1.0d0-x)**(c-a-b)
              hf=0.0d0
              r0=c0
              r1=c1
              do 90 k=1,250
                 r0=r0*(a+k-1.0d0)*(b+k-1.0)/(k*(a+b-c+k))*(1.0-x)
                 r1=r1*(c-a+k-1.0d0)*(c-b+k-1.0)/(k*(c-a-b+k))
     &              *(1.0-x)
                 hf=hf+r0+r1
                 if (dabs(hf-hw).lt.dabs(hf)*eps) go to 95
90               hw=hf
95            hf=hf+c0+c1
           endif
        else
           a0=1.0d0
           if (c.gt.a.and.c.lt.2.0d0*a.and.
     &         c.gt.b.and.c.lt.2.0d0*b) then
              a0=(1.0d0-x)**(c-a-b)
              a=c-a
              b=c-b
           endif
           hf=1.0d0
           r=1.0d0
           do 100 k=1,250
              r=r*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*x
              hf=hf+r
              if (dabs(hf-hw).le.dabs(hf)*eps) go to 105
100           hw=hf
105        hf=a0*hf
        endif
        if (x1.lt.0.0d0) then
           x=x1
           c0=1.0d0/(1.0d0-x)**aa
           hf=c0*hf
        endif
        a=aa
        b=bb
        if (k.gt.120) write(*,115)
115     format(1x,'warning! you should check the accuracy')
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


        subroutine psi(x,ps)
c
c       ======================================
c       purpose: compute psi function
c       input :  x  --- argument of psi(x)
c       output:  ps --- psi(x)
c       ======================================
c
        implicit double precision (a-h,o-z)
        xa=dabs(x)
        pi=3.141592653589793d0
        el=.5772156649015329d0
        s=0.0d0
        if (x.eq.int(x).and.x.le.0.0) then
           ps=1.0d+300
           return
        else if (xa.eq.int(xa)) then
           n=xa
           do 10 k=1 ,n-1
10            s=s+1.0d0/k
           ps=-el+s
        else if (xa+.5.eq.int(xa+.5)) then
           n=xa-.5
           do 20 k=1,n
20            s=s+1.0/(2.0d0*k-1.0d0)
           ps=-el+2.0d0*s-1.386294361119891d0
        else
           if (xa.lt.10.0) then
              n=10-int(xa)
              do 30 k=0,n-1
30               s=s+1.0d0/(xa+k)
              xa=xa+n
           endif
           x2=1.0d0/(xa*xa)
           a1=-.8333333333333d-01
           a2=.83333333333333333d-02
           a3=-.39682539682539683d-02
           a4=.41666666666666667d-02
           a5=-.75757575757575758d-02
           a6=.21092796092796093d-01
           a7=-.83333333333333333d-01
           a8=.4432598039215686d0
           ps=dlog(xa)-.5d0/xa+x2*(((((((a8*x2+a7)*x2+
     &        a6)*x2+a5)*x2+a4)*x2+a3)*x2+a2)*x2+a1)
           ps=ps-s
        endif
        if (x.lt.0.0) ps=ps-pi*dcos(pi*x)/dsin(pi*x)-1.0d0/x
        return
        end
