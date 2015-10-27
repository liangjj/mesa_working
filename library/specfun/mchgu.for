        program mchgu
c
c       =======================================================
c       purpose: this program computes the confluent
c                hypergeometric function u(a,b,x) using
c                subroutine chgu
c       input  : a  --- parameter
c                b  --- parameter
c                x  --- argument  ( x ע 0 )
c       output:  hu --- u(a,b,x)
c                md --- method code
c       example:
c                a       b       x        u(a,b,x)
c             --------------------------------------
c              -2.5     2.5     5.0     -9.02812446
c              -1.5     2.5     5.0      2.15780560
c               -.5     2.5     5.0      1.76649370
c                .0     2.5     5.0      1.00000000
c                .5     2.5     5.0       .49193496
c               1.5     2.5     5.0       .08944272
c               2.5     2.5     5.0       .01239387
c
c                a       b       x        u(a,b,x)
c             --------------------------------------
c              -2.5     5.0    10.0     -2.31982196
c              -1.5     5.0    10.0      8.65747115
c               -.5     5.0    10.0      2.37997143
c                .0     5.0    10.0      1.00000000
c                .5     5.0    10.0       .38329536
c               1.5     5.0    10.0       .04582817
c               2.5     5.0    10.0       .00444535
c       =======================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'pleas enter a, b and x '
        read(*,*)a,b,x
        write(*,*)'   a       b       x        u(a,b,x)'
        write(*,*)'--------------------------------------'
        call chgu(a,b,x,hu,md)
        write(*,10)a,b,x,hu
10      format(1x,f5.1,3x,f5.1,3x,f5.1,e15.8)
        end


        subroutine chgu(a,b,x,hu,md)
c
c       =======================================================
c       purpose: compute the confluent hypergeometric function
c                u(a,b,x)
c       input  : a  --- parameter
c                b  --- parameter
c                x  --- argument  ( x > 0 )
c       output:  hu --- u(a,b,x)
c                md --- method code
c       routines called:
c            (1) chgus for small x ( md=1 )
c            (2) chgul for large x ( md=2 )
c            (3) chgubi for integer b ( md=3 )
c            (4) chguit for numerical integration ( md=4 )
c       =======================================================
c
        implicit double precision (a-h,o-z)
        logical il1,il2,il3,bl1,bl2,bl3,bn
        aa=a-b+1.0d0
        il1=a.eq.int(a).and.a.le.0.0
        il2=aa.eq.int(aa).and.aa.le.0.0
        il3=abs(a*(a-b+1.0))/x.le.2.0
        bl1=x.le.5.0.or.(x.le.10.0.and.a.le.2.0)
        bl2=(x.gt.5.0.and.x.le.12.5).and.(a.ge.1.0.and.b.ge.a+4.0)
        bl3=x.gt.12.5.and.a.ge.5.0.and.b.ge.a+5.0
        bn=b.eq.int(b).and.b.ne.0.0
        id1=-100
        if (b.ne.int(b)) then
           call chgus(a,b,x,hu,id1)
           md=1
           if (id1.ge.6) return
           hu1=hu
        endif
        if (il1.or.il2.or.il3) then
           call chgul(a,b,x,hu,id)
           md=2
           if (id.ge.6) return
           if (id1.gt.id) then
              md=1
              id=id1
              hu=hu1
           endif
        endif
        if (a.ge.0.0) then
           if (bn.and.(bl1.or.bl2.or.bl3)) then
              call chgubi(a,b,x,hu,id)
              md=3
           else
              call chguit(a,b,x,hu,id)
              md=4
           endif
        else
           if (b.le.a) then
              a00=a
              b00=b
              a=a-b+1.0d0
              b=2.0d0-b
              call chguit(a,b,x,hu,id)
              hu=x**(1.0d0-b00)*hu
              a=a00
              b=b00
              md=4
           else if (bn.and.(.not.il1)) then
              call chgubi(a,b,x,hu,id)
              md=3
           endif
        endif
        if (id.lt.6) write(*,*)'no accurate result obtained'
        return
        end


        subroutine chgus(a,b,x,hu,id)
c
c       ======================================================
c       purpose: compute confluent hypergeometric function
c                u(a,b,x) for small argument x
c       input  : a  --- parameter
c                b  --- parameter ( b <> 0,-1,-2,...)
c                x  --- argument
c       output:  hu --- u(a,b,x)
c                id --- estimated number of significant digits
c       routine called: gamma for computing gamma function
c       ======================================================
c
        implicit double precision (a-h,o-z)
        id=-100
        pi=3.141592653589793d0
        call gamma(a,ga)
        call gamma(b,gb)
        xg1=1.0d0+a-b
        call gamma(xg1,gab)
        xg2=2.0d0-b
        call gamma(xg2,gb2)
        hu0=pi/dsin(pi*b)
        r1=hu0/(gab*gb)
        r2=hu0*x**(1.0d0-b)/(ga*gb2)
        hu=r1-r2
        hmax=0.0d0
        hmin=1.0d+300
        do 10 j=1,150
           r1=r1*(a+j-1.0d0)/(j*(b+j-1.0d0))*x
           r2=r2*(a-b+j)/(j*(1.0d0-b+j))*x
           hu=hu+r1-r2
           hua=dabs(hu)
           if (hua.gt.hmax) hmax=hua
           if (hua.lt.hmin) hmin=hua
           if (dabs(hu-h0).lt.dabs(hu)*1.0d-15) go to 15
10         h0=hu
15      d1=log10(hmax)
        if (hmin.ne.0.0) d2=log10(hmin)
        id=15-abs(d1-d2)
        return
        end


        subroutine chgul(a,b,x,hu,id)
c
c       =======================================================
c       purpose: compute the confluent hypergeometric function
c                u(a,b,x) for large argument x
c       input  : a  --- parameter
c                b  --- parameter
c                x  --- argument
c       output:  hu --- u(a,b,x)
c                id --- estimated number of significant digits
c       =======================================================
c
        implicit double precision (a-h,o-z)
        logical il1,il2
        id=-100
        aa=a-b+1.0d0
        il1=a.eq.int(a).and.a.le.0.0
        il2=aa.eq.int(aa).and.aa.le.0.0
        if (il1) nm=abs(a)
        if (il2) nm=abs(aa)
        if (il1.or.il2) then
           hu=1.0d0
           r=1.0d0
           do 10 k=1,nm
              r=-r*(a+k-1.0d0)*(a-b+k)/(k*x)
              hu=hu+r
10         continue
           hu=x**(-a)*hu
           id=10
        else
           hu=1.0d0
           r=1.0d0
           do 15 k=1,25
              r=-r*(a+k-1.0d0)*(a-b+k)/(k*x)
              ra=dabs(r)
              if (k.gt.5.and.ra.ge.r0.or.ra.lt.1.0d-15) go to 20
              r0=ra
15            hu=hu+r
20         id=abs(log10(ra))
           hu=x**(-a)*hu
        endif
        return
        end


        subroutine chgubi(a,b,x,hu,id)
c
c       ======================================================
c       purpose: compute confluent hypergeometric function
c                u(a,b,x) with integer b ( b = ס1,ס2,... )
c       input  : a  --- parameter
c                b  --- parameter
c                x  --- argument
c       output:  hu --- u(a,b,x)
c                id --- estimated number of significant digits
c       routines called:
c            (1) gamma for computing gamma function ג(x)
c            (2) psi for computing psi function
c       ======================================================
c
        implicit double precision (a-h,o-z)
        id=-100
        el=0.5772156649015329d0
        n=abs(b-1)
        rn1=1.0d0
        rn=1.0d0
        do 10 j=1,n
           rn=rn*j
           if (j.eq.n-1) rn1=rn
10      continue
        call psi(a,ps)
        call gamma(a,ga)
        if (b.gt.0.0) then
           a0=a
           a1=a-n
           a2=a1
           call gamma(a1,ga1)
           ua=(-1)**(n-1)/(rn*ga1)
           ub=rn1/ga*x**(-n)
        else
           a0=a+n
           a1=a0
           a2=a
           call gamma(a1,ga1)
           ua=(-1)**(n-1)/(rn*ga)*x**n
           ub=rn1/ga1
        endif
        hm1=1.0d0
        r=1.0d0
        hmax=0.0d0
        hmin=1.0d+300
        do 15 k=1,150
           r=r*(a0+k-1.0d0)*x/((n+k)*k)
           hm1=hm1+r
           hu1=dabs(hm1)
           if (hu1.gt.hmax) hmax=hu1
           if (hu1.lt.hmin) hmin=hu1
           if (dabs(hm1-h0).lt.dabs(hm1)*1.0d-15) go to 20
15         h0=hm1
20      da1=log10(hmax)
        if (hmin.ne.0.0) da2=log10(hmin)
        id=15-abs(da1-da2)
        hm1=hm1*dlog(x)
        s0=0.0d0
        do 25 m=1,n
           if (b.ge.0.0) s0=s0-1.0d0/m
25         if (b.lt.0.0) s0=s0+(1.0d0-a)/(m*(a+m-1.0d0))
        hm2=ps+2.0d0*el+s0
        r=1.0d0
        hmax=0.0d0
        hmin=1.0d+300
        do 50 k=1,150
           s1=0.0d0
           s2=0.0d0
           if (b.gt.0.0) then
              do 30 m=1,k
30               s1=s1-(m+2.0d0*a-2.0d0)/(m*(m+a-1.0d0))
              do 35 m=1,n
35               s2=s2+1.0d0/(k+m)
           else
              do 40 m=1,k+n
40               s1=s1+(1.0d0-a)/(m*(m+a-1.0d0))
              do 45 m=1,k
45               s2=s2+1.0d0/m
           endif
           hw=2.0d0*el+ps+s1-s2
           r=r*(a0+k-1.0d0)*x/((n+k)*k)
           hm2=hm2+r*hw
           hu2=dabs(hm2)
           if (hu2.gt.hmax) hmax=hu2
           if (hu2.lt.hmin) hmin=hu2
           if (dabs((hm2-h0)/hm2).lt.1.0d-15) go to 55
50         h0=hm2
55      db1=log10(hmax)
        if (hmin.ne.0.0) db2=log10(hmin)
        id1=15-abs(db1-db2)
        if (id1.lt.id) id=id1
        hm3=1.0d0
        if (n.eq.0) hm3=0.0d0
        r=1.0d0
        do 60 k=1,n-1
           r=r*(a2+k-1.0d0)/((k-n)*k)*x
60         hm3=hm3+r
        sa=ua*(hm1+hm2)
        sb=ub*hm3
        hu=sa+sb
        if (sa.ne.0.0) id1=int(log10(abs(sa)))
        if (hu.ne.0.0) id2=int(log10(abs(hu)))
        if (sa*sb.lt.0.0) id=id-abs(id1-id2)
        return
        end


        subroutine chguit(a,b,x,hu,id)
c
c       ======================================================
c       purpose: compute hypergeometric function u(a,b,x) by
c                using gaussian-legendre integration (n=60)
c       input  : a  --- parameter ( a > 0 )
c                b  --- parameter
c                x  --- argument ( x > 0 )
c       output:  hu --- u(a,b,z)
c                id --- estimated number of significant digits
c       routine called: gamma for computing ג(x)
c       ======================================================
c
        implicit double precision (a-h,o-z)
        dimension t(30),w(30)
        data t/ .259597723012478d-01, .778093339495366d-01,
     &          .129449135396945d+00, .180739964873425d+00,
     &          .231543551376029d+00, .281722937423262d+00,
     &          .331142848268448d+00, .379670056576798d+00,
     &          .427173741583078d+00, .473525841761707d+00,
     &          .518601400058570d+00, .562278900753945d+00,
     &          .604440597048510d+00, .644972828489477d+00,
     &          .683766327381356d+00, .720716513355730d+00,
     &          .755723775306586d+00, .788693739932264d+00,
     &          .819537526162146d+00, .848171984785930d+00,
     &          .874519922646898d+00, .898510310810046d+00,
     &          .920078476177628d+00, .939166276116423d+00,
     &          .955722255839996d+00, .969701788765053d+00,
     &          .981067201752598d+00, .989787895222222d+00,
     &          .995840525118838d+00, .999210123227436d+00/
        data w/ .519078776312206d-01, .517679431749102d-01,
     &          .514884515009810d-01, .510701560698557d-01,
     &          .505141845325094d-01, .498220356905502d-01,
     &          .489955754557568d-01, .480370318199712d-01,
     &          .469489888489122d-01, .457343797161145d-01,
     &          .443964787957872d-01, .429388928359356d-01,
     &          .413655512355848d-01, .396806954523808d-01,
     &          .378888675692434d-01, .359948980510845d-01,
     &          .340038927249464d-01, .319212190192963d-01,
     &          .297524915007890d-01, .275035567499248d-01,
     &          .251804776215213d-01, .227895169439978d-01,
     &          .203371207294572d-01, .178299010142074d-01,
     &          .152746185967848d-01, .126781664768159d-01,
     &          .100475571822880d-01, .738993116334531d-02,
     &          .471272992695363d-02, .202681196887362d-02/
        id=7
        a1=a-1.0d0
        b1=b-a-1.0d0
        c=12.0/x
        do 20 m=10,100,5
           hu1=0.0d0
           g=0.5d0*c/m
           d=g
           do 15 j=1,m
              s=0.0d0
              do 10 k=1,30
                 t1=d+g*t(k)
                 t2=d-g*t(k)
                 f1=dexp(-x*t1)*t1**a1*(1.0d0+t1)**b1
                 f2=dexp(-x*t2)*t2**a1*(1.0d0+t2)**b1
                 s=s+w(k)*(f1+f2)
10            continue
              hu1=hu1+s*g
              d=d+2.0d0*g
15         continue
           if (dabs(1.0d0-hu0/hu1).lt.1.0d-7) go to 25
           hu0=hu1
20      continue
25      call gamma(a,ga)
        hu1=hu1/ga
        do 40 m=2,10,2
           hu2=0.0d0
           g=0.5d0/m
           d=g
           do 35 j=1,m
              s=0.0d0
              do 30 k=1,30
                 t1=d+g*t(k)
                 t2=d-g*t(k)
                 t3=c/(1.0d0-t1)
                 t4=c/(1.0d0-t2)
                 f1=t3*t3/c*dexp(-x*t3)*t3**a1*(1.0d0+t3)**b1
                 f2=t4*t4/c*dexp(-x*t4)*t4**a1*(1.0d0+t4)**b1
                 s=s+w(k)*(f1+f2)
30            continue
              hu2=hu2+s*g
              d=d+2.0d0*g
35         continue
           if (dabs(1.0d0-hu0/hu2).lt.1.0d-7) go to 45
           hu0=hu2
40      continue
45      call gamma(a,ga)
        hu2=hu2/ga
        hu=hu1+hu2
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
