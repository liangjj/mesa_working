        program mhygfz
c
c       ============================================================
c       purpose: this program computes hypergeometric function for
c                a complex argument, f(a,b,c,z), using subroutine
c                hygfz
c       input :  a --- parameter
c                b --- parameter
c                c --- parameter,  c <> 0,-1,-2,...
c                x --- real part of complex argument z
c                y --- imaginary part of complex argument z
c                      ( z = x+iy )
c       output:  zhf --- f(a,b,c,z)
c       examples:
c     a     b     c    z = x+ iy             f(a,b,c,z)
c   --------------------------------------------------------------
c    3.2   1.8   6.7   1.0+0.0 i    .54689992d+01+.00000000d+00 i
c    3.2  -1.8   6.7   1.0+0.0 i    .33750635d+00+.00000000d+00 i
c   -5.0   3.3   6.7   5.2+4.8 i    .11682745d+03+.60389104d+03 i
c    3.3  -6.0   3.7   5.2-4.8 i    .17620425d+05+.38293812d+05 i
c   -7.0   3.3  -3.7   5.2-4.8 i   -.11772779d+11-.14382286d+11 i
c    4.3  -8.0  -3.7   5.2+4.8 i    .13161188d+13-.10129870d+12 i
c    3.3   5.8   6.7   0.2+0.1 i    .17330557d+01+.63401030d+00 i
c    3.5  -2.4   6.7   0.2+0.5 i    .64762241d+00-.52110507d+00 i
c    3.3   4.3   6.7   0.8+0.3 i   -.14830086d+01+.83744258d+01 i
c    7.0   5.0   4.1   3.0-1.0 i   -.40376095d-02-.29566326d-02 i
c    5.0   7.0   4.1   3.0-1.0 i   -.40376095d-02-.29566326d-02 i
c    3.5   1.2   9.7   0.6+0.9 i    .10343044d+01+.54473814d+00 i
c    2.1   5.4   9.7   0.5+0.7 i    .68850442d+00+.12274187d+01 i
c    8.7   3.2   6.7   0.5+0.7 i   -.90046505d+00-.11198900d+01 i
c    8.7   2.7   6.7   0.6+0.9 i   -.46083890d+00-.54575701d+00 i
c       ============================================================
c
        implicit double precision (a-h,o-y)
        implicit complex *16 (z)
        write(*,*)'please enter a,b,c,x and y '
        read(*,*)a,b,c,x,y
        z=cmplx(x,y)
        call hygfz(a,b,c,z,zhf)
        write(*,*)
        write(*,*)'     a      b      c      x      y',
     &            '          re[f]           im[f]'
        write(*,*)'   --------------------------------',
     &            '-----------------------------------'
        write(*,10)a,b,c,x,y,zhf
10      format(1x,5f7.1,2x,2d16.8)
        end


        subroutine hygfz(a,b,c,z,zhf)
c
c       ======================================================
c       purpose: compute the hypergeometric function for a 
c                complex argument, f(a,b,c,z)
c       input :  a --- parameter
c                b --- parameter
c                c --- parameter,  c <> 0,-1,-2,...
c                z --- complex argument
c       output:  zhf --- f(a,b,c,z)
c       routines called:
c            (1) gamma for computing gamma function
c            (2) psi for computing psi function
c       ======================================================
c
        implicit double precision (a-h,o-y)
        implicit complex *16 (z)
        logical l0,l1,l2,l3,l4,l5,l6
        x=real(z)
        y=dimag(z)
        eps=1.0d-15
        l0=c.eq.int(c).and.c.lt.0.0d0
        l1=dabs(1.0d0-x).lt.eps.and.y.eq.0.0d0.and.c-a-b.le.0.0d0
        l2=cdabs(z+1.0d0).lt.eps.and.dabs(c-a+b-1.0d0).lt.eps
        l3=a.eq.int(a).and.a.lt.0.0d0
        l4=b.eq.int(b).and.b.lt.0.0d0
        l5=c-a.eq.int(c-a).and.c-a.le.0.0d0
        l6=c-b.eq.int(c-b).and.c-b.le.0.0d0
        aa=a
        bb=b
        a0=cdabs(z)
        if (a0.gt.0.95d0) eps=1.0d-8
        pi=3.141592653589793d0
        el=.5772156649015329d0
        if (l0.or.l1) then
           write(*,*)'the hypergeometric series is divergent'
           return
        endif
        if (a0.eq.0.0d0.or.a.eq.0.0d0.or.b.eq.0.0d0) then
           zhf=(1.0d0,0.0d0)
        else if (z.eq.1.0d0.and.c-a-b.gt.0.0d0) then
           call gamma(c,gc)
           call gamma(c-a-b,gcab)
           call gamma(c-a,gca)
           call gamma(c-b,gcb)
           zhf=gc*gcab/(gca*gcb)
        else if (l2) then
           g0=dsqrt(pi)*2.0d0**(-a)
           call gamma(c,g1)
           call gamma(1.0d0+a/2.0d0-b,g2)
           call gamma(0.5d0+0.5d0*a,g3)
           zhf=g0*g1/(g2*g3)
        else if (l3.or.l4) then
           if (l3) nm=int(abs(a))
           if (l4) nm=int(abs(b))
           zhf=(1.0d0,0.0d0)
           zr=(1.0d0,0.0d0)
           do 10 k=1,nm
              zr=zr*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*z
10            zhf=zhf+zr
        else if (l5.or.l6) then
           if (l5) nm=int(abs(c-a))
           if (l6) nm=int(abs(c-b))
           zhf=(1.0d0,0.0d0)
           zr=(1.0d0,0.0d0)
           do 15 k=1,nm
              zr=zr*(c-a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c+k-1.0d0))*z
15            zhf=zhf+zr
           zhf=(1.0d0-z)**(c-a-b)*zhf
        else if (a0.le.1.0d0) then
           if (x.lt.0.0d0) then
              z1=z/(z-1.0d0)
              if (c.gt.a.and.b.lt.a.and.b.gt.0.0) then  
                 a=bb
                 b=aa
              endif
              zc0=1.0d0/((1.0d0-z)**a)
              zhf=(1.0d0,0.0d0)
              zr0=(1.0d0,0.0d0)
              do 20 k=1,500
                 zr0=zr0*(a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c+k-1.0d0))*z1
                 zhf=zhf+zr0
                 if (cdabs(zhf-zw).lt.cdabs(zhf)*eps) go to 25
20               zw=zhf
25            zhf=zc0*zhf
           else if (a0.ge.0.90d0) then
              gm=0.0d0
              mcab=int(c-a-b+eps*dsign(1.0d0,c-a-b))
              if (dabs(c-a-b-mcab).lt.eps) then
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
30                  gm=gm*j
                 rm=1.0d0
                 do 35 j=1,abs(m)
35                  rm=rm*j
                 zf0=(1.0d0,0.0d0)
                 zr0=(1.0d0,0.0d0)
                 zr1=(1.0d0,0.0d0)
                 sp0=0.d0
                 sp=0.0d0
                 if (m.ge.0) then
                    zc0=gm*gc/(gam*gbm)
                    zc1=-gc*(z-1.0d0)**m/(ga*gb*rm)
                    do 40 k=1,m-1
                       zr0=zr0*(a+k-1.d0)*(b+k-1.d0)/(k*(k-m))*(1.d0-z)
40                     zf0=zf0+zr0
                    do 45 k=1,m
45                     sp0=sp0+1.0d0/(a+k-1.0d0)+1.0/(b+k-1.0d0)-1.d0/k
                    zf1=pa+pb+sp0+2.0d0*el+cdlog(1.0d0-z)
                    do 55 k=1,500
                       sp=sp+(1.0d0-a)/(k*(a+k-1.0d0))+(1.0d0-b)/
     &                    (k*(b+k-1.0d0))
                       sm=0.0d0
                       do 50 j=1,m
                          sm=sm+(1.0d0-a)/((j+k)*(a+j+k-1.0d0))
     &                       +1.0d0/(b+j+k-1.0d0)
50                     continue
                       zp=pa+pb+2.0d0*el+sp+sm+cdlog(1.0d0-z)
                       zr1=zr1*(a+m+k-1.0d0)*(b+m+k-1.0d0)/(k*(m+k))
     &                     *(1.0d0-z)
                       zf1=zf1+zr1*zp
                       if (cdabs(zf1-zw).lt.cdabs(zf1)*eps) go to 60
55                     zw=zf1
60                  zhf=zf0*zc0+zf1*zc1
                 else if (m.lt.0) then
                    m=-m
                    zc0=gm*gc/(ga*gb*(1.0d0-z)**m)
                    zc1=-(-1)**m*gc/(gam*gbm*rm)
                    do 65 k=1,m-1
                       zr0=zr0*(a-m+k-1.0d0)*(b-m+k-1.0d0)/(k*(k-m))
     &                     *(1.0d0-z)
65                     zf0=zf0+zr0
                    do 70 k=1,m
70                     sp0=sp0+1.0d0/k
                    zf1=pa+pb-sp0+2.0d0*el+cdlog(1.0d0-z)
                    do 80 k=1,500
                       sp=sp+(1.0d0-a)/(k*(a+k-1.0d0))+(1.0d0-b)/(k*
     &                    (b+k-1.0d0))
                       sm=0.0d0
                       do 75 j=1,m
75                        sm=sm+1.0d0/(j+k)
                       zp=pa+pb+2.0d0*el+sp-sm+cdlog(1.0d0-z)
                       zr1=zr1*(a+k-1.d0)*(b+k-1.d0)/(k*(m+k))*(1.d0-z)
                       zf1=zf1+zr1*zp
                       if (cdabs(zf1-zw).lt.cdabs(zf1)*eps) go to 85
80                     zw=zf1
85                  zhf=zf0*zc0+zf1*zc1
                 endif
              else
                 call gamma(a,ga)
                 call gamma(b,gb)
                 call gamma(c,gc)
                 call gamma(c-a,gca)
                 call gamma(c-b,gcb)
                 call gamma(c-a-b,gcab)
                 call gamma(a+b-c,gabc)
                 zc0=gc*gcab/(gca*gcb)
                 zc1=gc*gabc/(ga*gb)*(1.0d0-z)**(c-a-b)
                 zhf=(0.0d0,0.0d0)
                 zr0=zc0
                 zr1=zc1
                 do 90 k=1,500
                    zr0=zr0*(a+k-1.d0)*(b+k-1.d0)/(k*(a+b-c+k))*(1.d0-z)
                    zr1=zr1*(c-a+k-1.0d0)*(c-b+k-1.0d0)/(k*(c-a-b+k))
     &                  *(1.0d0-z)
                    zhf=zhf+zr0+zr1
                    if (cdabs(zhf-zw).lt.cdabs(zhf)*eps) go to 95
90                  zw=zhf
95               zhf=zhf+zc0+zc1
              endif
           else
              z00=(1.0d0,0.0d0)
              if (c-a.lt.a.and.c-b.lt.b) then
                  z00=(1.0d0-z)**(c-a-b)
                  a=c-a
                  b=c-b
              endif
              zhf=(1.0d0,0.d0)
              zr=(1.0d0,0.0d0)
              do 100 k=1,1500
                 zr=zr*(a+k-1.0d0)*(b+k-1.0d0)/(k*(c+k-1.0d0))*z
                 zhf=zhf+zr
                 if (cdabs(zhf-zw).le.cdabs(zhf)*eps) go to 105
100              zw=zhf
105           zhf=z00*zhf
           endif
        else if (a0.gt.1.0d0) then
           mab=int(a-b+eps*dsign(1.0d0,a-b))
           if (dabs(a-b-mab).lt.eps.and.a0.le.1.1d0) b=b+eps
           if (dabs(a-b-mab).gt.eps) then
              call gamma(a,ga)
              call gamma(b,gb)
              call gamma(c,gc)
              call gamma(a-b,gab)
              call gamma(b-a,gba)
              call gamma(c-a,gca)
              call gamma(c-b,gcb)
              zc0=gc*gba/(gca*gb*(-z)**a)
              zc1=gc*gab/(gcb*ga*(-z)**b)
              zr0=zc0
              zr1=zc1
              zhf=(0.0d0,0.0d0)
              do 110 k=1,500
                 zr0=zr0*(a+k-1.0d0)*(a-c+k)/((a-b+k)*k*z)
                 zr1=zr1*(b+k-1.0d0)*(b-c+k)/((b-a+k)*k*z)
                 zhf=zhf+zr0+zr1
                 if (cdabs((zhf-zw)/zhf).le.eps) go to 115
110              zw=zhf
115           zhf=zhf+zc0+zc1
           else
              if (a-b.lt.0.0d0) then
                 a=bb
                 b=aa
              endif
              ca=c-a
              cb=c-b
              nca=int(ca+eps*dsign(1.0d0,ca))
              ncb=int(cb+eps*dsign(1.0d0,cb))
              if (dabs(ca-nca).lt.eps.or.dabs(cb-ncb).lt.eps) c=c+eps
              call gamma(a,ga)
              call gamma(c,gc)
              call gamma(c-b,gcb)
              call psi(a,pa)
              call psi(c-a,pca)
              call psi(a-c,pac)
              mab=int(a-b+eps)
              zc0=gc/(ga*(-z)**b)
              call gamma(a-b,gm)
              zf0=gm/gcb*zc0
              zr=zc0
              do 120 k=1,mab-1
                 zr=zr*(b+k-1.0d0)/(k*z)
                 t0=a-b-k
                 call gamma(t0,g0)
                 call gamma(c-b-k,gcbk)
120              zf0=zf0+zr*g0/gcbk
              if (mab.eq.0) zf0=(0.0d0,0.0d0)
              zc1=gc/(ga*gcb*(-z)**a)
              sp=-2.0d0*el-pa-pca
              do 125 j=1,mab
125              sp=sp+1.0d0/j
              zp0=sp+cdlog(-z)
              sq=1.0d0
              do 130 j=1,mab
130              sq=sq*(b+j-1.0d0)*(b-c+j)/j
              zf1=(sq*zp0)*zc1
              zr=zc1
              rk1=1.0d0
              sj1=0.0d0
              do 145 k=1,10000
                 zr=zr/z
                 rk1=rk1*(b+k-1.0d0)*(b-c+k)/(k*k)
                 rk2=rk1
                 do 135 j=k+1,k+mab
135                 rk2=rk2*(b+j-1.0d0)*(b-c+j)/j
                 sj1=sj1+(a-1.0d0)/(k*(a+k-1.0d0))+(a-c-1.0d0)/
     &               (k*(a-c+k-1.0d0))
                 sj2=sj1
                 do 140 j=k+1,k+mab
140                 sj2=sj2+1.0d0/j
                 zp=-2.0d0*el-pa-pac+sj2-1.0d0/(k+a-c)
     &              -pi/dtan(pi*(k+a-c))+cdlog(-z)
                 zf1=zf1+rk2*zr*zp
                 ws=cdabs(zf1)
                 if (dabs((ws-w0)/ws).lt.eps) go to 150
145              w0=ws
150           zhf=zf0+zf1
           endif
        endif
155     a=aa
        b=bb
        if (k.gt.150) write(*,160)
160     format(1x,'warning! you should check the accuracy')
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
