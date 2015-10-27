        program mcyzo
c
c       ===========================================================
c       purpose : this program evaluates the complex zeros of 
c                 y0(z), y0'(z), y1(z) and y1'(z), and their 
c                 associated values at the zeros using the
c                 modified newton's iteration method
c       input:    nt --- total number of roots/zeros
c                 kf --- function choice code
c                        kf=0 for  y0(z) & y1(z0)
c                        kf=1 for  y1(z) & y0(z1)
c                        kf=2 for  y1'(z) & y1(z1')
c                 kc --- choice code
c                        kc=0 for complex roots
c                        kc=1 for real roots
c       output:   zo(l) --- l-th zero of y0(z) or y1(z) or y1'(z)
c                 zv(l) --- value of y0'(z) or y1'(z) or y1(z)
c                           at the l-th zero
c       examples: nt = 5
c
c   no.      z0, zeros of y0(z)                y1(z0)
c  -----------------------------------------------------------------
c    1   -2.403016632 + i .5398823130   .1007476893 - i .8819677101
c    2   -5.519876702 + i .5471800106  -.0292464182 + i .5871695027
c    3   -8.653672403 + i .5484120673   .0149080637 - i .4694587524
c    4  -11.791512030 + i .5488191184  -.0093736817 + i .4023045429
c    5  -14.930906564 + i .5490008289   .0065788031 - i .3575673214
c
c   no.      z1, zeros of y1(z)                 y0(z1)
c  -----------------------------------------------------------------
c    1    -.502743273 + i .7862437145  -.4595276847 + i1.3171019361
c    2   -3.833535193 + i .5623565382   .0483019087 - i .6925128842
c    3   -7.015903683 + i .5533930459  -.0201269494 + i .5186425332
c    4  -10.173573834 + i .5512733877   .0116140017 - i .4320329636
c    5  -13.323739307 + i .5504585830  -.0077719300 + i .3779698048
c
c   no.      z1', zeros of y1'(z)                y1(z1')
c   ----------------------------------------------------------------
c    1     .576785129 + i .9039847922  -.7634970879 + i .5892448647
c    2   -1.940477342 - i .7211859189   .1620640057 + i .9520278864
c    3   -5.333478617 - i .5672196368  -.0317940081 - i .5968536736
c    4   -8.536768577 - i .5560607040   .0154177166 + i .4726011652
c    5  -11.706175219 - i .5528590607  -.0095443768 - i .4037533396
c       ============================================================
c
        implicit complex *16 (z)
        dimension zo(50),zv(50)
        write(*,*)'please enter nt, kf and kc'
        write(*,*)'  nt --- total number of the roots'
        write(*,*)'  kf  --- function choice code'
        write(*,*)'          kf=0 for y0(z) & y1(z0)'
        write(*,*)'          kf=1 for y1(z) & y0(z1)'
        write(*,*)'          kf=2 for y1''(z) & y1(z1'')'
        write(*,*)'  kc  --- choice code'
        write(*,*)'          kc=0 for complex roots'
        write(*,*)'          kc=1 for real roots'
        read(*,*)nt,kf,kc
        write(*,20)nt,kf,kc
        write(*,*)
        write(*,15)
        call cyzo(nt,kf,kc,zo,zv)
        write(*,*)
        if (kf.eq.0) then
           write(*,*)' no.          z0, zeros of y0(z)',
     &               '                 y1(z0)'
        else if (kf.eq.1) then
           write(*,*)' no.          z1, zeros of y1(z)',
     &               '                 y0(z1)'
        else if (kf.eq.2) then
           write(*,*)' no.        z1'', zeros of y1''(z)',
     &               '                y1(z1'')'
        endif
        write(*,*)'--------------------------------------',
     &            '----------------------------'
        do 10 i=1,nt
10         write(*,25)i,zo(i),zv(i)
15      format(20x,'*****    please wait !    *****')
20      format(2x,'nt=',i3,',  ','kf=',i3,',  ','kc=',i3)
25      format(1x,i3,2x,f15.9,3f15.10)
        end


        subroutine cyzo(nt,kf,kc,zo,zv)
c
c       ===========================================================
c       purpose : compute the complex zeros of y0(z), y1(z) and
c                 y1'(z), and their associated values at the zeros 
c                 using the modified newton's iteration method
c       input:    nt --- total number of zeros/roots
c                 kf --- function choice code
c                        kf=0 for  y0(z) & y1(z0)
c                        kf=1 for  y1(z) & y0(z1)
c                        kf=2 for  y1'(z) & y1(z1')
c                 kc --- choice code
c                        kc=0 for complex roots
c                        kc=1 for real roots
c       output:   zo(l) --- l-th zero of y0(z) or y1(z) or y1'(z)
c                 zv(l) --- value of y0'(z) or y1'(z) or y1(z)
c                           at the l-th zero
c       routine called: cy01 for computing y0(z) and y1(z), and
c                       their derivatives
c       ===========================================================
        implicit double precision (h,o-y)
        implicit complex*16 (c,z)
        dimension zo(nt),zv(nt)
        if (kc.eq.0) then
           x=-2.4d0
           y=0.54d0
           h=3.14
        else if (kc.eq.1) then
           x=0.89
           y=0.0
           h=-3.14
        endif
        if (kf.eq.1) x=-0.503
        if (kf.eq.2) x=0.577
        zero=cmplx(x,y)
        z=zero
        do 35 nr=1,nt
10         if (nr.ne.1) z=zo(nr-1)-h
           it=0
15         it=it+1
           call cy01(kf,z,zf,zd)
           zp=(1.0d0,0.0d0)
           do 20 i=1,nr-1
20            zp=zp*(z-zo(i))
           zfd=zf/zp
           zq=(0.0d0,0.0d0)
           do 30 i=1,nr-1
              zw=(1.0d0,0.0d0)
              do 25 j=1,nr-1
                 if (j.eq.i) go to 25
                 zw=zw*(z-zo(j))
25            continue
              zq=zq+zw
30         continue
           zgd=(zd-zq*zfd)/zp
           z=z-zfd/zgd
           w0=w
           w=cdabs(z)
           if (it.le.50.and.dabs((w-w0)/w).gt.1.0d-12) go to 15
           zo(nr)=z
35      continue
        do 40 i=1,nt
           z=zo(i)
           if (kf.eq.0.or.kf.eq.2) then
              call cy01(1,z,zf,zd)
              zv(i)=zf
           else if (kf.eq.1) then
              call cy01(0,z,zf,zd)
              zv(i)=zf
           endif
40      continue
        return
        end


        subroutine cy01(kf,z,zf,zd)
c
c       ===========================================================
c       purpose: compute complex bessel functions y0(z), y1(z)
c                and their derivatives
c       input :  z  --- complex argument of yn(z) ( n=0,1 )
c                kf --- function choice code
c                    kf=0 for zf=y0(z) and zd=y0'(z)
c                    kf=1 for zf=y1(z) and zd=y1'(z)
c                    kf=2 for zf=y1'(z) and zd=y1''(z)
c       output:  zf --- y0(z) or y1(z) or y1'(z)
c                zd --- y0'(z) or y1'(z) or y1''(z)
c       ===========================================================
c
        implicit double precision (a,b,e,p,r,w)
        implicit complex*16 (c,z)
        dimension a(12),b(12),a1(12),b1(12)
        pi=3.141592653589793d0
        el=0.5772156649015329d0
        rp2=2.0d0/pi
        ci=(0.0d0,1.0d0)
        a0=cdabs(z)
        z2=z*z
        z1=z
        if (a0.eq.0.0d0) then
           cbj0=(1.0d0,0.0d0)
           cbj1=(0.0d0,0.0d0)
           cby0=-(1.0d300,0.0d0)
           cby1=-(1.0d300,0.0d0)
           cdy0=(1.0d300,0.0d0)
           cdy1=(1.0d300,0.0d0)
           go to 70
        endif
        if (real(z).lt.0.0) z1=-z
        if (a0.le.12.0) then
           cbj0=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 10 k=1,40
              cr=-0.25d0*cr*z2/(k*k)
              cbj0=cbj0+cr
              if (cdabs(cr).lt.cdabs(cbj0)*1.0d-15) go to 15
10         continue
15         cbj1=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 20 k=1,40
              cr=-0.25d0*cr*z2/(k*(k+1.0d0))
              cbj1=cbj1+cr
              if (cdabs(cr).lt.cdabs(cbj1)*1.0d-15) go to 25
20         continue
25         cbj1=0.5d0*z1*cbj1
           w0=0.0d0
           cr=(1.0d0,0.0d0)
           cs=(0.0d0,0.0d0)
           do 30 k=1,40
              w0=w0+1.0d0/k
              cr=-0.25d0*cr/(k*k)*z2
              cp=cr*w0
              cs=cs+cp
              if (cdabs(cp).lt.cdabs(cs)*1.0d-15) go to 35
30         continue
35         cby0=rp2*(cdlog(z1/2.0d0)+el)*cbj0-rp2*cs
           w1=0.0d0
           cr=(1.0d0,0.0d0)
           cs=(1.0d0,0.0d0)
           do 40 k=1,40
              w1=w1+1.0d0/k
              cr=-0.25d0*cr/(k*(k+1))*z2
              cp=cr*(2.0d0*w1+1.0d0/(k+1.0d0))
              cs=cs+cp
              if (cdabs(cp).lt.cdabs(cs)*1.0d-15) go to 45
40         continue
45         cby1=rp2*((cdlog(z1/2.0d0)+el)*cbj1-1.0d0/z1-.25d0*z1*cs)
        else
           data a/-.703125d-01,.112152099609375d+00,
     &            -.5725014209747314d+00,.6074042001273483d+01,
     &            -.1100171402692467d+03,.3038090510922384d+04,
     &            -.1188384262567832d+06,.6252951493434797d+07,
     &            -.4259392165047669d+09,.3646840080706556d+11,
     &            -.3833534661393944d+13,.4854014686852901d+15/
           data b/ .732421875d-01,-.2271080017089844d+00,
     &             .1727727502584457d+01,-.2438052969955606d+02,
     &             .5513358961220206d+03,-.1825775547429318d+05,
     &             .8328593040162893d+06,-.5006958953198893d+08,
     &             .3836255180230433d+10,-.3649010818849833d+12,
     &             .4218971570284096d+14,-.5827244631566907d+16/
           data a1/.1171875d+00,-.144195556640625d+00,
     &             .6765925884246826d+00,-.6883914268109947d+01,
     &             .1215978918765359d+03,-.3302272294480852d+04,
     &             .1276412726461746d+06,-.6656367718817688d+07,
     &             .4502786003050393d+09,-.3833857520742790d+11,
     &             .4011838599133198d+13,-.5060568503314727d+15/
           data b1/-.1025390625d+00,.2775764465332031d+00,
     &             -.1993531733751297d+01,.2724882731126854d+02,
     &             -.6038440767050702d+03,.1971837591223663d+05,
     &             -.8902978767070678d+06,.5310411010968522d+08,
     &             -.4043620325107754d+10,.3827011346598605d+12,
     &             -.4406481417852278d+14,.6065091351222699d+16/
           k0=12
           if (a0.ge.35.0) k0=10
           if (a0.ge.50.0) k0=8
           ct1=z1-.25d0*pi
           cp0=(1.0d0,0.0d0)
           do 50 k=1,k0
50            cp0=cp0+a(k)*z1**(-2*k)
           cq0=-0.125d0/z1
           do 55 k=1,k0
55            cq0=cq0+b(k)*z1**(-2*k-1)
           cu=cdsqrt(rp2/z1)
           cbj0=cu*(cp0*cdcos(ct1)-cq0*cdsin(ct1))
           cby0=cu*(cp0*cdsin(ct1)+cq0*cdcos(ct1))
           ct2=z1-.75d0*pi
           cp1=(1.0d0,0.0d0)
           do 60 k=1,k0
60            cp1=cp1+a1(k)*z1**(-2*k)
           cq1=0.375d0/z1
           do 65 k=1,k0
65            cq1=cq1+b1(k)*z1**(-2*k-1)
           cbj1=cu*(cp1*cdcos(ct2)-cq1*cdsin(ct2))
           cby1=cu*(cp1*cdsin(ct2)+cq1*cdcos(ct2))
        endif
        if (real(z).lt.0.0) then
           if (dimag(z).lt.0.0) cby0=cby0-2.0d0*ci*cbj0
           if (dimag(z).gt.0.0) cby0=cby0+2.0d0*ci*cbj0
           if (dimag(z).lt.0.0) cby1=-(cby1-2.0d0*ci*cbj1)
           if (dimag(z).gt.0.0) cby1=-(cby1+2.0d0*ci*cbj1)
           cbj1=-cbj1
        endif
        cdy0=-cby1
        cdy1=cby0-1.0d0/z*cby1
70      if (kf.eq.0) then
           zf=cby0
           zd=cdy0
        else if (kf.eq.1) then
           zf=cby1
           zd=cdy1
        else if (kf.eq.2) then
           zf=cdy1
           zd=-cdy1/z-(1.0d0-1.0d0/(z*z))*cby1
        endif
        return
        end

