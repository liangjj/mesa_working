        program mjy01a
c
c       =========================================================
c       purpose: this program computes the bessel functions  
c                jn(x) and yn(x) ( n=0,1 ) and their derivatives 
c                using subroutine jy01a
c       input :  x   --- argument of jn(x) & yn(x) ( x � 0 )
c       output:  bj0 --- j0(x)
c                dj0 --- j0'(x)
c                bj1 --- j1(x)
c                dj1 --- j1'(x)
c                by0 --- y0(x)
c                dy0 --- y0'(x)
c                by1 --- y1(x)
c                dy1 --- y1'(x)
c       example:
c
c        x       j0(x)        j0'(x)       j1(x)        j1'(x)
c       ---------------------------------------------------------
c        1     .76519769   -.44005059    .44005059    .32514710
c        5    -.17759677    .32757914   -.32757914   -.11208094
c       10    -.24593576   -.04347275    .04347275   -.25028304
c       20     .16702466   -.06683312    .06683312    .16368301
c       30    -.08636798    .11875106   -.11875106   -.08240961
c       40     .00736689   -.12603832    .12603832    .00421593
c       50     .05581233    .09751183   -.09751183    .05776256
c
c        x       y0(x)        y0'(x)       y1(x)        y1'(x)
c      ---------------------------------------------------------
c        1     .08825696    .78121282   -.78121282    .86946979
c        5    -.30851763   -.14786314    .14786314   -.33809025
c       10     .05567117   -.24901542    .24901542    .03076962
c       20     .06264060    .16551161   -.16551161    .07091618
c       30    -.11729573   -.08442557    .08442557   -.12010992
c       40     .12593642    .00579351   -.00579351    .12608125
c       50    -.09806500    .05679567   -.05679567   -.09692908
c       =========================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter x'
        read(*,*)x
        write(*,20)x
        write(*,*)'  x          j0(x)          j0''(x)         j1(x)',
     &            '          j1''(x)'
        write(*,*)'------------------------------------------',
     &            '----------------------------'
        call jy01a(x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1)
        write(*,10)x,bj0,dj0,bj1,dj1
        write(*,*)
        write(*,*)'  x          y0(x)          y0''(x)         y1(x)',
     &            '          y1''(x)'
        write(*,*)'------------------------------------------',
     &            '----------------------------'
        write(*,10)x,by0,dy0,by1,dy1
10      format(1x,f5.1,4e16.8)
20      format(3x,'x =',f5.1)
        end


        subroutine jy01a(x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1)
c
c       =======================================================
c       purpose: compute bessel functions j0(x), j1(x), y0(x),
c                y1(x), and their derivatives
c       input :  x   --- argument of jn(x) & yn(x) ( x � 0 )
c       output:  bj0 --- j0(x)
c                dj0 --- j0'(x)
c                bj1 --- j1(x)
c                dj1 --- j1'(x)
c                by0 --- y0(x)
c                dy0 --- y0'(x)
c                by1 --- y1(x)
c                dy1 --- y1'(x)
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension a(12),b(12),a1(12),b1(12)
        pi=3.141592653589793d0
        rp2=0.63661977236758d0
        x2=x*x
        if (x.eq.0.0d0) then
           bj0=1.0d0
           bj1=0.0d0
           dj0=0.0d0
           dj1=0.5d0
           by0=-1.0d+300
           by1=-1.0d+300
           dy0=1.0d+300
           dy1=1.0d+300
           return
        endif
        if (x.le.12.0d0) then
           bj0=1.0d0
           r=1.0d0
           do 5 k=1,30
              r=-0.25d0*r*x2/(k*k)
              bj0=bj0+r
              if (dabs(r).lt.dabs(bj0)*1.0d-15) go to 10
5          continue
10         bj1=1.0d0
           r=1.0d0
           do 15 k=1,30
              r=-0.25d0*r*x2/(k*(k+1.0d0))
              bj1=bj1+r
              if (dabs(r).lt.dabs(bj1)*1.0d-15) go to 20
15         continue
20         bj1=0.5d0*x*bj1
           ec=dlog(x/2.0d0)+0.5772156649015329d0
           cs0=0.0d0
           w0=0.0d0
           r0=1.0d0
           do 25 k=1,30
              w0=w0+1.0d0/k
              r0=-0.25d0*r0/(k*k)*x2
              r=r0*w0
              cs0=cs0+r
              if (dabs(r).lt.dabs(cs0)*1.0d-15) go to 30
25         continue
30         by0=rp2*(ec*bj0-cs0)
           cs1=1.0d0
           w1=0.0d0
           r1=1.0d0
           do 35 k=1,30
              w1=w1+1.0d0/k
              r1=-0.25d0*r1/(k*(k+1))*x2
              r=r1*(2.0d0*w1+1.0d0/(k+1.0d0))
              cs1=cs1+r
              if (dabs(r).lt.dabs(cs1)*1.0d-15) go to 40
35         continue
40         by1=rp2*(ec*bj1-1.0d0/x-0.25d0*x*cs1)
        else
           data a/-.7031250000000000d-01,.1121520996093750d+00,
     &            -.5725014209747314d+00,.6074042001273483d+01,
     &            -.1100171402692467d+03,.3038090510922384d+04,
     &            -.1188384262567832d+06,.6252951493434797d+07,
     &            -.4259392165047669d+09,.3646840080706556d+11,
     &            -.3833534661393944d+13,.4854014686852901d+15/
           data b/ .7324218750000000d-01,-.2271080017089844d+00,
     &             .1727727502584457d+01,-.2438052969955606d+02,
     &             .5513358961220206d+03,-.1825775547429318d+05,
     &             .8328593040162893d+06,-.5006958953198893d+08,
     &             .3836255180230433d+10,-.3649010818849833d+12,
     &             .4218971570284096d+14,-.5827244631566907d+16/
           data a1/.1171875000000000d+00,-.1441955566406250d+00,
     &             .6765925884246826d+00,-.6883914268109947d+01,
     &             .1215978918765359d+03,-.3302272294480852d+04,
     &             .1276412726461746d+06,-.6656367718817688d+07,
     &             .4502786003050393d+09,-.3833857520742790d+11,
     &             .4011838599133198d+13,-.5060568503314727d+15/
           data b1/-.1025390625000000d+00,.2775764465332031d+00,
     &             -.1993531733751297d+01,.2724882731126854d+02,
     &             -.6038440767050702d+03,.1971837591223663d+05,
     &             -.8902978767070678d+06,.5310411010968522d+08,
     &             -.4043620325107754d+10,.3827011346598605d+12,
     &             -.4406481417852278d+14,.6065091351222699d+16/
           k0=12
           if (x.ge.35.0) k0=10
           if (x.ge.50.0) k0=8
           t1=x-0.25d0*pi
           p0=1.0d0
           q0=-0.125d0/x
           do 45 k=1,k0
              p0=p0+a(k)*x**(-2*k)
45            q0=q0+b(k)*x**(-2*k-1)
           cu=dsqrt(rp2/x)
           bj0=cu*(p0*dcos(t1)-q0*dsin(t1))
           by0=cu*(p0*dsin(t1)+q0*dcos(t1))
           t2=x-0.75d0*pi
           p1=1.0d0
           q1=0.375d0/x
           do 50 k=1,k0
              p1=p1+a1(k)*x**(-2*k)
50            q1=q1+b1(k)*x**(-2*k-1)
           cu=dsqrt(rp2/x)
           bj1=cu*(p1*dcos(t2)-q1*dsin(t2))
           by1=cu*(p1*dsin(t2)+q1*dcos(t2))
        endif
        dj0=-bj1
        dj1=bj0-bj1/x
        dy0=-by1
        dy1=by0-by1/x
        return
        end
