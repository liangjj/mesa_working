        program mcjy01
c
c       ================================================================
c       purpose: this program computes bessel functions j0(z), j1(z),
c                y0(z), y1(z), and their derivatives for a complex
c                argument using subroutine cjy01
c       input :  z --- complex argument
c       output:  cbj0 --- j0(z)
c                cdj0 --- j0'(z)
c                cbj1 --- j1(z)
c                cdj1 --- j1'(z)
c                cby0 --- y0(z)
c                cdy0 --- y0'(z)
c                cby1 --- y1(z)
c                cdy1 --- y1'(z)
c       example: z =  4.0 + i  2.0
c
c     n     re[jn(z)]       im[jn(z)]       re[jn'(z)]      im[jn'(z)]
c   --------------------------------------------------------------------
c     0  -.13787022d+01   .39054236d+00   .50735255d+00   .12263041d+01
c     1  -.50735255d+00  -.12263041d+01  -.11546013d+01   .58506793d+00
c
c     n     re[yn(z)]       im[yn(z)]       re[yn'(z)]      im[yn'(z)]
c   --------------------------------------------------------------------
c     0  -.38145893d+00  -.13291649d+01  -.12793101d+01   .51220420d+00
c     1   .12793101d+01  -.51220420d+00  -.58610052d+00  -.10987930d+01
c       ================================================================
c
        implicit double precision (x,y)
        implicit complex*16 (c,z)
        write(*,*)'  please enter x,y (z=x+iy) '
        read(*,*)x,y
        z=cmplx(x,y)
        write(*,30)x,y
        call cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1)
        write(*,*)
        write(*,*)'  n      re[jn(z)]       im[jn(z)]',
     &            '       re[jn''(z)]      im[jn''(z)]'
        write(*,*)' -------------------------------------',
     &            '-------------------------------'
        write(*,10)cbj0,cdj0
        write(*,20)cbj1,cdj1
        write(*,*)
        write(*,*)'  n      re[yn(z)]       im[yn(z)]',
     &            '       re[yn''(z)]      im[yn''(z)]'
        write(*,*)' -------------------------------------',
     &            '-------------------------------'
        write(*,10)cby0,cdy0
        write(*,20)cby1,cdy1
10      format(3x,'0',2x,4d16.8)
20      format(3x,'1',2x,4d16.8)
30      format(3x,3hz =,f6.2,' + i',f6.2)
        end


        subroutine cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1)
c
c       =======================================================
c       purpose: compute bessel functions j0(z), j1(z), y0(z), 
c                y1(z), and their derivatives for a complex
c                argument
c       input :  z --- complex argument
c       output:  cbj0 --- j0(z)
c                cdj0 --- j0'(z)
c                cbj1 --- j1(z)
c                cdj1 --- j1'(z)
c                cby0 --- y0(z)
c                cdy0 --- y0'(z)
c                cby1 --- y1(z)
c                cdy1 --- y1'(z)
c       =======================================================
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
           cdj0=(0.0d0,0.0d0)
           cdj1=(0.5d0,0.0d0)
           cby0=-(1.0d300,0.0d0)
           cby1=-(1.0d300,0.0d0)
           cdy0=(1.0d300,0.0d0)
           cdy1=(1.0d300,0.0d0)
           return
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
        cdj0=-cbj1
        cdj1=cbj0-1.0d0/z*cbj1
        cdy0=-cby1
        cdy1=cby0-1.0d0/z*cby1
        return
        end
