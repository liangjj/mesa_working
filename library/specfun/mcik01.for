        program mcik01
c
c       =============================================================
c       purpose: this program computes the modified bessel functions  
c                i0(z), i1(z), k0(z), k1(z), and their derivatives 
c                for a complex argument using subroutine cik01
c       input :  z --- complex argument
c       output:  cbi0 --- i0(z)
c                cdi0 --- i0'(z)
c                cbi1 --- i1(z)
c                cdi1 --- i1'(z)
c                cbk0 --- k0(z)
c                cdk0 --- k0'(z)
c                cbk1 --- k1(z)
c                cdk1 --- k1'(z)
c       example: z = 20.0 + i 10.0
c
c     n      re[in(z)]      im[in(z)]      re[in'(z)]     im[in'(z)]
c    -----------------------------------------------------------------
c     0   -.38773811d+08 -.13750292d+08 -.37852037d+08 -.13869150d+08
c     1   -.37852037d+08 -.13869150d+08 -.36982347d+08 -.13952566d+08
c
c     n      re[kn(z)]      im[kn(z)]      re[kn'(z)]     im[kn'(z)]
c    -----------------------------------------------------------------
c     0   -.37692389d-09  .39171613d-09  .38056380d-09 -.40319029d-09
c     1   -.38056380d-09  .40319029d-09  .38408264d-09 -.41545502d-09
c       =============================================================
c
        implicit double precision (x,y)
        implicit complex*16 (c,z)
        write(*,*)'  please enter x and y (z=x+iy) '
        read(*,*)x,y
        z=cmplx(x,y)
        write(*,30)x,y
        call cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1)
        write(*,*)
        write(*,*)'  n      re[in(z)]      im[in(z)]',
     &            '      re[in''(z)]     im[in''(z)]'
        write(*,*)' -------------------------------',
     &            '----------------------------------'
        write(*,10)cbi0,cdi0
        write(*,20)cbi1,cdi1
        write(*,*)
        write(*,*)'  n      re[kn(z)]      im[kn(z)]',
     &            '      re[kn''(z)]     im[kn''(z)]'
        write(*,*)' -------------------------------',
     &            '----------------------------------'
        write(*,10)cbk0,cdk0
        write(*,20)cbk1,cdk1
10      format(3x,'0',2x,4d15.7)
20      format(3x,'1',2x,4d15.7)
30      format(3x,3hz =,f7.2,' + i',f7.2)
        end


        subroutine cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1)
c
c       ==========================================================
c       purpose: compute modified bessel functions i0(z), i1(z), 
c                k0(z), k1(z), and their derivatives for a 
c                complex argument
c       input :  z --- complex argument
c       output:  cbi0 --- i0(z)
c                cdi0 --- i0'(z)
c                cbi1 --- i1(z)
c                cdi1 --- i1'(z)
c                cbk0 --- k0(z)
c                cdk0 --- k0'(z)
c                cbk1 --- k1(z)
c                cdk1 --- k1'(z)
c       ==========================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension a(12),b(12),a1(10)
        pi=3.141592653589793d0
        ci=(0.0d0,1.0d0)
        a0=cdabs(z)
        z2=z*z
        z1=z
        if (a0.eq.0.0d0) then
           cbi0=(1.0d0,0.0d0)
           cbi1=(0.0d0,0.0d0)
           cdi0=(0.0d0,0.0d0)
           cdi1=(0.5d0,0.0d0)
           cbk0=(1.0d+300,0.0d0)
           cbk1=(1.0d+300,0.0d0)
           cdk0=-(1.0d+300,0.0d0)
           cdk1=-(1.0d+300,0.0d0)
           return
        endif
        if (real(z).lt.0.0) z1=-z
        if (a0.le.18.0) then
           cbi0=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 10 k=1,50
              cr=0.25d0*cr*z2/(k*k)
              cbi0=cbi0+cr
              if (cdabs(cr/cbi0).lt.1.0d-15) go to 15
10         continue
15         cbi1=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 20 k=1,50
              cr=0.25d0*cr*z2/(k*(k+1))
              cbi1=cbi1+cr
              if (cdabs(cr/cbi1).lt.1.0d-15) go to 25
20         continue
25         cbi1=0.5d0*z1*cbi1
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
           if (a0.ge.35.0) k0=9
           if (a0.ge.50.0) k0=7
           ca=cdexp(z1)/cdsqrt(2.0d0*pi*z1)
           cbi0=(1.0d0,0.0d0)
           zr=1.0d0/z1
           do 30 k=1,k0
30            cbi0=cbi0+a(k)*zr**k
           cbi0=ca*cbi0
           cbi1=(1.0d0,0.0d0)
           do 35 k=1,k0
35            cbi1=cbi1+b(k)*zr**k
           cbi1=ca*cbi1
        endif
        if (a0.le.9.0) then
           cs=(0.0d0,0.0d0)
           ct=-cdlog(0.5d0*z1)-0.5772156649015329d0
           w0=0.0d0
           cr=(1.0d0,0.0d0)
           do 40 k=1,50
              w0=w0+1.0d0/k
              cr=0.25d0*cr/(k*k)*z2
              cs=cs+cr*(w0+ct)
              if (cdabs((cs-cw)/cs).lt.1.0d-15) go to 45
40            cw=cs
45         cbk0=ct+cs
        else
           data a1/0.125d0,0.2109375d0,
     &             1.0986328125d0,1.1775970458984d01,
     &             2.1461706161499d02,5.9511522710323d03,
     &             2.3347645606175d05,1.2312234987631d07,
     &             8.401390346421d08,7.2031420482627d10/
           cb=0.5d0/z1
           zr2=1.0d0/z2
           cbk0=(1.0d0,0.0d0)
           do 50 k=1,10
50            cbk0=cbk0+a1(k)*zr2**k
           cbk0=cb*cbk0/cbi0
        endif
        cbk1=(1.0d0/z1-cbi1*cbk0)/cbi0
        if (real(z).lt.0.0) then
           if (dimag(z).lt.0.0) cbk0=cbk0+ci*pi*cbi0
           if (dimag(z).gt.0.0) cbk0=cbk0-ci*pi*cbi0
           if (dimag(z).lt.0.0) cbk1=-cbk1+ci*pi*cbi1
           if (dimag(z).gt.0.0) cbk1=-cbk1-ci*pi*cbi1
           cbi1=-cbi1
        endif
        cdi0=cbi1
        cdi1=cbi0-1.0d0/z*cbi1
        cdk0=-cbk1
        cdk1=-cbk0-1.0d0/z*cbk1
        return
        end
