        program mcikna
c
c       =============================================================
c       purpose: this program computes the modified bessel functions 
c                in(z) and kn(z), and their derivatives for a  
c                complex argument using subroutine cikna
c       input :  z --- complex argument of in(z) and kn(z)
c                n --- order of in(z) and kn(z)
c                      ( n = 0,1,תתת, n ף 250 )
c       output:  cbi(n) --- in(z)
c                cdi(n) --- in'(z)
c                cbk(n) --- kn(z)
c                cdk(n) --- kn'(z)
c       example: z = 4.0 + i 2.0 ,      nmax = 5
c
c     n     re[in(z)]      im[in(z)]      re[in'(z)]     im[in'(z)]
c   -----------------------------------------------------------------
c     0  -.19056142d+01  .10403505d+02 -.23059657d+01  .92222463d+01
c     1  -.23059657d+01  .92222463d+01 -.23666457d+01  .83284588d+01
c     2  -.28276772d+01  .62534130d+01 -.24255774d+01  .61553456d+01
c     3  -.25451891d+01  .30884450d+01 -.22270972d+01  .36367893d+01
c     4  -.16265172d+01  .10201656d+01 -.16520416d+01  .16217056d+01
c     5  -.75889410d+00  .15496632d+00 -.94510625d+00  .48575220d+00
c
c     n     re[kn(z)]      im[kn(z)]      re[kn'(z)]     im[kn'(z)]
c   -----------------------------------------------------------------
c     0  -.64221754d-02 -.84393648d-02  .74307276d-02  .89585853d-02
c     1  -.74307276d-02 -.89585853d-02  .88041795d-02  .94880091d-02
c     2  -.11186184d-01 -.10536653d-01  .14012532d-01  .10936010d-01
c     3  -.20594336d-01 -.12913435d-01  .27416815d-01  .12106413d-01
c     4  -.43647447d-01 -.13676173d-01  .60982763d-01  .63953943d-02
c     5  -.10137119d+00  .12264588d-03  .14495731d+00 -.37132068d-01
c       =============================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        common cbi(0:250),cdi(0:250),cbk(0:250),cdk(0:250)
        write(*,*)'  please input n, x,y (z=x+iy)=?'
        read(*,*)n,x,y
        z=cmplx(x,y)
        write(*,40)x,y,n
        if (n.le.8) then
           ns=1
        else
           write(*,*)' please enter order step ns'
           read(*,*)ns
        endif
        call cikna(n,z,nm,cbi,cdi,cbk,cdk)
        write(*,*)
        write(*,*)'   n      re[in(z)]       im[in(z)]',
     &            '       re[in''(z)]      im[in''(z)]'
        write(*,*)' -----------------------------------',
     &            '----------------------------------'
        do 10 k=0,nm,ns
           write(*,30)k,cbi(k),cdi(k)
10      continue
        write(*,*)
        write(*,*)'   n      re[kn(z)]       im[kn(z)]',
     &            '       re[kn''(z)]      im[kn''(z)]'
        write(*,*)' -----------------------------------',
     &            '----------------------------------'
        do 20 k=0,nm,ns
           write(*,30)k,cbk(k),cdk(k)
20      continue
30      format(1x,i4,1x,4d16.8)
40      format(3x,3hz =,f7.1,' + i',f7.1,' ,',6x,6hnmax =,i4)
        end


        subroutine cikna(n,z,nm,cbi,cdi,cbk,cdk)
c
c       ========================================================
c       purpose: compute modified bessel functions in(z), kn(x)
c                and their derivatives for a complex argument
c       input :  z --- complex argument of in(z) and kn(z)
c                n --- order of in(z) and kn(z)
c       output:  cbi(n) --- in(z)
c                cdi(n) --- in'(z)
c                cbk(n) --- kn(z)
c                cdk(n) --- kn'(z)
c                nm --- highest order computed
c       routines called:
c             (1) cik01 to compute i0(z), i1(z) k0(z) & k1(z)
c             (2) msta1 and msta2 to compute the starting
c                 point for backward recurrence
c       ========================================================
c
        implicit double precision (a,b,p,w,x,y)
        implicit complex*16 (c,z)
        dimension cbi(0:n),cdi(0:n),cbk(0:n),cdk(0:n)
        a0=cdabs(z)
        nm=n
        if (a0.lt.1.0d-100) then
           do 10 k=0,n
              cbi(k)=(0.0d0,0.0d0)
              cdi(k)=(0.0d0,0.0d0)
              cbk(k)=-(1.0d+300,0.0d0)
10            cdk(k)=(1.0d+300,0.0d0)
           cbi(0)=(1.0d0,0.0d0)
           cdi(1)=(0.5d0,0.0d0)
           return
        endif
        call cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1)
        cbi(0)=cbi0
        cbi(1)=cbi1
        cbk(0)=cbk0
        cbk(1)=cbk1
        cdi(0)=cdi0
        cdi(1)=cdi1
        cdk(0)=cdk0
        cdk(1)=cdk1
        if (n.le.1) return
        m=msta1(a0,200)
        if (m.lt.n) then
           nm=m
        else
           m=msta2(a0,n,15)
        endif
        cf2=(0.0d0,0.0d0)
        cf1=(1.0d-100,0.0d0)
        do 45 k=m,0,-1
           cf=2.0d0*(k+1.0d0)/z*cf1+cf2
           if (k.le.nm) cbi(k)=cf
           cf2=cf1
45         cf1=cf
        cs=cbi0/cf
        do 50 k=0,nm
50         cbi(k)=cs*cbi(k)
        do 60 k=2,nm
           if (cdabs(cbi(k-1)).gt.cdabs(cbi(k-2))) then
              ckk=(1.0d0/z-cbi(k)*cbk(k-1))/cbi(k-1)
           else
              ckk=(cbi(k)*cbk(k-2)+2.0d0*(k-1.0d0)/(z*z))/cbi(k-2)
           endif
60         cbk(k)=ckk
        do 70 k=2,nm
           cdi(k)=cbi(k-1)-k/z*cbi(k)
70         cdk(k)=-cbk(k-1)-k/z*cbk(k)
        return
        end


        subroutine cik01(z,cbi0,cdi0,cbi1,cdi1,cbk0,cdk0,cbk1,cdk1)
c
c       ==========================================================
c       purpose: compute modified complex bessel functions i0(z),
c                i1(z), k0(z), k1(z), and their derivatives
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
