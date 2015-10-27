        program mcjyna
c
c       ================================================================
c       purpose: this program computes bessel functions jn(z), yn(z)  
c                and their derivatives for a complex argument using
c                subroutine cjyna
c       input :  z --- complex argument of jn(z) and yn(z)
c                n --- order of jn(z) and yn(z)
c                      ( n = 0,1,תתת, n ף 250 )
c       output:  cbj(n) --- jn(z)
c                cdj(n) --- jn'(z)
c                cby(n) --- yn(z)
c                cdy(n) --- yn'(z)
c       eaxmple: z = 4.0 + i 2.0
c
c     n     re[jn(z)]       im[jn(z)]       re[jn'(z)]      im[jn'(z)]
c    -------------------------------------------------------------------
c     0  -.13787022d+01   .39054236d+00   .50735255d+00   .12263041d+01
c     1  -.50735255d+00  -.12263041d+01  -.11546013d+01   .58506793d+00
c     2   .93050039d+00  -.77959350d+00  -.72363400d+00  -.72836666d+00
c     3   .93991546d+00   .23042918d+00   .29742236d+00  -.63587637d+00
c     4   .33565567d+00   .49215925d+00   .47452722d+00  -.29035945d-01
c     5  -.91389835d-02   .28850107d+00   .20054412d+00   .19908868d+00
c
c     n     re[yn(z)]       im[yn(z)]       re[yn'(z)]      im[yn'(z)]
c   --------------------------------------------------------------------
c     0  -.38145893d+00  -.13291649d+01  -.12793101d+01   .51220420d+00
c     1   .12793101d+01  -.51220420d+00  -.58610052d+00  -.10987930d+01
c     2   .79074211d+00   .86842120d+00   .78932897d+00  -.70142425d+00
c     3  -.29934789d+00   .89064431d+00   .70315755d+00   .24423024d+00
c     4  -.61557299d+00   .37996071d+00   .41126221d-01   .34044655d+00
c     5  -.38160033d+00   .20975121d+00  -.33884827d+00  -.20590670d-01
c
c                z = 20.0 + i 10.0 ,      nmax = 5
c
c     n     re[jn(z)]       im[jn(z)]       re[jn'(z)]      im[jn'(z)]
c   --------------------------------------------------------------------
c     0   .15460268d+04  -.10391216d+04  -.10601232d+04  -.15098284d+04
c     1   .10601232d+04   .15098284d+04   .14734253d+04  -.10783122d+04
c     2  -.14008238d+04   .11175029d+04   .11274890d+04   .13643952d+04
c     3  -.11948548d+04  -.12189620d+04  -.11843035d+04   .11920871d+04
c     4   .96778325d+03  -.12666712d+04  -.12483664d+04  -.93887194d+03
c     5   .13018781d+04   .65878188d+03   .64152944d+03  -.12682398d+04
c
c     n     re[yn(z)]       im[yn(z)]       re[yn'(z)]      im[yn'(z)]
c   --------------------------------------------------------------------
c     0   .10391216d+04   .15460268d+04   .15098284d+04  -.10601232d+04
c     1  -.15098284d+04   .10601232d+04   .10783122d+04   .14734253d+04
c     2  -.11175029d+04  -.14008238d+04  -.13643952d+04   .11274890d+04
c     3   .12189620d+04  -.11948548d+04  -.11920871d+04  -.11843035d+04
c     4   .12666712d+04   .96778324d+03   .93887194d+03  -.12483664d+04
c     5  -.65878189d+03   .13018781d+04   .12682398d+04   .64152944d+03
c       ================================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        common cbj(0:251),cdj(0:251),cby(0:251),cdy(0:251)
        write(*,*)'  please enter n, x,y (z=x+iy) '
        read(*,*)n,x,y
        z=cmplx(x,y)
        write(*,40)x,y,n
        if (n.le.8) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call cjyna(n,z,nm,cbj,cdj,cby,cdy)
        write(*,*)
        write(*,*)'   n     re[jn(z)]       im[jn(z)]',
     &            '       re[jn''(z)]      im[jn''(z)]'
        write(*,*)' -------------------------------------',
     &            '-------------------------------'
        do 10 k=0,nm,ns
10         write(*,30)k,cbj(k),cdj(k)
        write(*,*)
        write(*,*)'   n     re[yn(z)]       im[yn(z)]',
     &            '       re[yn''(z)]      im[yn''(z)]'
        write(*,*)' -------------------------------------',
     &            '-------------------------------'
        do 20 k=0,nm,ns
20         write(*,30)k,cby(k),cdy(k)
30      format(1x,i4,4d16.8)
40      format(3x,3hz =,f5.1,' + i ',f5.1,' ,',6x,6hnmax =,i3)
        end


        subroutine cjyna(n,z,nm,cbj,cdj,cby,cdy)
c
c       =======================================================
c       purpose: compute bessel functions jn(z), yn(z) and
c                their derivatives for a complex argument
c       input :  z --- complex argument of jn(z) and yn(z)
c                n --- order of jn(z) and yn(z)
c       output:  cbj(n) --- jn(z)
c                cdj(n) --- jn'(z)
c                cby(n) --- yn(z)
c                cdy(n) --- yn'(z)
c                nm --- highest order computed
c       rouitines called:
c            (1) cjy01 to calculate j0(z), j1(z), y0(z), y1(z)
c            (2) msta1 and msta2 to calculate the starting 
c                point for backward recurrence
c       =======================================================
c
        implicit double precision (a,b,e,p,r,w,y)
        implicit complex*16 (c,z)
        dimension cbj(0:n),cdj(0:n),cby(0:n),cdy(0:n)
        pi=3.141592653589793d0
        a0=cdabs(z)
        nm=n
        if (a0.lt.1.0d-100) then
           do 5 k=0,n
              cbj(k)=(0.0d0,0.0d0)
              cdj(k)=(0.0d0,0.0d0)
              cby(k)=-(1.0d+300,0.0d0)
5             cdy(k)=(1.0d+300,0.0d0)
           cbj(0)=(1.0d0,0.0d0)
           cdj(1)=(0.5d0,0.0d0)
           return
        endif
        call cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1)
        cbj(0)=cbj0
        cbj(1)=cbj1
        cby(0)=cby0
        cby(1)=cby1
        cdj(0)=cdj0
        cdj(1)=cdj1
        cdy(0)=cdy0
        cdy(1)=cdy1
        if (n.le.1) return
        if (n.lt.int(0.25*a0)) then
           cj0=cbj0
           cj1=cbj1
           do 70 k=2,n
              cjk=2.0d0*(k-1.0d0)/z*cj1-cj0
              cbj(k)=cjk
              cj0=cj1
70            cj1=cjk
        else
           m=msta1(a0,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(a0,n,15)
           endif
           cf2=(0.0d0,0.0d0)
           cf1=(1.0d-100,0.0d0)
           do 75 k=m,0,-1
              cf=2.0d0*(k+1.0d0)/z*cf1-cf2
              if (k.le.nm) cbj(k)=cf
              cf2=cf1
75            cf1=cf
           if (cdabs(cbj0).gt.cdabs(cbj1)) then
              cs=cbj0/cf
           else
              cs=cbj1/cf2
           endif
           do 80 k=0,nm
80            cbj(k)=cs*cbj(k)
        endif
        do 85 k=2,nm
85         cdj(k)=cbj(k-1)-k/z*cbj(k)
        ya0=cdabs(cby0)
        lb=0
        cg0=cby0
        cg1=cby1
        do 90 k=2,nm
           cyk=2.0d0*(k-1.0d0)/z*cg1-cg0
           if (cdabs(cyk).gt.1.0d+290) go to 90            
           yak=cdabs(cyk)
           ya1=cdabs(cg0)
           if (yak.lt.ya0.and.yak.lt.ya1) lb=k
           cby(k)=cyk
           cg0=cg1
           cg1=cyk
90      continue
        if (lb.le.4.or.dimag(z).eq.0.0d0) go to 125
95      if (lb.eq.lb0) go to 125
        ch2=(1.0d0,0.0d0)
        ch1=(0.0d0,0.0d0)
        lb0=lb
        do 100 k=lb,1,-1
           ch0=2.0d0*k/z*ch1-ch2
           ch2=ch1
100        ch1=ch0
        cp12=ch0
        cp22=ch2
        ch2=(0.0d0,0.0d0)
        ch1=(1.0d0,0.0d0)
        do 105 k=lb,1,-1
           ch0=2.0d0*k/z*ch1-ch2
           ch2=ch1
105        ch1=ch0
        cp11=ch0
        cp21=ch2
        if (lb.eq.nm) cbj(lb+1)=2.0d0*lb/z*cbj(lb)-cbj(lb-1)
        if (cdabs(cbj(0)).gt.cdabs(cbj(1))) then
           cby(lb+1)=(cbj(lb+1)*cby0-2.0d0*cp11/(pi*z))/cbj(0)
           cby(lb)=(cbj(lb)*cby0+2.0d0*cp12/(pi*z))/cbj(0)
        else
           cby(lb+1)=(cbj(lb+1)*cby1-2.0d0*cp21/(pi*z))/cbj(1)
           cby(lb)=(cbj(lb)*cby1+2.0d0*cp22/(pi*z))/cbj(1)
        endif
        cyl2=cby(lb+1)
        cyl1=cby(lb)
        do 110 k=lb-1,0,-1
           cylk=2.0d0*(k+1.0d0)/z*cyl1-cyl2
           cby(k)=cylk
           cyl2=cyl1
110        cyl1=cylk
        cyl1=cby(lb)
        cyl2=cby(lb+1)
        do 115 k=lb+1,nm-1
           cylk=2.0d0*k/z*cyl2-cyl1
           cby(k+1)=cylk
           cyl1=cyl2
115        cyl2=cylk
        do 120 k=2,nm
           wa=cdabs(cby(k))
           if (wa.lt.cdabs(cby(k-1))) lb=k
120     continue
        go to 95
125     continue
        do 130 k=2,nm
130        cdy(k)=cby(k-1)-k/z*cby(k)
        return
        end


        subroutine cjy01(z,cbj0,cdj0,cbj1,cdj1,cby0,cdy0,cby1,cdy1)
c
c       ===========================================================
c       purpose: compute complex bessel functions j0(z), j1(z)
c                y0(z), y1(z), and their derivatives
c       input :  z --- complex argument
c       output:  cbj0 --- j0(z)
c                cdj0 --- j0'(z)
c                cbj1 --- j1(z)
c                cdj1 --- j1'(z)
c                cby0 --- y0(z)
c                cdy0 --- y0'(z)
c                cby1 --- y1(z)
c                cdy1 --- y1'(z)
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
              if (cdabs(cr/cbj0).lt.1.0d-15) go to 15
10         continue
15         cbj1=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 20 k=1,40
              cr=-0.25d0*cr*z2/(k*(k+1.0d0))
              cbj1=cbj1+cr
              if (cdabs(cr/cbj1).lt.1.0d-15) go to 25
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
              if (cdabs(cp/cs).lt.1.0d-15) go to 35
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
              if (cdabs(cp/cs).lt.1.0d-15) go to 45
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
           ct1=z1-0.25d0*pi
           cp0=(1.0d0,0.0d0)
           do 50 k=1,k0
50            cp0=cp0+a(k)*z1**(-2*k)
           cq0=-0.125d0/z1
           do 55 k=1,k0
55            cq0=cq0+b(k)*z1**(-2*k-1)
           cu=cdsqrt(rp2/z1)
           cbj0=cu*(cp0*cdcos(ct1)-cq0*cdsin(ct1))
           cby0=cu*(cp0*cdsin(ct1)+cq0*cdcos(ct1))
           ct2=z1-0.75d0*pi
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
