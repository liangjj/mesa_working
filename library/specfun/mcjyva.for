        program mcjyva
c
c       ===============================================================
c       purpose: this program computes bessel functions jv(z), yv(z),
c                and their derivatives for a complex argument using
c                subroutine cjyva
c       input :  z --- complex argument
c                v --- order of jv(z) and yv(z)
c                      ( v = n+v0, 0 ף n ף 250, 0 ף v0 < 1 )
c       output:  cbj(n) --- jn+v0(z)
c                cdj(n) --- jn+v0'(z)
c                cby(n) --- yn+v0(z)
c                cdy(n) --- yn+v0'(z)
c       example:
c                v = n +v0,  v0 = 1/3,   z = 4.0 + i 2.0
c
c     n     re[jv(z)]       im[jv(z)]      re[jv'(z)]      im[jv'(z)]
c    ------------------------------------------------------------------
c     0  -.13829878d+01  -.30855145d+00  -.18503756d+00   .13103689d+01
c     1   .82553327d-01  -.12848394d+01  -.12336901d+01   .45079506d-01
c     2   .10843924d+01  -.39871046d+00  -.33046401d+00  -.84574964d+00
c     3   .74348135d+00   .40665987d+00   .45318486d+00  -.42198992d+00
c     4   .17802266d+00   .44526939d+00   .39624497d+00   .97902890d-01
c     5  -.49008598d-01   .21085409d+00   .11784299d+00   .19422044d+00
c
c     n     re[yv(z)]      im[yv(z)]       re[yv'(z)]      im[yv'(z)]
c    ------------------------------------------------------------------
c     0   .34099851d+00  -.13440666d+01  -.13544477d+01  -.15470699d+00
c     1   .13323787d+01   .53735934d-01  -.21467271d-01  -.11807457d+01
c     2   .38393305d+00   .10174248d+01   .91581083d+00  -.33147794d+00
c     3  -.49924295d+00   .71669181d+00   .47786442d+00   .37321597d+00
c     4  -.57179578d+00   .27099289d+00  -.12111686d+00   .23405313d+00
c     5  -.25700924d+00   .24858555d+00  -.43023156d+00  -.13123662d+00
c       ===============================================================
c
        implicit double precision (v,x,y)
        implicit complex*16 (c,z)
        common cbj(0:251),cdj(0:251),cby(0:251),cdy(0:251)
        write(*,*)'  please enter v, x and y ( z=x+iy )'
        read(*,*)v,x,y
        z=cmplx(x,y)
        n=int(v)
        v0=v-n
        write(*,25)v0,x,y
        if (n.le.8) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call cjyva(v,z,vm,cbj,cdj,cby,cdy)
        nm=int(vm)
        write(*,*)
        write(*,*)'  n       re[jv(z)]       im[jv(z)]',
     &            '       re[jv''(z)]      im[jv''(z)]'
        write(*,*)' ----------------------------------',
     &            '-----------------------------------'
        do 10 k=0,nm,ns
10         write(*,20) k,cbj(k),cdj(k)
        write(*,*)
        write(*,*)'  n       re[yv(z)]       im[yv(z)]',
     &            '       re[yv''(z)]      im[yv''(z)]'
        write(*,*)' ----------------------------------',
     &            '-----------------------------------'
        do 15 k=0,nm,ns
15         write(*,20) k,cby(k),cdy(k)
20      format(1x,i3,2x,4d16.8)
25      format(8x,'v = n+v0',',  v0 =',f5.2,',  z =',f7.2,' +',f7.2,'i')
        end


        subroutine cjyva(v,z,vm,cbj,cdj,cby,cdy)
c
c       ===========================================================
c       purpose: compute bessel functions jv(z), yv(z) and their 
c                derivatives for a complex argument
c       input :  z --- complex argument
c                v --- order of jv(z) and yv(z)
c                      ( v = n+v0, n = 0,1,2,..., 0 ף v0 < 1 )
c       output:  cbj(n) --- jn+v0(z)
c                cdj(n) --- jn+v0'(z)
c                cby(n) --- yn+v0(z)
c                cdy(n) --- yn+v0'(z)
c                vm --- highest order computed
c       routines called:
c            (1) gamma for computing the gamma function
c            (2) msta1 and msta2 for computing the starting
c                point for backward recurrence
c       ===========================================================
c
        implicit double precision (a,b,g,o-y)
        implicit complex*16 (c,z)
        dimension cbj(0:*),cdj(0:*),cby(0:*),cdy(0:*)
        pi=3.141592653589793d0
        rp2=.63661977236758d0
        ci=(0.0d0,1.0d0)
        a0=cdabs(z)
        z1=z
        z2=z*z
        n=int(v)
        v0=v-n
        pv0=pi*v0
        pv1=pi*(1.0d0+v0)
        if (a0.lt.1.0d-100) then
           do 10 k=0,n
              cbj(k)=(0.0d0,0.0d0)
              cdj(k)=(0.0d0,0.0d0)
              cby(k)=-(1.0d+300,0.0d0)
10            cdy(k)=(1.0d+300,0.0d0)
           if (v0.eq.0.0) then
              cbj(0)=(1.0d0,0.0d0)
              cdj(1)=(0.5d0,0.0d0)
           else
              cdj(0)=(1.0d+300,0.0d0)
           endif
           vm=v                     
           return
        endif
        if (real(z).lt.0.0) z1=-z
        if (a0.le.12.0) then
           do 25 l=0,1
              vl=v0+l
              cjvl=(1.0d0,0.0d0)
              cr=(1.0d0,0.0d0)
              do 15 k=1,40
                 cr=-0.25d0*cr*z2/(k*(k+vl))
                 cjvl=cjvl+cr
                 if (cdabs(cr).lt.cdabs(cjvl)*1.0d-15) go to 20
15            continue
20            vg=1.0d0+vl
              call gamma(vg,ga)
              ca=(0.5d0*z1)**vl/ga
              if (l.eq.0) cjv0=cjvl*ca
              if (l.eq.1) cjv1=cjvl*ca
25         continue
        else
           k0=11
           if (a0.ge.35.0) k0=10
           if (a0.ge.50.0) k0=8
           do 40 j=0,1
              vv=4.0d0*(j+v0)*(j+v0)
              cpz=(1.0d0,0.0d0)
              crp=(1.0d0,0.0d0)
              do 30 k=1,k0
                 crp=-0.78125d-2*crp*(vv-(4.0*k-3.0)**2.0)*(vv-
     &               (4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*z2)
30               cpz=cpz+crp
              cqz=(1.0d0,0.0d0)
              crq=(1.0d0,0.0d0)
              do 35 k=1,k0
                 crq=-0.78125d-2*crq*(vv-(4.0*k-1.0)**2.0)*(vv-
     &               (4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*z2)
35               cqz=cqz+crq
              cqz=0.125d0*(vv-1.0)*cqz/z1
              zk=z1-(0.5d0*(j+v0)+0.25d0)*pi
              ca0=cdsqrt(rp2/z1)
              cck=cdcos(zk)
              csk=cdsin(zk)
              if (j.eq.0) then
                 cjv0=ca0*(cpz*cck-cqz*csk)
                 cyv0=ca0*(cpz*csk+cqz*cck)
              else if (j.eq.1) then
                 cjv1=ca0*(cpz*cck-cqz*csk)
                 cyv1=ca0*(cpz*csk+cqz*cck)
              endif
40         continue
        endif
        if (a0.le.12.0) then
           if (v0.ne.0.0) then
              do 55 l=0,1
                 vl=v0+l
                 cjvl=(1.0d0,0.0d0)
                 cr=(1.0d0,0.0d0)
                 do 45 k=1,40
                    cr=-0.25d0*cr*z2/(k*(k-vl))
                    cjvl=cjvl+cr
                    if (cdabs(cr).lt.cdabs(cjvl)*1.0d-15) go to 50
45               continue
50               vg=1.0d0-vl
                 call gamma(vg,gb)
                 cb=(2.0d0/z1)**vl/gb
                 if (l.eq.0) cju0=cjvl*cb
                 if (l.eq.1) cju1=cjvl*cb
55            continue
              cyv0=(cjv0*dcos(pv0)-cju0)/dsin(pv0)
              cyv1=(cjv1*dcos(pv1)-cju1)/dsin(pv1)
           else
              cec=cdlog(z1/2.0d0)+.5772156649015329d0
              cs0=(0.0d0,0.0d0)
              w0=0.0d0
              cr0=(1.0d0,0.0d0)
              do 60 k=1,30
                 w0=w0+1.0d0/k
                 cr0=-0.25d0*cr0/(k*k)*z2
60               cs0=cs0+cr0*w0
              cyv0=rp2*(cec*cjv0-cs0)
              cs1=(1.0d0,0.0d0)
              w1=0.0d0
              cr1=(1.0d0,0.0d0)
              do 65 k=1,30
                 w1=w1+1.0d0/k
                 cr1=-0.25d0*cr1/(k*(k+1))*z2
65               cs1=cs1+cr1*(2.0d0*w1+1.0d0/(k+1.0d0))
              cyv1=rp2*(cec*cjv1-1.0d0/z1-0.25d0*z1*cs1)
           endif
        endif
        if (real(z).lt.0.0d0) then
           cfac0=cdexp(pv0*ci)
           cfac1=cdexp(pv1*ci)
           if (dimag(z).lt.0.0d0) then
              cyv0=cfac0*cyv0-2.0d0*ci*dcos(pv0)*cjv0
              cyv1=cfac1*cyv1-2.0d0*ci*dcos(pv1)*cjv1
              cjv0=cjv0/cfac0
              cjv1=cjv1/cfac1
           else if (dimag(z).gt.0.0d0) then
              cyv0=cyv0/cfac0+2.0d0*ci*dcos(pv0)*cjv0
              cyv1=cyv1/cfac1+2.0d0*ci*dcos(pv1)*cjv1
              cjv0=cfac0*cjv0
              cjv1=cfac1*cjv1
           endif
        endif
        cbj(0)=cjv0
        cbj(1)=cjv1
        if (n.ge.2.and.n.le.int(0.25*a0)) then
           cf0=cjv0
           cf1=cjv1
           do 70 k=2,n
              cf=2.0d0*(k+v0-1.0d0)/z*cf1-cf0
              cbj(k)=cf
              cf0=cf1
70            cf1=cf
        else if (n.ge.2) then
           m=msta1(a0,200)
           if (m.lt.n) then
              n=m
           else
              m=msta2(a0,n,15)
           endif
           cf2=(0.0d0,0.0d0)
           cf1=(1.0d-100,0.0d0)
           do 75 k=m,0,-1
              cf=2.0d0*(v0+k+1.0d0)/z*cf1-cf2
              if (k.le.n) cbj(k)=cf
              cf2=cf1
75            cf1=cf
           if (cdabs(cjv0).gt.cdabs(cjv1)) cs=cjv0/cf
           if (cdabs(cjv0).le.cdabs(cjv1)) cs=cjv1/cf2
           do 80 k=0,n
80            cbj(k)=cs*cbj(k)
        endif
        cdj(0)=v0/z*cbj(0)-cbj(1)
        do 85 k=1,n
85         cdj(k)=-(k+v0)/z*cbj(k)+cbj(k-1)
        cby(0)=cyv0
        cby(1)=cyv1
        ya0=cdabs(cyv0)
        lb=0
        cg0=cyv0
        cg1=cyv1
        do 90 k=2,n
           cyk=2.0d0*(v0+k-1.0d0)/z*cg1-cg0
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
           ch0=2.0d0*(k+v0)/z*ch1-ch2
           ch2=ch1
100        ch1=ch0
        cp12=ch0
        cp22=ch2
        ch2=(0.0d0,0.0d0)
        ch1=(1.0d0,0.0d0)
        do 105 k=lb,1,-1
           ch0=2.0d0*(k+v0)/z*ch1-ch2
           ch2=ch1
105        ch1=ch0
        cp11=ch0
        cp21=ch2
        if (lb.eq.n) cbj(lb+1)=2.0d0*(lb+v0)/z*cbj(lb)-cbj(lb-1)
        if (cdabs(cbj(0)).gt.cdabs(cbj(1))) then
           cby(lb+1)=(cbj(lb+1)*cyv0-2.0d0*cp11/(pi*z))/cbj(0)
           cby(lb)=(cbj(lb)*cyv0+2.0d0*cp12/(pi*z))/cbj(0)
        else
           cby(lb+1)=(cbj(lb+1)*cyv1-2.0d0*cp21/(pi*z))/cbj(1)
           cby(lb)=(cbj(lb)*cyv1+2.0d0*cp22/(pi*z))/cbj(1)
        endif
        cyl2=cby(lb+1)
        cyl1=cby(lb)
        do 110 k=lb-1,0,-1
           cylk=2.0d0*(k+v0+1.0d0)/z*cyl1-cyl2
           cby(k)=cylk
           cyl2=cyl1
110        cyl1=cylk
        cyl1=cby(lb)
        cyl2=cby(lb+1)
        do 115 k=lb+1,n-1
           cylk=2.0d0*(k+v0)/z*cyl2-cyl1
           cby(k+1)=cylk
           cyl1=cyl2
115        cyl2=cylk
        do 120 k=2,n
           wa=cdabs(cby(k))
           if (wa.lt.cdabs(cby(k-1))) lb=k
120     continue
        go to 95
125     cdy(0)=v0/z*cby(0)-cby(1)
        do 130 k=1,n
130        cdy(k)=cby(k-1)-(k+v0)/z*cby(k)
        vm=n+v0
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
