        program mcjyvb
c
c       =============================================================
c       purpose: this program computes bessel functions jv(z), yv(z),
c                and their derivatives for a complex argument using
c                subroutine cjyvb
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
c    -------------------------------------------------------------------
c     0  -.13829878d+01  -.30855145d+00  -.18503756d+00   .13103689d+01
c     1   .82553327d-01  -.12848394d+01  -.12336901d+01   .45079506d-01
c     2   .10843924d+01  -.39871046d+00  -.33046401d+00  -.84574964d+00
c     3   .74348135d+00   .40665987d+00   .45318486d+00  -.42198992d+00
c     4   .17802266d+00   .44526939d+00   .39624497d+00   .97902890d-01
c     5  -.49008598d-01   .21085409d+00   .11784299d+00   .19422044d+00
c
c     n     re[yv(z)]      im[yv(z)]       re[yv'(z)]      im[yv'(z)]
c    -------------------------------------------------------------------
c     0   .34099851d+00  -.13440666d+01  -.13544477d+01  -.15470699d+00
c     1   .13323787d+01   .53735934d-01  -.21467271d-01  -.11807457d+01
c     2   .38393305d+00   .10174248d+01   .91581083d+00  -.33147794d+00
c     3  -.49924295d+00   .71669181d+00   .47786442d+00   .37321597d+00
c     4  -.57179578d+00   .27099289d+00  -.12111686d+00   .23405313d+00
c     5  -.25700924d+00   .24858555d+00  -.43023156d+00  -.13123662d+00
c       =============================================================
c
        implicit double precision (v,x,y)
        implicit complex*16 (c,z)
        dimension cbj(0:250),cdj(0:250),cby(0:250),cdy(0:250)
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
        call cjyvb(v,z,vm,cbj,cdj,cby,cdy)
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


        subroutine cjyvb(v,z,vm,cbj,cdj,cby,cdy)
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
        if (real(z).lt.0.0d0) z1=-z
        if (a0.le.12.0) then
           cjv0=(1.0d0,0.0d0)
           cr=(1.0d0,0.0d0)
           do 15 k=1,40
              cr=-0.25d0*cr*z2/(k*(k+v0))
              cjv0=cjv0+cr
              if (cdabs(cr).lt.cdabs(cjv0)*1.0d-15) go to 20
15         continue
20         vg=1.0d0+v0
           call gamma(vg,ga)
           ca=(0.5d0*z1)**v0/ga
           cjv0=cjv0*ca
        else
           k0=11
           if (a0.ge.35.0) k0=10
           if (a0.ge.50.0) k0=8
           vv=4.0d0*v0*v0
           cpz=(1.0d0,0.0d0)
           crp=(1.0d0,0.0d0)
           do 25 k=1,k0
              crp=-0.78125d-2*crp*(vv-(4.0*k-3.0)**2.0)*(vv-
     &            (4.0*k-1.0)**2.0)/(k*(2.0*k-1.0)*z2)
25            cpz=cpz+crp
           cqz=(1.0d0,0.0d0)
           crq=(1.0d0,0.0d0)
           do 30 k=1,k0
              crq=-0.78125d-2*crq*(vv-(4.0*k-1.0)**2.0)*(vv-
     &            (4.0*k+1.0)**2.0)/(k*(2.0*k+1.0)*z2)
30            cqz=cqz+crq
           cqz=0.125d0*(vv-1.0)*cqz/z1
           zk=z1-(0.5d0*v0+0.25d0)*pi
           ca0=cdsqrt(rp2/z1)
           cck=cdcos(zk)
           csk=cdsin(zk)
           cjv0=ca0*(cpz*cck-cqz*csk)
           cyv0=ca0*(cpz*csk+cqz*cck)
        endif
        if (a0.le.12.0) then
           if (v0.ne.0.0) then
              cjvn=(1.0d0,0.0d0)
              cr=(1.0d0,0.0d0)
              do 35 k=1,40
                 cr=-0.25d0*cr*z2/(k*(k-v0))
                 cjvn=cjvn+cr
                 if (cdabs(cr).lt.cdabs(cjvn)*1.0d-15) go to 40
35            continue
40            vg=1.0d0-v0
              call gamma(vg,gb)
              cb=(2.0d0/z1)**v0/gb
              cju0=cjvn*cb
              cyv0=(cjv0*dcos(pv0)-cju0)/dsin(pv0)
           else
              cec=cdlog(z1/2.0d0)+.5772156649015329d0
              cs0=(0.0d0,0.0d0)
              w0=0.0d0
              cr0=(1.0d0,0.0d0)
              do 45 k=1,30
                 w0=w0+1.0d0/k
                 cr0=-0.25d0*cr0/(k*k)*z2
45               cs0=cs0+cr0*w0
              cyv0=rp2*(cec*cjv0-cs0)
           endif
        endif
        if (n.eq.0) n=1
        m=msta1(a0,200)
        if (m.lt.n) then
           n=m
        else
           m=msta2(a0,n,15)
        endif
        cf2=(0.0d0,0.0d0)
        cf1=(1.0d-100,0.0d0)
        do 50 k=m,0,-1
           cf=2.0d0*(v0+k+1.0d0)/z1*cf1-cf2
           if (k.le.n) cbj(k)=cf
           cf2=cf1
50         cf1=cf
        cs=cjv0/cf
        do 55 k=0,n
55         cbj(k)=cs*cbj(k)
        if (real(z).lt.0.0d0) then
           cfac0=cdexp(pv0*ci)
           if (dimag(z).lt.0.0d0) then
              cyv0=cfac0*cyv0-2.0d0*ci*dcos(pv0)*cjv0
           else if (dimag(z).gt.0.0d0) then
              cyv0=cyv0/cfac0+2.0d0*ci*dcos(pv0)*cjv0
           endif
           do 60 k=0,n
              if (dimag(z).lt.0.0d0) then
                 cbj(k)=cdexp(-pi*(k+v0)*ci)*cbj(k)
              else if (dimag(z).gt.0.0d0) then
                 cbj(k)=cdexp(pi*(k+v0)*ci)*cbj(k)
              endif
60         continue
           z1=z1
        endif
        cby(0)=cyv0
        do 65 k=1,n
           cyy=(cbj(k)*cby(k-1)-2.0d0/(pi*z))/cbj(k-1)
           cby(k)=cyy
65      continue
        cdj(0)=v0/z*cbj(0)-cbj(1)
        do 70 k=1,n
70         cdj(k)=-(k+v0)/z*cbj(k)+cbj(k-1)
        cdy(0)=v0/z*cby(0)-cby(1)
        do 75 k=1,n
75         cdy(k)=cby(k-1)-(k+v0)/z*cby(k)
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
