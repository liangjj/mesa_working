        program mch12n
c
c       =====================================================
c       purpose: this program computes hankel functions of
c                the first and second kinds and their 
c                derivatives for a complex argument using
c                subroutine ch12n
c       input :  z --- complex argument
c                n --- order of hn(1)(z) and hn(2)(z)
c                      ( n = 0,1,תתת, n ף 250 )
c       output:  chf1(n) --- hn(1)(z)
c                chd1(n) --- hn(1)'(z)
c                chf2(n) --- hn(2)(z)
c                chd2(n) --- hn(2)'(z)
c       =====================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension chf1(0:250),chd1(0:250),chf2(0:250),chd2(0:250)
        write(*,*)'  please enter n, x and y (z=x+iy) '
        read(*,*)n,x,y
        write(*,45)x,y,n
        z=cmplx(x,y)
        if (n.le.8) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call ch12n(n,z,nm,chf1,chd1,chf2,chd2)
        write(*,*)
        write(*,*)'   n     re[hn(1)(z)]     im[hn(1)(z)]',
     &            '      re[hn(1)''(z)]     im[hn(1)''(z)]'
        write(*,*)' -------------------------------------',
     &               '---------------------------------------'
        do 30 k=0,nm,ns
30         write(*,40)k,chf1(k),chd1(k)
        write(*,*)
        write(*,*)'   n     re[hn(2)(z)]     im[hn(2)(z)]',
     &            '      re[hn(2)''(z)]     im[hn(2)''(z)]'
        write(*,*)' -------------------------------------',
     &            '---------------------------------------'
        do 35 k=0,nm,ns
35         write(*,40)k,chf2(k),chd2(k)
40      format(1x,i4,4d18.10)
45      format(3x,3hz =,f8.3,' + i ',f8.3,' ,',6x,6hnmax =,i4)
        end


        subroutine ch12n(n,z,nm,chf1,chd1,chf2,chd2)
c
c       ====================================================
c       purpose: compute hankel functions of the first and
c                second kinds and their derivatives for a
c                complex argument
c       input :  z --- complex argument
c                n --- order of hn(1)(z) and hn(2)(z)
c       output:  chf1(n) --- hn(1)(z)
c                chd1(n) --- hn(1)'(z)
c                chf2(n) --- hn(2)(z)
c                chd2(n) --- hn(2)'(z)
c                nm --- highest order computed
c       routines called:
c             (1) cjynb for computing jn(z) and yn(z)
c             (2) ciknb for computing in(z) and kn(z)
c       ====================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension cbj(0:250),cdj(0:250),cby(0:250),cdy(0:250),
     &            cbi(0:250),cdi(0:250),cbk(0:250),cdk(0:250)
        dimension chf1(0:n),chd1(0:n),chf2(0:n),chd2(0:n)
        ci=(0.0d0,1.0d0)
        pi=3.141592653589793d0
        if (dimag(z).lt.0.0d0) then
           call cjynb(n,z,nm,cbj,cdj,cby,cdy)
           do 10 k=0,nm
              chf1(k)=cbj(k)+ci*cby(k)
10            chd1(k)=cdj(k)+ci*cdy(k)
           zi=ci*z
           call ciknb(n,zi,nm,cbi,cdi,cbk,cdk)
           cfac=-2.0d0/(pi*ci)
           do 15 k=0,nm
              chf2(k)=cfac*cbk(k)
              chd2(k)=cfac*ci*cdk(k)
15            cfac=cfac*ci
        else if (dimag(z).gt.0.0d0) then
           zi=-ci*z
           call ciknb(n,zi,nm,cbi,cdi,cbk,cdk)
           cf1=-ci
           cfac=2.0d0/(pi*ci)
           do 20 k=0,nm
              chf1(k)=cfac*cbk(k)
              chd1(k)=-cfac*ci*cdk(k)
20            cfac=cfac*cf1
           call cjynb(n,z,nm,cbj,cdj,cby,cdy)
           do 25 k=0,nm
              chf2(k)=cbj(k)-ci*cby(k)
25            chd2(k)=cdj(k)-ci*cdy(k)
        else
           call cjynb(n,z,nm,cbj,cdj,cby,cdy)
           do 30 k=0,nm
              chf1(k)=cbj(k)+ci*cby(k)
              chd1(k)=cdj(k)+ci*cdy(k)
              chf2(k)=cbj(k)-ci*cby(k)
30            chd2(k)=cdj(k)-ci*cdy(k)
        endif
        return
        end


        subroutine cjynb(n,z,nm,cbj,cdj,cby,cdy)
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
c       routines called:
c                msta1 and msta2 to calculate the starting
c                point for backward recurrence
c       =======================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension cbj(0:n),cdj(0:n),cby(0:n),cdy(0:n),
     &            a(4),b(4),a1(4),b1(4)
        el=0.5772156649015329d0
        pi=3.141592653589793d0
        r2p=.63661977236758d0
        y0=dabs(dimag(z))
        a0=cdabs(z)
        nm=n
        if (a0.lt.1.0d-100) then
           do 10 k=0,n
              cbj(k)=(0.0d0,0.0d0)
              cdj(k)=(0.0d0,0.0d0)
              cby(k)=-(1.0d+300,0.0d0)
10            cdy(k)=(1.0d+300,0.0d0)
           cbj(0)=(1.0d0,0.0d0)
           cdj(1)=(0.5d0,0.0d0)
           return
        endif
        if (a0.le.300.d0.or.n.gt.80) then
           if (n.eq.0) nm=1
           m=msta1(a0,200)
           if (m.lt.nm) then
              nm=m
           else
              m=msta2(a0,nm,15)
           endif
           cbs=(0.0d0,0.0d0)
           csu=(0.0d0,0.0d0)
           csv=(0.0d0,0.0d0)
           cf2=(0.0d0,0.0d0)
           cf1=(1.0d-100,0.0d0)
           do 15 k=m,0,-1
              cf=2.0d0*(k+1.0d0)/z*cf1-cf2
              if (k.le.nm) cbj(k)=cf
              if (k.eq.2*int(k/2).and.k.ne.0) then
                 if (y0.le.1.0d0) then
                    cbs=cbs+2.0d0*cf
                 else
                    cbs=cbs+(-1)**(k/2)*2.0d0*cf
                 endif
                 csu=csu+(-1)**(k/2)*cf/k
              else if (k.gt.1) then
                 csv=csv+(-1)**(k/2)*k/(k*k-1.0d0)*cf
              endif
              cf2=cf1
15            cf1=cf
           if (y0.le.1.0d0) then
              cs0=cbs+cf
           else
              cs0=(cbs+cf)/cdcos(z)
           endif
           do 20 k=0,nm
20            cbj(k)=cbj(k)/cs0
           ce=cdlog(z/2.0d0)+el
           cby(0)=r2p*(ce*cbj(0)-4.0d0*csu/cs0)
           cby(1)=r2p*(-cbj(0)/z+(ce-1.0d0)*cbj(1)-4.0d0*csv/cs0)
        else
           data a/-.7031250000000000d-01,.1121520996093750d+00,
     &            -.5725014209747314d+00,.6074042001273483d+01/
           data b/ .7324218750000000d-01,-.2271080017089844d+00,
     &             .1727727502584457d+01,-.2438052969955606d+02/
           data a1/.1171875000000000d+00,-.1441955566406250d+00,
     &             .6765925884246826d+00,-.6883914268109947d+01/
           data b1/-.1025390625000000d+00,.2775764465332031d+00,
     &             -.1993531733751297d+01,.2724882731126854d+02/
           ct1=z-0.25d0*pi
           cp0=(1.0d0,0.0d0)
           do 25 k=1,4
25            cp0=cp0+a(k)*z**(-2*k)
           cq0=-0.125d0/z
           do 30 k=1,4
30            cq0=cq0+b(k)*z**(-2*k-1)
           cu=cdsqrt(r2p/z)
           cbj0=cu*(cp0*cdcos(ct1)-cq0*cdsin(ct1))
           cby0=cu*(cp0*cdsin(ct1)+cq0*cdcos(ct1))
           cbj(0)=cbj0
           cby(0)=cby0
           ct2=z-0.75d0*pi
           cp1=(1.0d0,0.0d0)
           do 35 k=1,4
35            cp1=cp1+a1(k)*z**(-2*k)
           cq1=0.375d0/z
           do 40 k=1,4
40            cq1=cq1+b1(k)*z**(-2*k-1)
           cbj1=cu*(cp1*cdcos(ct2)-cq1*cdsin(ct2))
           cby1=cu*(cp1*cdsin(ct2)+cq1*cdcos(ct2))
           cbj(1)=cbj1
           cby(1)=cby1
           do 45 k=2,nm
              cbjk=2.0d0*(k-1.0d0)/z*cbj1-cbj0
              cbj(k)=cbjk
              cbj0=cbj1
45            cbj1=cbjk
        endif
        cdj(0)=-cbj(1)
        do 50 k=1,nm
50         cdj(k)=cbj(k-1)-k/z*cbj(k)
        if (cdabs(cbj(0)).gt.1.0d0) then
           cby(1)=(cbj(1)*cby(0)-2.0d0/(pi*z))/cbj(0)
        endif
        do 55 k=2,nm
           if (cdabs(cbj(k-1)).ge.cdabs(cbj(k-2))) then
              cyy=(cbj(k)*cby(k-1)-2.0d0/(pi*z))/cbj(k-1)
           else
              cyy=(cbj(k)*cby(k-2)-4.0d0*(k-1.0d0)/(pi*z*z))/cbj(k-2)
           endif
           cby(k)=cyy
55      continue
        cdy(0)=-cby(1)
        do 60 k=1,nm
60         cdy(k)=cby(k-1)-k/z*cby(k)
        return
        end


        subroutine ciknb(n,z,nm,cbi,cdi,cbk,cdk)
c
c       ============================================================
c       purpose: compute modified bessel functions in(z) and kn(z),
c                and their derivatives for a complex argument
c       input:   z --- complex argument
c                n --- order of in(z) and kn(z)
c       output:  cbi(n) --- in(z)
c                cdi(n) --- in'(z)
c                cbk(n) --- kn(z)
c                cdk(n) --- kn'(z)
c                nm --- highest order computed
c       routones called:
c                msta1 and msta2 to compute the starting point for
c                backward recurrence
c       ===========================================================
c
        implicit double precision (a,b,d-h,o-y)
        implicit complex*16 (c,z)
        dimension cbi(0:n),cdi(0:n),cbk(0:n),cdk(0:n)
        pi=3.141592653589793d0
        el=0.57721566490153d0
        a0=cdabs(z)
        nm=n
        if (a0.lt.1.0d-100) then
           do 10 k=0,n
              cbi(k)=(0.0d0,0.0d0)
              cbk(k)=(1.0d+300,0.0d0)
              cdi(k)=(0.0d0,0.0d0)
10            cdk(k)=-(1.0d+300,0.0d0)
           cbi(0)=(1.0d0,0.0d0)
           cdi(1)=(0.5d0,0.0d0)
           return
        endif
        z1=z
        ci=(0.0d0,1.0d0)
        if (real(z).lt.0.0) z1=-z
        if (n.eq.0) nm=1
        m=msta1(a0,200)
        if (m.lt.nm) then
           nm=m
        else
           m=msta2(a0,nm,15)
        endif
        cbs=0.0d0
        csk0=0.0d0
        cf0=0.0d0
        cf1=1.0d-100
        do 15 k=m,0,-1
           cf=2.0d0*(k+1.0d0)*cf1/z1+cf0
           if (k.le.nm) cbi(k)=cf
           if (k.ne.0.and.k.eq.2*int(k/2)) csk0=csk0+4.0d0*cf/k
           cbs=cbs+2.0d0*cf
           cf0=cf1
15         cf1=cf
        cs0=cdexp(z1)/(cbs-cf)
        do 20 k=0,nm
20         cbi(k)=cs0*cbi(k)
        if (a0.le.9.0) then
           cbk(0)=-(cdlog(0.5d0*z1)+el)*cbi(0)+cs0*csk0
           cbk(1)=(1.0d0/z1-cbi(1)*cbk(0))/cbi(0)
        else
           ca0=cdsqrt(pi/(2.0d0*z1))*cdexp(-z1)
           k0=16
           if (a0.ge.25.0) k0=10
           if (a0.ge.80.0) k0=8
           if (a0.ge.200.0) k0=6
           do 30 l=0,1
              cbkl=1.0d0
              vt=4.0d0*l
              cr=(1.0d0,0.0d0)
              do 25 k=1,k0
                 cr=0.125d0*cr*(vt-(2.0*k-1.0)**2)/(k*z1)
25               cbkl=cbkl+cr
              cbk(l)=ca0*cbkl
30         continue
        endif
        cg0=cbk(0)
        cg1=cbk(1)
        do 35 k=2,nm
           cg=2.0d0*(k-1.0d0)/z1*cg1+cg0
           cbk(k)=cg
           cg0=cg1
35         cg1=cg
        if (real(z).lt.0.0) then
           fac=1.0d0
           do 45 k=0,nm
              if (dimag(z).lt.0.0) then
                 cbk(k)=fac*cbk(k)+ci*pi*cbi(k)
              else
                 cbk(k)=fac*cbk(k)-ci*pi*cbi(k)
              endif
              cbi(k)=fac*cbi(k)
              fac=-fac
45         continue
        endif
        cdi(0)=cbi(1)
        cdk(0)=-cbk(1)
        do 50 k=1,nm
           cdi(k)=cbi(k-1)-k/z*cbi(k)
50         cdk(k)=-cbk(k-1)-k/z*cbk(k)
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
