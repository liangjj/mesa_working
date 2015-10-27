        program mcsphik
c
c       =============================================================
c       purpose: this program computes the modified spherical bessel 
c                functions and their derivatives for a complex
c                argument using subroutine csphik
c       input :  z --- complex argument
c                n --- order of in(z) & kn(z) ( 0 ó n ó 250 )
c       output:  csi(n) --- in(z)
c                cdi(n) --- in'(z)
c                csk(n) --- kn(z)
c                cdk(n) --- kn'(z)
c       example: z =4.0+i 2.0
c
c     n     re[in(z)]      im[in(z)]     re[in'(z)]     im[in'(z)]
c    ---------------------------------------------------------------
c     0   .2118080d+00   .6101922d+01  -.4439356d+00   .4900150d+01
c     1  -.4439356d+00   .4900150d+01  -.5906477d+00   .4053075d+01
c     2  -.9918756d+00   .3028652d+01  -.7574058d+00   .2785396d+01
c     3  -.9663859d+00   .1375561d+01  -.7689911d+00   .1541649d+01
c     4  -.6018277d+00   .4263967d+00  -.5777565d+00   .6482500d+00
c     5  -.2668530d+00   .6640148d-01  -.3214450d+00   .1866032d+00
c
c     n     re[kn(z)]      im[kn(z)]     re[kn'(z)]     im[kn'(z)]
c    ---------------------------------------------------------------
c     0  -.5010582d-02  -.4034862d-02   .6416184d-02   .4340777d-02
c     1  -.6416184d-02  -.4340777d-02   .8445211d-02   .4487936d-02
c     2  -.1016253d-01  -.4714473d-02   .1392804d-01   .4120703d-02
c     3  -.1893595d-01  -.3973987d-02   .2690088d-01   .3192843d-03
c     4  -.3945464d-01   .2977107d-02   .5690203d-01  -.1873044d-01
c     5  -.8727490d-01   .3689398d-01   .1220481d+00  -.9961483d-01
c       =============================================================
c
        implicit complex*16 (c,z)
        double precision x,y
        dimension csi(0:250),cdi(0:250),csk(0:250),cdk(0:250)
        write(*,*)'please enter n,x,y (z=x+iy) '
        read(*,*)n,x,y
        write(*,30)n,x,y
        z=cmplx(x,y)
        if (n.le.8) then
           ns=1
        else
           write(*,*)'please enter order step ns '
           read(*,*)ns
        endif
        call csphik(n,z,nm,csi,cdi,csk,cdk)
        write(*,*)
        write(*,*)'  n      re[in(z)]        im[in(z)]',
     &  '        re[in''(z)]       im[in''(z)]'
        write(*,*)'--------------------------------------------',
     &  '----------------------------'
        do 10 k=0,nm,ns
10         write(*,20)k,csi(k),cdi(k)
        write(*,*)
        write(*,*)'  n      re[kn(z)]        im[kn(z)]',
     &  '        re[kn''(z)]       im[kn''(z)]'
        write(*,*)'--------------------------------------------',
     &  '----------------------------'
        do 15 k=0,nm,ns
15         write(*,20)k,csk(k),cdk(k)
20      format(1x,i3,4d17.8)
30      format(3x,'nmaz =',i3,',     ','z = ',f8.1,'+ i',f8.1)
        end


        subroutine csphik(n,z,nm,csi,cdi,csk,cdk)
c
c       =======================================================
c       purpose: compute modified spherical bessel functions
c                and their derivatives for a complex argument
c       input :  z --- complex argument
c                n --- order of in(z) & kn(z) ( n = 0,1,2,... )
c       output:  csi(n) --- in(z)
c                cdi(n) --- in'(z)
c                csk(n) --- kn(z)
c                cdk(n) --- kn'(z)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the starting
c                point for backward recurrence
c       =======================================================
c
        implicit complex*16 (c,z)
        double precision a0,pi
        dimension csi(0:n),cdi(0:n),csk(0:n),cdk(0:n)
        pi=3.141592653589793d0
        a0=cdabs(z)            
        nm=n
        if (a0.lt.1.0d-60) then
           do 10 k=0,n
              csi(k)=0.0d0
              cdi(k)=0.0d0
              csk(k)=1.0d+300
10            cdk(k)=-1.0d+300
           csi(0)=1.0d0
           cdi(1)=0.3333333333333333d0
           return
        endif
        ci=cmplx(0.0d0,1.0d0)
        csinh=cdsin(ci*z)/ci
        ccosh=cdcos(ci*z)
        csi0=csinh/z
        csi1=(-csinh/z+ccosh)/z
        csi(0)=csi0
        csi(1)=csi1
        if (n.ge.2) then
           m=msta1(a0,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(a0,n,15)
           endif
           cf0=0.0d0
           cf1=1.0d0-100
           do 15 k=m,0,-1
              cf=(2.0d0*k+3.0d0)*cf1/z+cf0
              if (k.le.nm) csi(k)=cf
              cf0=cf1
15            cf1=cf
           if (cdabs(csi0).gt.cdabs(csi1)) cs=csi0/cf
           if (cdabs(csi0).le.cdabs(csi1)) cs=csi1/cf0
           do 20 k=0,nm
20            csi(k)=cs*csi(k)
        endif
        cdi(0)=csi(1)
        do 25 k=1,nm
25         cdi(k)=csi(k-1)-(k+1.0d0)*csi(k)/z
        csk(0)=0.5d0*pi/z*cdexp(-z)
        csk(1)=csk(0)*(1.0d0+1.0d0/z)
        do 30 k=2,nm
           if (cdabs(csi(k-1)).gt.cdabs(csi(k-2))) then
              csk(k)=(0.5d0*pi/(z*z)-csi(k)*csk(k-1))/csi(k-1)
           else
              csk(k)=(csi(k)*csk(k-2)+(k-0.5d0)*pi/z**3)/csi(k-2)
           endif
30      continue
        cdk(0)=-csk(1)
        do 35 k=1,nm
35         cdk(k)=-csk(k-1)-(k+1.0d0)*csk(k)/z
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
