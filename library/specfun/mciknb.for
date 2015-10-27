        program mciknb
c
c       =============================================================
c       purpose: this program computes the modified bessel functions 
c                in(z) and kn(z), and their derivatives for a
c                complex argument using subroutine ciknb
c       input:   z --- complex argument
c                n --- order of in(z) and kn(z)
c                      ( n = 0,1,תתת, n ף 250 )
c       output:  cbi(n) --- in(z)
c                cdi(n) --- in'(z)
c                cbk(n) --- kn(z)
c                cdk(n) --- kn'(z)
c       example: nmax = 5,   z = 4.0 + i 2.0
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
        dimension cbi(0:250),cdi(0:250),cbk(0:250),cdk(0:250)
        write(*,*)'  please enter n, x and y (z = x+iy) '
        read(*,*)n,x,y
        write(*,25)n,x,y
        z=cmplx(x,y)
        write(*,*)
        if (n.le.8) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call ciknb(n,z,nm,cbi,cdi,cbk,cdk)
        write(*,*)'   n      re[in(z)]       im[in(z)]',
     &            '      re[in''(z)]       im[in''(z)]'
        write(*,*)' ---------------------------------',
     &            '------------------------------------'
        do 10 k=0,nm,ns
           write(*,20)k,cbi(k),cdi(k)
10      continue
        write(*,*)
        write(*,*)'   n      re[kn(z)]       im[kn(z)]',
     &            '      re[kn''(z)]       im[kn''(z)]'
        write(*,*)' ---------------------------------',
     &            '------------------------------------'
        do 15 k=0,nm,ns
           write(*,20)k,cbk(k),cdk(k)
15      continue
20      format(1x,1x,i3,1x,4d16.8)
25      format(3x,6hnmaz =,i3,',    ','z =', f6.1,' + i',f6.1)
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
