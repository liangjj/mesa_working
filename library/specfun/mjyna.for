*deck bessel
        program bessel
c
c       ====================================================
c       purpose: this program computes bessel functions  
c                jn(x) and yn(x), and their derivatives 
c                using subroutine jyna
c       input :  x --- argument of jn(x) & yn(x)  ( x ע 0 )
c                n --- order of jn(x) & yn(x)
c                      ( n = 0,1,2,תתת, n ף 250 )
c       output:  bj(n) --- jn(x)
c                dj(n) --- jn'(x)
c                by(n) --- yn(x)
c                dy(n) --- yn'(x)
c       example:
c                x = 10.0
c
c                n        jn(x)           jn'(x)
c              -------------------------------------
c                0    -.2459358d+00   -.4347275d-01
c               10     .2074861d+00    .8436958d-01
c               20     .1151337d-04    .2011954d-04
c               30     .1551096d-11    .4396479d-11
c
c                n        yn(x)           yn'(x)
c              -------------------------------------
c                0     .5567117d-01   -.2490154d+00
c               10    -.3598142d+00    .1605149d+00
c               20    -.1597484d+04    .2737803d+04
c               30    -.7256142d+10    .2047617d+11
c       ====================================================
c
        implicit double precision (a-h,o-z)
        dimension bj(0:250),by(0:250),dj(0:250),dy(0:250)
        write(*,*)'  please enter n, x'
        read(*,*)n,x
        write(*,30)n,x
        if (n.le.8) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        write(*,*)
        call jyna(n,x,nm,bj,dj,by,dy)
        write(*,*)'  n        jn(x)           jn''(x)'
        write(*,*)'--------------------------------------'
        do 10 k=0,nm,ns
10         write(*,40)k,bj(k),dj(k)
        write(*,*)
        write(*,*)'  n        yn(x)           yn''(x)'
        write(*,*)'--------------------------------------'
        do 20 k=0,nm,ns
20         write(*,40)k,by(k),dy(k)
30      format(3x,6hnmax =,i3,',    ',3hx =,f6.1)
40      format(1x,i3,1x,2d16.7)
        end


        subroutine jyna(n,x,nm,bj,dj,by,dy)
c
c       ==========================================================
c       purpose: compute bessel functions jn(x) & yn(x) and
c                their derivatives
c       input :  x --- argument of jn(x) & yn(x)  ( x ע 0 )
c                n --- order of jn(x) & yn(x)
c       output:  bj(n) --- jn(x)
c                dj(n) --- jn'(x)
c                by(n) --- yn(x)
c                dy(n) --- yn'(x)
c                nm --- highest order computed
c       routines called:
c            (1) jy01b to calculate j0(x), j1(x), y0(x) & y1(x)
c            (2) msta1 and msta2 to calculate the starting 
c                point for backward recurrence
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension bj(0:n),by(0:n),dj(0:n),dy(0:n)
        nm=n
        if (x.lt.1.0d-100) then
           do 10 k=0,n
              bj(k)=0.0d0
              dj(k)=0.0d0
              by(k)=-1.0d+300
10            dy(k)=1.0d+300
           bj(0)=1.0d0
           dj(1)=0.5d0
           return
        endif
        call jy01b(x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1)
        bj(0)=bj0
        bj(1)=bj1
        by(0)=by0
        by(1)=by1
        dj(0)=dj0
        dj(1)=dj1
        dy(0)=dy0
        dy(1)=dy1
        if (n.le.1) return
        if (n.lt.int(0.9*x)) then
           do 20 k=2,n
              bjk=2.0d0*(k-1.0d0)/x*bj1-bj0
              bj(k)=bjk
              bj0=bj1
20            bj1=bjk
        else
           m=msta1(x,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(x,n,15)
           endif
           f2=0.0d0
           f1=1.0d-100
           do 30 k=m,0,-1
              f=2.0d0*(k+1.0d0)/x*f1-f2
              if (k.le.nm) bj(k)=f
              f2=f1
30            f1=f
           if (dabs(bj0).gt.dabs(bj1)) then
              cs=bj0/f
           else
              cs=bj1/f2
           endif
           do 40 k=0,nm
40            bj(k)=cs*bj(k)
        endif
        do 50 k=2,nm
50         dj(k)=bj(k-1)-k/x*bj(k)
        f0=by(0)
        f1=by(1)
        do 60 k=2,nm
           f=2.0d0*(k-1.0d0)/x*f1-f0
           by(k)=f
           f0=f1
60         f1=f
        do 70 k=2,nm
70         dy(k)=by(k-1)-k*by(k)/x
        return
        end


        subroutine jy01b(x,bj0,dj0,bj1,dj1,by0,dy0,by1,dy1)
c
c       =======================================================
c       purpose: compute bessel functions j0(x), j1(x), y0(x),
c                y1(x), and their derivatives
c       input :  x   --- argument of jn(x) & yn(x) ( x ע 0 )
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
        pi=3.141592653589793d0
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
        else if (x.le.4.0d0) then
           t=x/4.0d0
           t2=t*t
           bj0=((((((-.5014415d-3*t2+.76771853d-2)*t2
     &         -.0709253492d0)*t2+.4443584263d0)*t2
     &         -1.7777560599d0)*t2+3.9999973021d0)
     &         *t2-3.9999998721d0)*t2+1.0d0
           bj1=t*(((((((-.1289769d-3*t2+.22069155d-2)
     &         *t2-.0236616773d0)*t2+.1777582922d0)*t2
     &         -.8888839649d0)*t2+2.6666660544d0)*t2
     &         -3.9999999710d0)*t2+1.9999999998d0)
           by0=(((((((-.567433d-4*t2+.859977d-3)*t2
     &         -.94855882d-2)*t2+.0772975809d0)*t2
     &         -.4261737419d0)*t2+1.4216421221d0)*t2
     &         -2.3498519931d0)*t2+1.0766115157)*t2
     &         +.3674669052d0
           by0=2.0d0/pi*dlog(x/2.0d0)*bj0+by0
           by1=((((((((.6535773d-3*t2-.0108175626d0)*t2
     &         +.107657606d0)*t2-.7268945577d0)*t2
     &         +3.1261399273d0)*t2-7.3980241381d0)*t2
     &         +6.8529236342d0)*t2+.3932562018d0)*t2
     &         -.6366197726d0)/x
           by1=2.0d0/pi*dlog(x/2.0d0)*bj1+by1
        else
           t=4.0d0/x
           t2=t*t
           a0=dsqrt(2.0d0/(pi*x))
           p0=((((-.9285d-5*t2+.43506d-4)*t2-.122226d-3)*t2
     &        +.434725d-3)*t2-.4394275d-2)*t2+.999999997d0
           q0=t*(((((.8099d-5*t2-.35614d-4)*t2+.85844d-4)*t2
     &        -.218024d-3)*t2+.1144106d-2)*t2-.031249995d0)
           ta0=x-.25d0*pi
           bj0=a0*(p0*dcos(ta0)-q0*dsin(ta0))
           by0=a0*(p0*dsin(ta0)+q0*dcos(ta0))
           p1=((((.10632d-4*t2-.50363d-4)*t2+.145575d-3)*t2
     &        -.559487d-3)*t2+.7323931d-2)*t2+1.000000004d0
           q1=t*(((((-.9173d-5*t2+.40658d-4)*t2-.99941d-4)*t2
     &        +.266891d-3)*t2-.1601836d-2)*t2+.093749994d0)
           ta1=x-.75d0*pi
           bj1=a0*(p1*dcos(ta1)-q1*dsin(ta1))
           by1=a0*(p1*dsin(ta1)+q1*dcos(ta1))
        endif
        dj0=-bj1
        dj1=bj0-bj1/x
        dy0=-by1
        dy1=by0-by1/x
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
