        program mrctj
c
c       =======================================================
c       purpose: this program computes the riccati-bessel 
c                functions of the first kind, and their
c                derivatives using subroutine rctj
c       input:   x --- argument of riccati-bessel function
c                n --- order of jn(x)  ( 0 ó n ó 250 )
c       output:  rj(n) --- xújn(x)
c                dj(n) --- [xújn(x)]'
c       example: x = 10.0
c                  n        xújn(x)             [xújn(x)]'
c                --------------------------------------------
c                  0    -.5440211109d+00    -.8390715291d+00
c                  1     .7846694180d+00    -.6224880527d+00
c                  2     .7794219363d+00     .6287850307d+00
c                  3    -.3949584498d+00     .8979094712d+00
c                  4    -.1055892851d+01     .2739869063d-01
c                  5    -.5553451162d+00    -.7782202931d+00
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension rj(0:250),dj(0:250)
        write(*,*)'  please enter n and x '
        read(*,*)n,x
        write(*,30)n,x
        if (n.le.10) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        write(*,*)
        call rctj(n,x,nm,rj,dj)
        write(*,*)
        write(*,*)'  n        xújn(x)             [xújn(x)]'''
        write(*,*)'--------------------------------------------'
        do 10 k=0,nm,ns
10         write(*,20)k,rj(k),dj(k)
20      format(1x,i3,2d20.10)
30      format(3x,6hnmax =,i3,',    ',3hx =,f7.2)
        end


        subroutine rctj(n,x,nm,rj,dj)
c
c       ========================================================
c       purpose: compute riccati-bessel functions of the first
c                kind and their derivatives
c       input:   x --- argument of riccati-bessel function
c                n --- order of jn(x)  ( n = 0,1,2,... )
c       output:  rj(n) --- xújn(x)
c                dj(n) --- [xújn(x)]'
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the starting
c                point for backward recurrence
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension rj(0:n),dj(0:n)
        nm=n
        if (dabs(x).lt.1.0d-100) then
           do 10 k=0,n
              rj(k)=0.0d0
10            dj(k)=0.0d0
           dj(0)=1.0d0
           return
        endif
        rj(0)=dsin(x)
        rj(1)=rj(0)/x-dcos(x)
        rj0=rj(0)
        rj1=rj(1)
        if (n.ge.2) then
           m=msta1(x,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(x,n,15)
           endif
           f0=0.0d0
           f1=1.0d-100
           do 15 k=m,0,-1
              f=(2.0d0*k+3.0d0)*f1/x-f0
              if (k.le.nm) rj(k)=f
              f0=f1
15            f1=f
           if (dabs(rj0).gt.dabs(rj1)) cs=rj0/f
           if (dabs(rj0).le.dabs(rj1)) cs=rj1/f0
           do 20 k=0,nm
20            rj(k)=cs*rj(k)
        endif
        dj(0)=dcos(x)
        do 25 k=1,nm
25         dj(k)=-k*rj(k)/x+rj(k-1)
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
