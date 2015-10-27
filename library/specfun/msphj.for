        program msphj
c
c       =======================================================
c       purpose: this program computes the spherical bessel 
c                functions jn(x) and jn'(x) using subroutine
c                sphj
c       input :  x --- argument of jn(x)
c                n --- order of jn(x)  ( n = 0,1,תתת,ף 250 )
c       output:  sj(n) --- jn(x)
c                dj(n) --- jn'(x)
c       example:   x =10.0
c                  n          jn(x)               jn(x)
c                --------------------------------------------
c                  0    -.5440211109d-01    -.7846694180d-01
c                  1     .7846694180d-01    -.7009549945d-01
c                  2     .7794219363d-01     .5508428371d-01
c                  3    -.3949584498d-01     .9374053162d-01
c                  4    -.1055892851d+00     .1329879757d-01
c                  5    -.5553451162d-01    -.7226857814d-01
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension sj(0:250),dj(0:250)
        write(*,*)'please enter n and x '
        read(*,*)n,x
        write(*,30)n,x
        if (n.le.10) then
           ns=1
        else
           write(*,*)'please enter order step ns'
           read(*,*)ns
        endif
        call sphj(n,x,nm,sj,dj)
        write(*,*)
        write(*,*)'  n          jn(x)               jn''(x)'
        write(*,*)'--------------------------------------------'
        do 10 k=0,nm,ns
10         write(*,20)k,sj(k),dj(k)
20      format(1x,i3,2d20.10)
30      format(3x,6hnmax =,i3,',     ',2hx=,f5.1)
        end


        subroutine sphj(n,x,nm,sj,dj)
c
c       =======================================================
c       purpose: compute spherical bessel functions jn(x) and
c                their derivatives
c       input :  x --- argument of jn(x)
c                n --- order of jn(x)  ( n = 0,1,תתת )
c       output:  sj(n) --- jn(x)
c                dj(n) --- jn'(x)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the starting
c                point for backward recurrence
c       =======================================================
c
        implicit double precision (a-h,o-z)
        dimension sj(0:n),dj(0:n)
        nm=n
        if (dabs(x).eq.1.0d-100) then
           do 10 k=0,n
              sj(k)=0.0d0
10            dj(k)=0.0d0
           sj(0)=1.0d0
           dj(1)=.3333333333333333d0
           return
        endif
        sj(0)=dsin(x)/x
        sj(1)=(sj(0)-dcos(x))/x
        if (n.ge.2) then
           sa=sj(0)
           sb=sj(1)
           m=msta1(x,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(x,n,15)
           endif
           f0=0.0d0
           f1=1.0d0-100
           do 15 k=m,0,-1
              f=(2.0d0*k+3.0d0)*f1/x-f0
              if (k.le.nm) sj(k)=f
              f0=f1
15            f1=f
           if (dabs(sa).gt.dabs(sb)) cs=sa/f
           if (dabs(sa).le.dabs(sb)) cs=sb/f0
           do 20 k=0,nm
20            sj(k)=cs*sj(k)
        endif      
        dj(0)=(dcos(x)-dsin(x)/x)/x
        do 25 k=1,nm
25         dj(k)=sj(k-1)-(k+1.0d0)*sj(k)/x
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
