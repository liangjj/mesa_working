        program msphi
c
c       ======================================================
c       purpose: this program computes the modified spherical 
c                bessel functions of the first kind in(x) and 
c                in'(x) using subroutine sphi
c       input :  x --- argument of in(x)
c                n --- order of in(x) ( 0 ó n ó 250 )
c       output:  si(n) --- in(x)
c                di(n) --- in'(x)
c       example: x = 10.0
c                  n          in(x)               in'(x)
c                --------------------------------------------
c                  0     .1101323287d+04     .9911909633d+03
c                  1     .9911909633d+03     .9030850948d+03
c                  2     .8039659985d+03     .7500011637d+03
c                  3     .5892079640d+03     .5682828129d+03
c                  4     .3915204237d+03     .3934477522d+03
c                  5     .2368395827d+03     .2494166741d+03
c       ======================================================
c
        implicit double precision (a-h,o-z)
        dimension si(0:250),di(0:250)
        write(*,*)'please enter n and x '
        read(*,*)n,x
        write(*,30)n,x
        if (n.le.10) then
           ns=1
        else
           write(*,*) 'please enter order step ns'
           read(*,*) ns
        endif
        call sphi(n,x,nm,si,di)
        write(*,*)
        write(*,*)'  n          in(x)               in''(x)'
        write(*,*)'--------------------------------------------'
        do 10 k=0,nm,ns
10         write(*,20)k,si(k),di(k)
20      format(1x,i3,2d20.10)
30      format(3x,'nmax =',i3,',     ','x =',f6.1)
        end


        subroutine sphi(n,x,nm,si,di)
c
c       ========================================================
c       purpose: compute modified spherical bessel functions
c                of the first kind, in(x) and in'(x)
c       input :  x --- argument of in(x)
c                n --- order of in(x) ( n = 0,1,2,... )
c       output:  si(n) --- in(x)
c                di(n) --- in'(x)
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the starting
c                point for backward recurrence
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension si(0:n),di(0:n)
        nm=n
        if (dabs(x).lt.1.0d-100) then
           do 10 k=0,n
              si(k)=0.0d0
10            di(k)=0.0d0
           si(0)=1.0d0
           di(1)=0.333333333333333d0
           return
        endif
        si(0)=dsinh(x)/x
        si(1)=-(dsinh(x)/x-dcosh(x))/x
        si0=si(0)
        if (n.ge.2) then
           m=msta1(x,200)
           if (m.lt.n) then
              nm=m
           else
              m=msta2(x,n,15)
           endif
           f0=0.0d0
           f1=1.0d0-100
           do 15 k=m,0,-1
              f=(2.0d0*k+3.0d0)*f1/x+f0
              if (k.le.nm) si(k)=f
              f0=f1
15            f1=f
           cs=si0/f
           do 20 k=0,nm
20            si(k)=cs*si(k)
        endif
        di(0)=si(1)
        do 25 k=1,nm
25         di(k)=si(k-1)-(k+1.0d0)/x*si(k)
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
