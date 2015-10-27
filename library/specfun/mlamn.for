        program mlamn
c
c       ====================================================
c       purpose: this program computes the lambda functions 
c                and their derivatives using subroutine
c                lamn
c       input:   x --- argument of lambda function
c                n --- order of lambda function
c                      ( n = 0,1,..., n ó 250 )
c       output:  bl(n) --- lambda function of order n
c                dl(n) --- derivative of lambda function
c       example: nmax = 5,  x = 10.00
c
c                 n       lambda(x)        lambda'(x)
c                ---------------------------------------
c                 0    -.24593576d+00    -.43472746d-01
c                 1     .86945492d-02    -.50926063d-01
c                 2     .20370425d-01    -.46703503d-02
c                 3     .28022102d-02     .10540929d-01
c                 4    -.84327431d-02     .89879627d-02
c                 5    -.89879627d-02     .55521954d-03
c       ====================================================
        implicit double precision (a-h,o-z)
        dimension bl(0:250),dl(0:250)
        write(*,*)'  please enter n,x = ?'
        read(*,*)n,x
        write(*,15)n,x
        if (n.le.10) then
           ns=1
        else
           write(*,*)'  please enter order step ns'
           read(*,*)ns
        endif
        call lamn(n,x,nm,bl,dl)
        write(*,*)
        write(*,*) '  n       lambda(x)        lambda''(x)'
        write(*,*)' ---------------------------------------'
        do 10 k=0,nm,ns
10         write(*,20)k,bl(k),dl(k)
15      format(1x,3hn =,i4,6x,3hx =,f8.2)
20      format(1x,i3,2d18.8)
        end


        subroutine lamn(n,x,nm,bl,dl)
c
c       =========================================================
c       purpose: compute lambda functions and their derivatives
c       input:   x --- argument of lambda function
c                n --- order of lambda function
c       output:  bl(n) --- lambda function of order n
c                dl(n) --- derivative of lambda function
c                nm --- highest order computed
c       routines called:
c                msta1 and msta2 for computing the start
c                point for backward recurrence
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension bl(0:n),dl(0:n)
        nm=n
        if (dabs(x).lt.1.0d-100) then
           do 10 k=0,n
              bl(k)=0.0d0
10            dl(k)=0.0d0
           bl(0)=1.0d0
           dl(1)=0.5d0
           return
        endif
        if (x.le.12.0d0) then
           x2=x*x
           do 25 k=0,n
              bk=1.0d0
              r=1.0d0
              do 15 i=1,50
                 r=-0.25d0*r*x2/(i*(i+k))
                 bk=bk+r
                 if (dabs(r).lt.dabs(bk)*1.0d-15) go to 20
15            continue
20            bl(k)=bk
25            if (k.ge.1) dl(k-1)=-0.5d0*x/k*bk
           uk=1.0d0
           r=1.0d0
           do 30 i=1,50
              r=-0.25d0*r*x2/(i*(i+n+1.0d0))
              uk=uk+r
              if (dabs(r).lt.dabs(uk)*1.0d-15) go to 35
30            continue
35         dl(n)=-0.5d0*x/(n+1.0d0)*uk
           return
        endif
        if (n.eq.0) nm=1          
        m=msta1(x,200)
        if (m.lt.nm) then
           nm=m
        else
           m=msta2(x,nm,15)
        endif
        bs=0.0d0
        f0=0.0d0
        f1=1.0d-100
        do 40 k=m,0,-1
           f=2.0d0*(k+1.0d0)*f1/x-f0
           if (k.le.nm) bl(k)=f
           if (k.eq.2*int(k/2)) bs=bs+2.0d0*f
           f0=f1
40         f1=f
        bg=bs-f
        do 45 k=0,nm
45         bl(k)=bl(k)/bg
        r0=1.0d0
        do 50 k=1,nm
           r0=2.0d0*r0*k/x
50         bl(k)=r0*bl(k)
        dl(0)=-0.5d0*x*bl(1)
        do 55 k=1,nm
55         dl(k)=2.0d0*k/x*(bl(k-1)-bl(k))
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
