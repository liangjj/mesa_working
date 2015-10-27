        program mlqnb
c
c       ===============================================================
c       purpose: this program computes the legendre functions qn(x) 
c                and qn'(x) using subroutine lqnb
c       input :  x  --- argument of qn(x)
c                n  --- degree of qn(x)  ( n = 0,1,תתת)
c       output:  qn(n) --- qn(x)
c                qd(n) --- qn'(x)
c       examples:     x1 = 0.50,    x2 = 2.50
c
c       n      qn(x1)        qn'(x1)       qn(x2)          qn'(x2)
c     ----------------------------------------------------------------
c       0     .54930614    1.33333333   .42364893d+00  -.19047619d+00
c       1    -.72534693    1.21597281   .59122325d-01  -.52541546d-01
c       2    -.81866327    -.84270745   .98842555d-02  -.13109214d-01
c       3    -.19865477   -2.87734353   .17695141d-02  -.31202687d-02
c       4     .44017453   -2.23329085   .32843271d-03  -.72261513d-03
c       5     .55508089    1.08422720   .62335892d-04  -.16437427d-03
c       ===============================================================
c
        implicit double precision (a-h,o-z)
        dimension qn(0:100),qd(0:100)
        write(*,*)'please enter nmax and x '
        read(*,*)n,x
        write(*,40)x
        write(*,*)
        write(*,*)'  n          qn(x)           qn''(x)'
        write(*,*)'--------------------------------------'
        call lqnb(n,x,qn,qd)
        do 10 k=0,n
           if (x.le.1.0) then
              write(*,20)k,qn(k),qd(k)
           else
              write(*,30)k,qn(k),qd(k)
           endif
10      continue
20      format(1x,i3,2f17.8)
30      format(1x,i3,2d17.8)
40      format(3x,'x =',f5.2)
        end


        subroutine lqnb(n,x,qn,qd)
c
c       ====================================================
c       purpose: compute legendre functions qn(x) & qn'(x)
c       input :  x  --- argument of qn(x)
c                n  --- degree of qn(x)  ( n = 0,1,2,תתת)
c       output:  qn(n) --- qn(x)
c                qd(n) --- qn'(x)
c       ====================================================
c
        implicit double precision (a-h,o-z)
        dimension qn(0:n),qd(0:n)
        eps=1.0d-14
        if (dabs(x).eq.1.0d0) then
           do 10 k=0,n
              qn(k)=1.0d+300
10            qd(k)=1.0d+300
           return
        endif
        if (x.le.1.021d0) then
           x2=dabs((1.0d0+x)/(1.0d0-x))
           q0=0.5d0*dlog(x2)
           q1=x*q0-1.0d0
           qn(0)=q0
           qn(1)=q1
           qd(0)=1.0d0/(1.0d0-x*x)
           qd(1)=qn(0)+x*qd(0)
           do 15 k=2,n
              qf=((2.0d0*k-1.0d0)*x*q1-(k-1.0d0)*q0)/k
              qn(k)=qf
              qd(k)=(qn(k-1)-x*qf)*k/(1.0d0-x*x)
              q0=q1
15            q1=qf
        else
           qc2=1.0d0/x
           do 20 j=1,n
              qc2=qc2*j/((2.0*j+1.0d0)*x)
              if (j.eq.n-1) qc1=qc2
20         continue
           do 35 l=0,1
              nl=n+l
              qf=1.0d0
              qr=1.0d0
              do 25 k=1,500
                 qr=qr*(0.5d0*nl+k-1.0d0)*(0.5d0*(nl-1)+k)
     &              /((nl+k-0.5d0)*k*x*x)
                 qf=qf+qr
                 if (dabs(qr/qf).lt.eps) go to 30
25            continue
30            if (l.eq.0) then
                 qn(n-1)=qf*qc1
              else
                 qn(n)=qf*qc2
              endif
35         continue
           qf2=qn(n)
           qf1=qn(n-1)
           do 40 k=n,2,-1
              qf0=((2*k-1.0d0)*x*qf1-k*qf2)/(k-1.0d0)
              qn(k-2)=qf0
              qf2=qf1
40            qf1=qf0
           qd(0)=1.0d0/(1.0d0-x*x)
           do 45 k=1,n
45            qd(k)=k*(qn(k-1)-x*qn(k))/(1.0d0-x*x)
        endif
        return
        end
