        program mlqna
c
c       ======================================================
c       purpose: this program computes the legendre functions
c                qn(x) and qn'(x) using subroutine lqna
c       input :  x  --- argument of qn(x)  ( -1 ף x ף 1 )
c                n  --- degree of qn(x)  ( n = 0,1,תתת )
c       output:  qn(n) --- qn(x)
c                qd(n) --- qn'(x)
c       example:  x = 0.50
c                 n        qn(x)         qn'(x)
c                ---------------------------------
c                 0      .54930614     1.33333333
c                 1     -.72534693     1.21597281
c                 2     -.81866327     -.84270745
c                 3     -.19865477    -2.87734353
c                 4      .44017453    -2.23329085
c                 5      .55508089     1.08422720
c       ======================================================
c
        double precision qn,qd,x
        dimension qn(0:100),qd(0:100)
        write(*,*)'  please enter nmax and x'
        read(*,*)n,x
        write(*,30) x
        write(*,*)
        call lqna(n,x,qn,qd)
        write(*,*)'  n        qn(x)         qn''(x)'
        write(*,*)' ---------------------------------'
        do 10 k=0,n
10         write(*,20)k,qn(k),qd(k)
20      format(1x,i3,2f15.8)
30      format(3x,'x =',f5.2)
        end


        subroutine lqna(n,x,qn,qd)
c
c       =====================================================
c       purpose: compute legendre functions qn(x) and qn'(x)
c       input :  x  --- argument of qn(x) ( -1 ף x ף 1 )
c                n  --- degree of qn(x) ( n = 0,1,2,תתת )
c       output:  qn(n) --- qn(x)
c                qd(n) --- qn'(x)
c                ( 1.0d+300 stands for infinity )
c       =====================================================
c
        implicit double precision (q,x)
        dimension qn(0:n),qd(0:n)
        if (dabs(x).eq.1.0d0) then
           do 10 k=0,n
              qn(k)=1.0d+300
              qd(k)=-1.0d+300
10         continue
        else if (dabs(x).lt.1.0d0) then
           q0=0.5d0*dlog((1.0d0+x)/(1.0d0-x))
           q1=x*q0-1.0d0
           qn(0)=q0
           qn(1)=q1
           qd(0)=1.0d0/(1.0d0-x*x)
           qd(1)=qn(0)+x*qd(0)
           do 15 k=2,n
              qf=((2*k-1)*x*q1-(k-1)*q0)/k
              qn(k)=qf
              qd(k)=(qn(k-1)-x*qf)*k/(1.0d0-x*x)
              q0=q1
15            q1=qf
        endif
        return
        end
