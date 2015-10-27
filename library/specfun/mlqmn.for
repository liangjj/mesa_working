        program mlqmn
c
c       ===============================================================
c       purpose: this program computes the associated legendre  
c                functions qmn(x) and their derivatives qmn'(x) using
c                subroutine lqmn
c       input :  x --- argument of qmn(x) 
c                m --- order of qmn(x)  ( m = 0,1,2,תתת )
c                n --- degree of qmn(x) ( n = 0,1,2,תתת )
c       output:  qm(m,n) --- qmn(x)
c                qd(m,n) --- qmn'(x)
c       examples:
c
c       qmn(x):  x = 0.5
c       n\m      0           1           2           3           4
c       ---------------------------------------------------------------
c        0     .549306   -1.154701    1.333333   -5.388603   26.666667
c        1    -.725347   -1.053063    2.666667   -6.158403   32.000000
c        2    -.818663     .729806    4.069272  -12.316806   42.666667
c        3    -.198655    2.491853    -.493486  -23.778868   85.333333
c        4     .440175    1.934087  -11.036781   -9.325204  186.818394
c
c       qmn'(x): x = 0.5
c       n\m      0           1           2           3           4
c       ---------------------------------------------------------------
c        0    1.333333    -.769800    4.444444  -20.014809  145.777778
c        1    1.215973   -2.377159    3.555556  -24.633611  156.444444
c        2    -.842707   -5.185328    8.796526  -24.633611  199.111111
c        3   -2.877344   -1.091406   28.115454  -50.976710  227.555556
c        4   -2.233291   11.454786   25.483527 -197.068892  412.039838
c
c       qmn(x): x = 2.0
c       n\m      0           1           2           3           4
c       ---------------------------------------------------------------
c        0     .549306    -.577350    1.333333   -5.003702   26.666667
c        1     .098612    -.203274     .666667   -3.079201   18.666667
c        2     .021184    -.064946     .277089   -1.539601   10.666667
c        3     .004871    -.019817     .104220    -.679543    5.333333
c        4     .001161    -.005887     .036816    -.276005    2.427640
c
c       qmn'(x): x = 2.0
c       n\m      0           1           2           3           4
c       ---------------------------------------------------------------
c        0    -.333333     .384900   -1.111111    5.388603  -36.444444
c        1    -.117361     .249384    -.888889    4.618802  -32.000000
c        2    -.037496     .116680    -.519437    3.079201  -23.111111
c        3    -.011442     .046960    -.253375    1.720114  -14.222222
c        4    -.003399     .017331    -.110263     .849589   -7.748516
c       ===============================================================
c
        implicit double precision (q,x)
        dimension qm(0:100,0:100),qd(0:100,0:100)
        write(*,*)'  please enter m, n and x'
        read(*,*) m,n,x
        write(*,*)
        write(*,*)'  m     n      x          qmn(x)         qmn''(x)'
        write(*,*)' ---------------------------------------------------'
        call lqmn(100,m,n,x,qm,qd)
        do 15 j=0,n
           write(*,10)m,j,x,qm(m,j),qd(m,j)
15      continue
10      format(1x,i3,3x,i3,3x,f5.1,2d17.8)
        end


        subroutine lqmn(mm,m,n,x,qm,qd)
c
c       ==========================================================
c       purpose: compute the associated legendre functions of the
c                second kind, qmn(x) and qmn'(x)
c       input :  x  --- argument of qmn(x) 
c                m  --- order of qmn(x)  ( m = 0,1,2,תתת )
c                n  --- degree of qmn(x) ( n = 0,1,2,תתת )
c                mm --- physical dimension of qm and qd
c       output:  qm(m,n) --- qmn(x)
c                qd(m,n) --- qmn'(x)
c       ==========================================================
c
        implicit double precision (q,x)
        dimension qm(0:mm,0:n),qd(0:mm,0:n)
        if (dabs(x).eq.1.0d0) then
           do 10 i=0,m
           do 10 j=0,n
              qm(i,j)=1.0d+300
              qd(i,j)=1.0d+300
10         continue
           return
        endif
        ls=1
        if (dabs(x).gt.1.0d0) ls=-1
        xs=ls*(1.0d0-x*x)
        xq=dsqrt(xs)
        q0=0.5d0*dlog(dabs((x+1.0d0)/(x-1.0d0)))
        if (dabs(x).lt.1.0001d0) then
           qm(0,0)=q0
           qm(0,1)=x*q0-1.0d0
           qm(1,0)=-1.0d0/xq
           qm(1,1)=-xq*(q0+x/(1.0d0-x*x))
           do 15 i=0,1
           do 15 j=2,n
              qm(i,j)=((2.0d0*j-1.0d0)*x*qm(i,j-1)
     &               -(j+i-1.0d0)*qm(i,j-2))/(j-i)
15         continue
           do 20 j=0,n
           do 20 i=2,m
              qm(i,j)=-2.0d0*(i-1.0d0)*x/xq*qm(i-1,j)-ls*
     &                (j+i-1.0d0)*(j-i+2.0d0)*qm(i-2,j)
20         continue
        else
           if (dabs(x).gt.1.1d0) then
              km=40+m+n
           else
              km=(40+m+n)*int(-1.0-1.8*log(x-1.0))
           endif
           qf2=0.0d0
           qf1=1.0d0
           do 25 k=km,0,-1
              qf0=((2*k+3.0d0)*x*qf1-(k+2.0d0)*qf2)/(k+1.0d0)
              if (k.le.n) qm(0,k)=qf0
              qf2=qf1
25            qf1=qf0
           do 30 k=0,n
30            qm(0,k)=q0*qm(0,k)/qf0
           qf2=0.0d0
           qf1=1.0d0
           do 35 k=km,0,-1
              qf0=((2*k+3.0d0)*x*qf1-(k+1.0d0)*qf2)/(k+2.0d0)
              if (k.le.n) qm(1,k)=qf0
              qf2=qf1
35            qf1=qf0
           q10=-1.0d0/xq
           do 40 k=0,n
40            qm(1,k)=q10*qm(1,k)/qf0
           do 45 j=0,n
              q0=qm(0,j)
              q1=qm(1,j)
              do 45 i=0,m-2
                 qf=-2.0d0*(i+1)*x/xq*q1+(j-i)*(j+i+1.0d0)*q0
                 qm(i+2,j)=qf
                 q0=q1
                 q1=qf
45         continue
        endif
        qd(0,0)=ls/xs
        do 50 j=1,n
50         qd(0,j)=ls*j*(qm(0,j-1)-x*qm(0,j))/xs
        do 55 j=0,n
        do 55 i=1,m
           qd(i,j)=ls*i*x/xs*qm(i,j)+(i+j)*(j-i+1.0d0)/xq*qm(i-1,j)
55      continue
        return
        end
