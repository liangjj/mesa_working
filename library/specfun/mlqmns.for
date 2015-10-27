        program mlqmns
c
c       =========================================================
c       purpose: this program computes the associated legendre 
c                functions qmn(x) and their derivatives qmn'(x) 
c                for a given order using subroutine lqmns
c       input :  x --- argument of qmn(x)
c                m --- order of qmn(x),  m = 0,1,2,...
c                n --- degree of qmn(x), n = 0,1,2,...
c       output:  qm(n) --- qmn(x)
c                qd(n) --- qmn'(x)
c       examples:
c                m = 1,  n = 5,  x = .5
c                n        qmn(x)           qmn'(x)
c               -------------------------------------
c                0    .11547005d+01    .76980036d+00
c                1    .10530633d+01    .23771592d+01
c                2   -.72980606d+00    .51853281d+01
c                3   -.24918526d+01    .10914062d+01
c                4   -.19340866d+01   -.11454786d+02
c                5    .93896830d+00   -.18602587d+02
c
c                m = 2,  n = 5,  x = 2.5
c                n        qmn(x)           qmn'(x)
c               -------------------------------------
c                0    .95238095d+00   -.52607710d+00
c                1    .38095238d+00   -.36281179d+00
c                2    .12485160d+00   -.17134314d+00
c                3    .36835513d-01   -.66284127d-01
c                4    .10181730d-01   -.22703958d-01
c                5    .26919481d-02   -.71662396d-02
c       =========================================================
c
        implicit double precision (q,x,y)
        dimension qm(0:200),qd(0:200)
        write(*,*)'please enter m, n, and x '
        read(*,*)m,n,x
        write(*,30)m,n,x
        call lqmns(m,n,x,qm,qd)
        write(*,*)
        write(*,*)'  n        qmn(x)           qmn''(x)'
        write(*,*)' -------------------------------------'
        do 10 j=0,n
        write(*,20)j,qm(j),qd(j)
10      continue
20      format(1x,i3,2d17.8)
30      format(1x,'m =',i2,',  ','n =',i2,',  ','x =',f5.1)
        end


        subroutine lqmns(m,n,x,qm,qd)
c
c       ========================================================
c       purpose: compute associated legendre functions qmn(x)
c                and qmn'(x) for a given order
c       input :  x --- argument of qmn(x)
c                m --- order of qmn(x),  m = 0,1,2,...
c                n --- degree of qmn(x), n = 0,1,2,...
c       output:  qm(n) --- qmn(x)
c                qd(n) --- qmn'(x)
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension qm(0:n),qd(0:n)
        do 10 k=0,n
           qm(k)=0.0d0
10         qd(k)=0.0d0
        if (dabs(x).eq.1.0d0) then
           do 15 k=0,n
              qm(k)=1.0d+300
15            qd(k)=1.0d+300
           return
        endif
        ls=1
        if (dabs(x).gt.1.0d0) ls=-1
        xq=dsqrt(ls*(1.0d0-x*x))
        q0=0.5d0*dlog(dabs((x+1.0)/(x-1.0)))
        q00=q0
        q10=-1.0d0/xq
        q01=x*q0-1.0d0
        q11=-ls*xq*(q0+x/(1.0d0-x*x))
        qf0=q00
        qf1=q10
        do 20 k=2,m
           qm0=-2.0d0*(k-1.0)/xq*x*qf1-ls*(k-1.0)*(2.0-k)*qf0
           qf0=qf1
20         qf1=qm0
        if (m.eq.0) qm0=q00
        if (m.eq.1) qm0=q10
        qm(0)=qm0
        if (dabs(x).lt.1.0001d0) then
           if (m.eq.0.and.n.gt.0) then
              qf0=q00
              qf1=q01
              do 25 k=2,n
                 qf2=((2.0*k-1.0d0)*x*qf1-(k-1.0)*qf0)/k
                 qm(k)=qf2
                 qf0=qf1
25               qf1=qf2
           endif
           qg0=q01
           qg1=q11
           do 30 k=2,m
              qm1=-2.0d0*(k-1.0)/xq*x*qg1-ls*k*(3.0-k)*qg0
              qg0=qg1
30            qg1=qm1
           if (m.eq.0) qm1=q01
           if (m.eq.1) qm1=q11
           qm(1)=qm1
           if (m.eq.1.and.n.gt.1) then
              qh0=q10
              qh1=q11
              do 35 k=2,n
                 qh2=((2.0*k-1.0d0)*x*qh1-k*qh0)/(k-1.0)
                 qm(k)=qh2
                 qh0=qh1
35               qh1=qh2
           else if (m.ge.2) then
              qg0=q00
              qg1=q01
              qh0=q10
              qh1=q11
              do 45 l=2,n
                 q0l=((2.0d0*l-1.0d0)*x*qg1-(l-1.0d0)*qg0)/l
                 q1l=((2.0*l-1.0d0)*x*qh1-l*qh0)/(l-1.0d0)
                 qf0=q0l
                 qf1=q1l
                 do 40 k=2,m
                    qmk=-2.0d0*(k-1.0)/xq*x*qf1-ls*(k+l-1.0)*
     &                  (l+2.0-k)*qf0
                    qf0=qf1
40                  qf1=qmk
                 qm(l)=qmk
                 qg0=qg1
                 qg1=q0l
                 qh0=qh1
45               qh1=q1l
           endif
        else
           if (dabs(x).gt.1.1) then
              km=40+m+n
           else
              km=(40+m+n)*int(-1.0-1.8*log(x-1.0))
           endif
           qf2=0.0d0
           qf1=1.0d0
           do 50 k=km,0,-1
              qf0=((2.0*k+3.0d0)*x*qf1-(k+2.0-m)*qf2)/(k+m+1.0)
              if (k.le.n) qm(k)=qf0
              qf2=qf1
50            qf1=qf0
           do 55 k=0,n
55            qm(k)=qm(k)*qm0/qf0
        endif
        if (dabs(x).lt.1.0d0) then
           do 60 k=0,n
60            qm(k)=(-1)**m*qm(k)
        endif
        qd(0)=((1.0d0-m)*qm(1)-x*qm(0))/(x*x-1.0)
        do 65 k=1,n
65         qd(k)=(k*x*qm(k)-(k+m)*qm(k-1))/(x*x-1.0)
        return
        end
