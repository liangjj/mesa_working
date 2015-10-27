        program mlpmns
c
c       ========================================================
c       purpose: this program computes the associated legendre 
c                functions pmn(x) and their derivatives pmn'(x) 
c                for a given order using subroutine lpmns
c       input :  x --- argument of pmn(x)
c                m --- order of pmn(x),  m = 0,1,2,...,n
c                n --- degree of pmn(x), n = 0,1,2,...,n
c       output:  pm(n) --- pmn(x)
c                pd(n) --- pmn'(x)
c       examples:
c                m = 1,  n = 5,  x = .5
c                n        pmn(x)           pmn'(x)
c               -------------------------------------
c                0    .00000000d+00    .00000000d+00
c                1    .86602540d+00   -.57735027d+00
c                2    .12990381d+01    .17320508d+01
c                3    .32475953d+00    .62786842d+01
c                4   -.13531647d+01    .57735027d+01
c                5   -.19282597d+01   -.43977853d+01
c
c                m = 2,  n = 6,  x = 2.5
c                n        pmn(x)           pmn'(x)
c               -------------------------------------
c                0    .00000000d+00    .00000000d+00
c                1    .00000000d+00    .00000000d+00
c                2    .15750000d+02    .15000000d+02
c                3    .19687500d+03    .26625000d+03
c                4    .16832813d+04    .29812500d+04
c                5    .12230859d+05    .26876719d+05
c                6    .81141416d+05    .21319512d+06
c       =======================================================
c
        implicit double precision (p,x,y)
        dimension pm(0:200),pd(0:200)
        write(*,*)'please enter m, n, and x '
        read(*,*)m,n,x
        write(*,30)m,n,x
        call lpmns(m,n,x,pm,pd)
        write(*,*)
        write(*,*)'  n        pmn(x)           pmn''(x)    '
        write(*,*)' -------------------------------------'
        do 10 j=0,n
        write(*,20)j,pm(j),pd(j)
10      continue
20      format(1x,i3,2d17.8)
30      format(1x,'m =',i2,',  ','n =',i2,',  ','x =',f5.1)
        end


        subroutine lpmns(m,n,x,pm,pd)
c
c       ========================================================
c       purpose: compute associated legendre functions pmn(x)
c                and pmn'(x) for a given order
c       input :  x --- argument of pmn(x)
c                m --- order of pmn(x),  m = 0,1,2,...,n
c                n --- degree of pmn(x), n = 0,1,2,...,n
c       output:  pm(n) --- pmn(x)
c                pd(n) --- pmn'(x)
c       ========================================================
c
        implicit double precision (a-h,o-z)
        dimension pm(0:n),pd(0:n)
        do 10 k=0,n
           pm(k)=0.0d0
10         pd(k)=0.0d0
        if (dabs(x).eq.1.0d0) then
           do 15 k=0,n
              if (m.eq.0) then
                 pm(k)=1.0d0
                 pd(k)=0.5d0*k*(k+1.0)
                 if (x.lt.0.0) then
                    pm(k)=(-1)**k*pm(k)
                    pd(k)=(-1)**(k+1)*pd(k)
                 endif
              else if (m.eq.1) then
                 pd(k)=1.0d+300
              else if (m.eq.2) then
                 pd(k)=-0.25d0*(k+2.0)*(k+1.0)*k*(k-1.0)
                 if (x.lt.0.0) pd(k)=(-1)**(k+1)*pd(k)
              endif
15         continue
           return
        endif
        x0=dabs(1.0d0-x*x)
        pm0=1.0d0
        pmk=pm0
        do 20 k=1,m
           pmk=(2.0d0*k-1.0d0)*dsqrt(x0)*pm0
20         pm0=pmk
        pm1=(2.0d0*m+1.0d0)*x*pm0
        pm(m)=pmk
        pm(m+1)=pm1
        do 25 k=m+2,n
           pm2=((2.0d0*k-1.0d0)*x*pm1-(k+m-1.0d0)*pmk)/(k-m)
           pm(k)=pm2
           pmk=pm1
25         pm1=pm2
        pd(0)=((1.0d0-m)*pm(1)-x*pm(0))/(x*x-1.0)  
        do 30 k=1,n
30          pd(k)=(k*x*pm(k)-(k+m)*pm(k-1))/(x*x-1.0d0)
        return
        end
