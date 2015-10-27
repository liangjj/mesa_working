        program msegv
c
c       ============================================================
c       purpose: this program computes a sequence of characteristic 
c                values of spheroidal prolate and oblate functions 
c                using subroutine segv
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                kd --- function code
c                       kd=1 for prolate; kd=-1 for oblate
c       output:  cv --- characteristic value for given m, n and c
c                eg(l) --- characteristic value for mode m and n'
c                          ( l = n' - m + 1 )
c       examples:
c                prolate: ( kd = 1 , m = 1, n = 5, c = 5.0 )
c
c                m      n       c        lambda mn(c)
c              ---------------------------------------
c                1      1      5.0        5.35042230
c                1      2      5.0       14.64295624
c                1      3      5.0       23.39761312
c                1      4      5.0       32.42194359
c                1      5      5.0       42.65818215
c
c                oblate: ( kd = -1 , m = 1, n = 5, c = 5.0 )
c
c                m      n       c      lambda mn(-ic)
c               --------------------------------------
c                1      1      5.0       -7.49338828
c                1      2      5.0       -7.12783752
c                1      3      5.0        2.75036721
c                1      4      5.0        8.69495925
c                1      5      5.0       18.43931577
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension eg(100)
        write(*,*)'please enter kd, m, n and c '
        read(*,*)kd,m,n,c
        write(*,15)kd,m,n,c
        write(*,*)
        call segv(m,n,c,kd,cv,eg)
        if (kd.eq.1) then
           write(*,*)'  m      n       c       lambda mn(c)'
        else if (kd.eq.-1) then
           write(*,*)'  m      n       c      lambda mn(-ic)'
        endif
        write(*,*)'---------------------------------------'
        do 10 l=1,n-m+1
           n1=m+l-1
10         write(*,20)m,n1,c,eg(l)
15      format(1x,'kd =',i2,',  m =',i3,',   n =',i3,',  c =',f5.1)
20      format(1x,i3,4x,i3,4x,f5.1,f18.8)
        end


        subroutine segv(m,n,c,kd,cv,eg)
c
c       =========================================================
c       purpose: compute the characteristic values of spheroidal
c                wave functions
c       input :  m  --- mode parameter
c                n  --- mode parameter
c                c  --- spheroidal parameter
c                kd --- function code
c                       kd=1 for prolate; kd=-1 for oblate
c       output:  cv --- characteristic value for given m, n and c
c                eg(l) --- characteristic value for mode m and n'
c                          ( l = n' - m + 1 )
c       =========================================================
c
        implicit double precision (a-h,o-z)
        dimension b(100),h(100),d(300),e(300),f(300),cv0(100),
     &            a(300),g(300),eg(200)
        if (c.lt.1.0d-10) then
           do 5 i=1,n
5             eg(i)=(i+m)*(i+m-1.0d0)
           go to 70
        endif                                           
        icm=(n-m+2)/2
        nm=10+int(0.5*(n-m)+c)
        cs=c*c*kd
        do 60 l=0,1
           do 10 i=1,nm
              if (l.eq.0) k=2*(i-1)
              if (l.eq.1) k=2*i-1
              dk0=m+k
              dk1=m+k+1
              dk2=2*(m+k)
              d2k=2*m+k
              a(i)=(d2k+2.0)*(d2k+1.0)/((dk2+3.0)*(dk2+5.0))*cs
              d(i)=dk0*dk1+(2.0*dk0*dk1-2.0*m*m-1.0)/((dk2-1.0)
     &             *(dk2+3.0))*cs
10            g(i)=k*(k-1.0)/((dk2-3.0)*(dk2-1.0))*cs
           do 15 k=2,nm
              e(k)=dsqrt(a(k-1)*g(k))
15            f(k)=e(k)*e(k)
           f(1)=0.0d0
           e(1)=0.0d0
           xa=d(nm)+dabs(e(nm))
           xb=d(nm)-dabs(e(nm))
           nm1=nm-1
           do 20 i=1,nm1
              t=dabs(e(i))+dabs(e(i+1))
              t1=d(i)+t
              if (xa.lt.t1) xa=t1
              t1=d(i)-t
              if (t1.lt.xb) xb=t1
20         continue
           do 25 i=1,icm
              b(i)=xa
25            h(i)=xb
           do 55 k=1,icm
              do 30 k1=k,icm
                 if (b(k1).lt.b(k)) then
                    b(k)=b(k1)
                    go to 35
                 endif
30            continue
35            if (k.ne.1.and.h(k).lt.h(k-1)) h(k)=h(k-1)
40            x1=(b(k)+h(k))/2.0d0
              cv0(k)=x1
              if (dabs((b(k)-h(k))/x1).lt.1.0d-14) go to 50
              j=0
              s=1.0d0
              do 45 i=1,nm
                 if (s.eq.0.0d0) s=s+1.0d-30
                 t=f(i)/s
                 s=d(i)-t-x1
                 if (s.lt.0.0d0) j=j+1
45            continue
              if (j.lt.k) then
                 h(k)=x1
              else
                 b(k)=x1
                 if (j.ge.icm) then
                    b(icm)=x1
                 else
                    if (h(j+1).lt.x1) h(j+1)=x1
                    if (x1.lt.b(j)) b(j)=x1
                 endif
              endif
              go to 40
50            cv0(k)=x1
              if (l.eq.0) eg(2*k-1)=cv0(k)
              if (l.eq.1) eg(2*k)=cv0(k)
55         continue
60      continue
70      cv=eg(n-m+1)
        return
        end
