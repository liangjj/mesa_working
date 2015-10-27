        program mcva1
c
c       ============================================================
c       purpose: this program computes a sequence of characteristic 
c                values of mathieu functions using subroutine cva1
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c                kd --- case code
c                       kd=1 for cem(x,q)  ( m = 0,2,4,...)
c                       kd=2 for cem(x,q)  ( m = 1,3,5,...)
c                       kd=3 for sem(x,q)  ( m = 1,3,5,...)
c                       kd=4 for sem(x,q)  ( m = 2,4,6,...)
c       output:  cv(i) --- characteristic values; i = 1,2,3,...
c                for kd=1, cv(1), cv(2), cv(3),..., correspond to
c                the characteristic values of cem for m = 0,2,4,...
c                for kd=2, cv(1), cv(2), cv(3),..., correspond to
c                the characteristic values of cem for m = 1,3,5,...
c                for kd=3, cv(1), cv(2), cv(3),..., correspond to
c                the characteristic values of sem for m = 1,3,5,...
c                for kd=4, cv(1), cv(2), cv(3),..., correspond to
c                the characteristic values of sem for m = 0,2,4,...
c
c       example: mmax = 12,    q = 25.00
c
c                characteristic values of mathieu functions
c
c                  m            a                  b
c                ------------------------------------------
c                  0      -40.256779547
c                  1      -21.314899691      -40.256778985
c                  2       -3.522164727      -21.314860622
c                  3       12.964079444       -3.520941527
c                  4       27.805240581       12.986489953
c                  5       40.050190986       28.062765899
c                  6       48.975786716       41.801071292
c                  7       57.534689001       55.002957151
c                  8       69.524065166       69.057988351
c                  9       85.076999882       85.023356505
c                 10      103.230204804      103.225680042
c                 11      123.643012376      123.642713667
c                 12      146.207690643      146.207674647
c       ============================================================
c
        implicit double precision (a-h,o-z)
        dimension cv1(200),cv2(200),cve(200),cvs(200)
        write(*,*)'please enter mmax,q =?'
        read(*,*)mmax,q
        write(*,25)mmax,q
        write(*,*)
        call cva1(1,mmax,q,cv1)
        call cva1(2,mmax,q,cv2)
        do 10 j=1,mmax/2+1
           cve(2*j-1)=cv1(j)
10         cve(2*j)=cv2(j)
        call cva1(3,mmax,q,cv1)
        call cva1(4,mmax,q,cv2)
        do 15 j=1,mmax/2+1
           cvs(2*j)=cv1(j)
15         cvs(2*j+1)=cv2(j)
        write(*,35)
        write(*,*)
        write(*,*)'  m            a                  b'
        write(*,*)'------------------------------------------'
        do 20 j=0,mmax
           if (j.eq.0) write(*,30)j,cve(j+1)
           if (j.ne.0) write(*,30)j,cve(j+1),cvs(j+1)
20      continue
25      format(3x,6hmmax =,i3,',    ',3hq =,f6.2)
30      format(1x,i3,2f19.9)
35      format(1x,'characteristic values of mathieu functions')
        end


        subroutine cva1(kd,m,q,cv)
c
c       ============================================================
c       purpose: compute a sequence of characteristic values of
c                mathieu functions 
c       input :  m  --- maximum order of mathieu functions
c                q  --- parameter of mathieu functions
c                kd --- case code
c                       kd=1 for cem(x,q)  ( m = 0,2,4,תתת )
c                       kd=2 for cem(x,q)  ( m = 1,3,5,תתת )
c                       kd=3 for sem(x,q)  ( m = 1,3,5,תתת )
c                       kd=4 for sem(x,q)  ( m = 2,4,6,תתת )
c       output:  cv(i) --- characteristic values; i = 1,2,3,...
c                for kd=1, cv(1), cv(2), cv(3),..., correspond to
c                the characteristic values of cem for m = 0,2,4,...
c                for kd=2, cv(1), cv(2), cv(3),..., correspond to
c                the characteristic values of cem for m = 1,3,5,...
c                for kd=3, cv(1), cv(2), cv(3),..., correspond to
c                the characteristic values of sem for m = 1,3,5,...
c                for kd=4, cv(1), cv(2), cv(3),..., correspond to
c                the characteristic values of sem for m = 0,2,4,...
c       ============================================================
c
        implicit double precision (a-h,o-z)
        dimension g(200),h(200),d(500),e(500),f(500),cv(200)
        eps=1.0d-14
        icm=int(m/2)+1
        if (kd.eq.4) icm=m/2
        if (q.eq.0.0d0) then
           if (kd.eq.1) then
              do 10 ic=1,icm
10               cv(ic)=4.0d0*(ic-1.0d0)**2
           else if (kd.ne.4) then
              do 15 ic=1,icm
15               cv(ic)=(2.0d0*ic-1.0d0)**2
           else
              do 20 ic=1,icm
20               cv(ic)=4.0d0*ic*ic
           endif
        else
           nm=int(10+1.5*m+0.5*q)
           e(1)=0.0d0
           f(1)=0.0d0
           if (kd.eq.1) then
              d(1)=0.0d0
              do 25 i=2,nm
                 d(i)=4.0d0*(i-1.0d0)**2
                 e(i)=q
25               f(i)=q*q
              e(2)=dsqrt(2.0d0)*q
              f(2)=2.0d0*q*q
           else if (kd.ne.4) then
              d(1)=1.0d0+(-1)**kd*q
              do 30 i=2,nm
                 d(i)=(2.0d0*i-1.0d0)**2
                 e(i)=q
30               f(i)=q*q
           else
              d(1)=4.0d0
              do 35 i=2,nm
                 d(i)=4.0d0*i*i
                 e(i)=q
35               f(i)=q*q
           endif
           xa=d(nm)+dabs(e(nm))
           xb=d(nm)-dabs(e(nm))
           nm1=nm-1
           do 40 i=1,nm1
              t=dabs(e(i))+dabs(e(i+1))
              t1=d(i)+t
              if (xa.lt.t1) xa=t1
              t1=d(i)-t
              if (t1.lt.xb) xb=t1
40         continue
           do 45 i=1,icm
              g(i)=xa
45            h(i)=xb
           do 75 k=1,icm
              do 50 k1=k,icm
                 if (g(k1).lt.g(k)) then
                    g(k)=g(k1)
                    go to 55
                 endif
50            continue
55            if (k.ne.1.and.h(k).lt.h(k-1)) h(k)=h(k-1)
60            x1=(g(k)+h(k))/2.0d0
              cv(k)=x1
              if (dabs((g(k)-h(k))/x1).lt.eps) go to 70
              j=0
              s=1.0d0
              do 65 i=1,nm
                 if (s.eq.0.0d0) s=s+1.0d-30
                 t=f(i)/s
                 s=d(i)-t-x1
                 if (s.lt.0.0) j=j+1
65            continue
              if (j.lt.k) then
                 h(k)=x1
              else
                 g(k)=x1
                 if (j.ge.icm) then
                    g(icm)=x1
                 else
                    if (h(j+1).lt.x1) h(j+1)=x1
                    if (x1.lt.g(j)) g(j)=x1
                 endif
              endif
              go to 60
70            cv(k)=x1
75         continue
        endif
        return
        end
