        program mmtu0
c
c       ============================================================
c       purpose: this program computes mathieu functions cem(x,q), 
c                sem(x,q) and their derivatives using subroutine
c                mtu0 ( q ò 0 )
c       input :  kf  --- function code
c                        kf=1 for computing cem(x,q) and cem'(x,q)
c                        kf=2 for computing sem(x,q) and sem'(x,q)
c                m   --- order of mathieu functions
c                q   --- parameter of mathieu functions
c                x   --- argument of mathieu functions (in degrees)
c       output:  csf --- cem(x,q) or sem(x,q)
c                csd --- cem'x,q) or sem'x,q)
c       example: x = 40
c           m     q    cem(x,q)   cem'(x,q)    sem(x,q)  sem'(x,q)
c          --------------------------------------------------------
c           0    5.0   .3025683    .9470247
c           1    5.0   .7669652   1.2873097    .2988052   .9606824
c           2    5.0   .9102723   -.3463855    .7549264  1.4743128
c           5    5.0  -.9810931   -.6328576    .1694850 -4.8676455
c           0   25.0   .0515371    .3823737
c           1   25.0   .2074402   1.2646301    .0515365   .3823777
c           2   25.0  -.5297051  -2.4292679    .2074275  1.2646996
c           5   25.0   .7507159  -3.9047012   1.1881232   .3258081
c       ============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter kf, m, q and x (in degrees)'
        read(*,*)kf,m,q,x
        write(*,10)kf,m,q,x
        write(*,*)
        if (kf.eq.1) write(*,*)' x(degs.)    cem(x,q)       cem''(x,q)'
        if (kf.eq.2) write(*,*)' x(degs.)    sem(x,q)       sem''(x,q)'
        write(*,*)' --------------------------------------'
        call mtu0(kf,m,q,x,csf,csd)
        write(*,20)x,csf,csd
10      format(1x,4hkf =,i2,',  ',3hm =,i2,',  ',
     &         3hq =,f5.1,',  ',3hx =,f5.1)
20      format(2x,f5.1,2f16.9)
        end


        subroutine mtu0(kf,m,q,x,csf,csd)
c
c       ===============================================================
c       purpose: compute mathieu functions cem(x,q) and sem(x,q)
c                and their derivatives ( q ò 0 )
c       input :  kf  --- function code
c                        kf=1 for computing cem(x,q) and cem'(x,q)
c                        kf=2 for computing sem(x,q) and sem'(x,q)
c                m   --- order of mathieu functions
c                q   --- parameter of mathieu functions
c                x   --- argument of mathieu functions (in degrees)
c       output:  csf --- cem(x,q) or sem(x,q)
c                csd --- cem'x,q) or sem'x,q)
c       routines called:
c            (1) cva2 for computing the characteristic values
c            (2) fcoef for computing the expansion coefficients
c       ===============================================================
c
        implicit double precision (a-h,o-z)
        dimension fg(251)
        eps=1.0d-14
        if (kf.eq.1.and.m.eq.2*int(m/2)) kd=1
        if (kf.eq.1.and.m.ne.2*int(m/2)) kd=2
        if (kf.eq.2.and.m.ne.2*int(m/2)) kd=3
        if (kf.eq.2.and.m.eq.2*int(m/2)) kd=4
        call cva2(kd,m,q,a)
        if (q.le.1.0d0) then
           qm=7.5+56.1*sqrt(q)-134.7*q+90.7*sqrt(q)*q
        else
           qm=17.0+3.1*sqrt(q)-.126*q+.0037*sqrt(q)*q
        endif
        km=int(qm+0.5*m)
        call fcoef(kd,m,q,a,fg)
        ic=int(m/2)+1
        rd=1.74532925199433d-2
        xr=x*rd
        csf=0.0d0
        do 10 k=1,km
           if (kd.eq.1) then
              csf=csf+fg(k)*dcos((2*k-2)*xr)
           else if (kd.eq.2) then
              csf=csf+fg(k)*dcos((2*k-1)*xr)
           else if (kd.eq.3) then
              csf=csf+fg(k)*dsin((2*k-1)*xr)
           else if (kd.eq.4) then
              csf=csf+fg(k)*dsin(2*k*xr)
           endif
           if (k.ge.ic.and.dabs(fg(k)).lt.dabs(csf)*eps) go to 15
10         continue
15      csd=0.0d0
        do 20 k=1,km
           if (kd.eq.1) then
              csd=csd-(2*k-2)*fg(k)*dsin((2*k-2)*xr)
           else if (kd.eq.2) then
              csd=csd-(2*k-1)*fg(k)*dsin((2*k-1)*xr)
           else if (kd.eq.3) then
              csd=csd+(2*k-1)*fg(k)*dcos((2*k-1)*xr)
           else if (kd.eq.4) then
              csd=csd+2.0d0*k*fg(k)*dcos(2*k*xr)
           endif
           if (k.ge.ic.and.dabs(fg(k)).lt.dabs(csd)*eps) go to 25
20         continue
25      return
        end


        subroutine fcoef(kd,m,q,a,fc)
c
c       =====================================================
c       purpose: compute expansion coefficients for mathieu
c                functions and modified mathieu functions
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c                kd --- case code
c                       kd=1 for cem(x,q)  ( m = 0,2,4,...)
c                       kd=2 for cem(x,q)  ( m = 1,3,5,...)
c                       kd=3 for sem(x,q)  ( m = 1,3,5,...)
c                       kd=4 for sem(x,q)  ( m = 2,4,6,...)
c                a  --- characteristic value of mathieu
c                       functions for given m and q
c       output:  fc(k) --- expansion coefficients of mathieu
c                       functions ( k= 1,2,...,km )
c                       fc(1),fc(2),fc(3),... correspond to
c                       a0,a2,a4,... for kd=1 case, a1,a3,
c                       a5,... for kd=2 case, b1,b3,b5,...
c                       for kd=3 case and b2,b4,b6,... for
c                       kd=4 case
c       =====================================================
c
        implicit double precision (a-h,o-z)
        dimension fc(251)
        if (q.le.1.0d0) then
           qm=7.5+56.1*sqrt(q)-134.7*q+90.7*sqrt(q)*q
        else
           qm=17.0+3.1*sqrt(q)-.126*q+.0037*sqrt(q)*q
        endif
        km=int(qm+0.5*m)
        if (q.eq.0.0d0) then
           do 10 k=1,km
10            fc(k)=0.0d0
           if (kd.eq.1) then
              fc((m+2)/2)=1.0d0
              if (m.eq.0) fc(1)=1.0d0/dsqrt(2.0d0)
           else if (kd.eq.4) then
              fc(m/2)=1.0d0
           else
              fc((m+1)/2)=1.0d0
           endif
           return
        endif
        kb=0
        s=0.0d0
        f=1.0d-100
        u=0.0d0
        fc(km)=0.0d0
        if (kd.eq.1) then
           do 25 k=km,3,-1
              v=u
              u=f
              f=(a-4.0d0*k*k)*u/q-v
              if (dabs(f).lt.dabs(fc(k+1))) then
                 kb=k
                 fc(1)=1.0d-100
                 sp=0.0d0
                 f3=fc(k+1)
                 fc(2)=a/q*fc(1)
                 fc(3)=(a-4.0d0)*fc(2)/q-2.0d0*fc(1)
                 u=fc(2)
                 f1=fc(3)
                 do 15 i=3,kb
                    v=u
                    u=f1
                    f1=(a-4.0d0*(i-1.0d0)**2)*u/q-v
                    fc(i+1)=f1
                    if (i.eq.kb) f2=f1
                    if (i.ne.kb) sp=sp+f1*f1
15               continue
                 sp=sp+2.0d0*fc(1)**2+fc(2)**2+fc(3)**2
                 ss=s+sp*(f3/f2)**2
                 s0=dsqrt(1.0d0/ss)
                 do 20 j=1,km
                    if (j.le.kb+1) then
                       fc(j)=s0*fc(j)*f3/f2
                    else
                       fc(j)=s0*fc(j)
                    endif
20               continue
                 go to 85
              else
                 fc(k)=f
                 s=s+f*f
              endif
25         continue
           fc(2)=q*fc(3)/(a-4.0d0-2.0d0*q*q/a)
           fc(1)=q/a*fc(2)
           s=s+2.0d0*fc(1)**2+fc(2)**2
           s0=dsqrt(1.0d0/s)
           do 30 k=1,km
30            fc(k)=s0*fc(k)
        else if (kd.eq.2.or.kd.eq.3) then
           do 35 k=km,3,-1
              v=u
              u=f
              f=(a-(2.0d0*k-1)**2)*u/q-v
              if (dabs(f).ge.dabs(fc(k))) then
                 fc(k-1)=f
                 s=s+f*f
              else
                 kb=k
                 f3=fc(k)
                 go to 45
              endif
35         continue
           fc(1)=q/(a-1.0d0-(-1)**kd*q)*fc(2)
           s=s+fc(1)*fc(1)
           s0=dsqrt(1.0d0/s)
           do 40 k=1,km
40            fc(k)=s0*fc(k)
           go to 85
45         fc(1)=1.0d-100
           fc(2)=(a-1.0d0-(-1)**kd*q)/q*fc(1)
           sp=0.0d0
           u=fc(1)
           f1=fc(2)
           do 50 i=2,kb-1
              v=u
              u=f1
              f1=(a-(2.0d0*i-1.0d0)**2)*u/q-v
              if (i.ne.kb-1) then
                 fc(i+1)=f1
                 sp=sp+f1*f1
              else
                 f2=f1
              endif
50         continue
           sp=sp+fc(1)**2+fc(2)**2
           ss=s+sp*(f3/f2)**2
           s0=1.0d0/dsqrt(ss)
           do 55 j=1,km
              if (j.lt.kb) fc(j)=s0*fc(j)*f3/f2
              if (j.ge.kb) fc(j)=s0*fc(j)
55         continue
        else if (kd.eq.4) then
           do 60 k=km,3,-1
              v=u
              u=f
              f=(a-4.0d0*k*k)*u/q-v
              if (dabs(f).ge.dabs(fc(k))) then
                 fc(k-1)=f
                 s=s+f*f
              else
                 kb=k
                 f3=fc(k)
                 go to 70
              endif
60         continue
           fc(1)=q/(a-4.0d0)*fc(2)
           s=s+fc(1)*fc(1)
           s0=dsqrt(1.0d0/s)
           do 65 k=1,km
65            fc(k)=s0*fc(k)
           go to 85
70         fc(1)=1.0d-100
           fc(2)=(a-4.0d0)/q*fc(1)
           sp=0.0d0
           u=fc(1)
           f1=fc(2)
           do 75 i=2,kb-1
              v=u
              u=f1
              f1=(a-4.0d0*i*i)*u/q-v
              if (i.ne.kb-1) then
                 fc(i+1)=f1
                 sp=sp+f1*f1
              else
                 f2=f1
              endif
75         continue
           sp=sp+fc(1)**2+fc(2)**2
           ss=s+sp*(f3/f2)**2
           s0=1.0d0/dsqrt(ss)
           do 80 j=1,km
              if (j.lt.kb) fc(j)=s0*fc(j)*f3/f2
              if (j.ge.kb) fc(j)=s0*fc(j)
80         continue
        endif
85      if (fc(1).lt.0.0d0) then
           do 90 j=1,km
90            fc(j)=-fc(j)
        endif
        return
        end


        subroutine cva2(kd,m,q,a)
c
c       ======================================================
c       purpose: calculate a specific characteristic value of
c                mathieu functions
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c                kd --- case code
c                       kd=1 for cem(x,q)  ( m = 0,2,4,...)
c                       kd=2 for cem(x,q)  ( m = 1,3,5,...)
c                       kd=3 for sem(x,q)  ( m = 1,3,5,...)
c                       kd=4 for sem(x,q)  ( m = 2,4,6,...)
c       output:  a  --- characteristic value
c       routines called:
c             (1) refine for finding accurate characteristic
c                 values using an iteration method
c             (2) cv0 for finding initial characteristic
c                 values using polynomial approximation
c             (3) cvqm for computing initial characteristic
c                 values for q ó 3*m
c             (3) cvql for computing initial characteristic
c                 values for q ò m*m
c       ======================================================
c
        implicit double precision (a-h,o-z)
        if (m.le.12.or.q.le.3.0*m.or.q.gt.m*m) then
            call cv0(kd,m,q,a)
            if (q.ne.0.0d0) call refine(kd,m,q,a,1)
        else
           ndiv=10
           delta=(m-3.0)*m/ndiv
           if ((q-3.0*m).le.(m*m-q)) then
5             nn=int((q-3.0*m)/delta)+1
              delta=(q-3.0*m)/nn
              q1=2.0*m
              call cvqm(m,q1,a1)
              q2=3.0*m
              call cvqm(m,q2,a2)
              qq=3.0*m
              do 10 i=1,nn
                 qq=qq+delta
                 a=(a1*q2-a2*q1+(a2-a1)*qq)/(q2-q1)
                 iflag=1
                 if (i.eq.nn) iflag=-1
                 call refine(kd,m,qq,a,iflag)
                 q1=q2
                 q2=qq
                 a1=a2
                 a2=a
10            continue
              if (iflag.eq.-10) then
                 ndiv=ndiv*2
                 delta=(m-3.0)*m/ndiv
                 go to 5
              endif
           else
15            nn=int((m*m-q)/delta)+1
              delta=(m*m-q)/nn
              q1=m*(m-1.0)
              call cvql(kd,m,q1,a1)
              q2=m*m
              call cvql(kd,m,q2,a2)
              qq=m*m
              do 20 i=1,nn
                 qq=qq-delta
                 a=(a1*q2-a2*q1+(a2-a1)*qq)/(q2-q1)
                 iflag=1
                 if (i.eq.nn) iflag=-1
                 call refine(kd,m,qq,a,iflag)
                 q1=q2
                 q2=qq
                 a1=a2
                 a2=a
20            continue
              if (iflag.eq.-10) then
                 ndiv=ndiv*2
                 delta=(m-3.0)*m/ndiv
                 go to 15
              endif
           endif
        endif
        return
        end


        subroutine refine(kd,m,q,a,iflag)
c
c       =====================================================
c       purpose: calculate the accurate characteristic value
c                by the secant method
c       input :  m --- order of mathieu functions
c                q --- parameter of mathieu functions
c                a --- initial characteristic value
c       output:  a --- refineed characteristic value
c       routine called:  cvf for computing the value of f for
c                        characteristic equation
c       ========================================================
c
        implicit double precision (a-h,o-z)
        eps=1.0d-14
        mj=10+m
        ca=a
        delta=0.0d0
        x0=a
        call cvf(kd,m,q,x0,mj,f0)
        x1=1.002*a
        call cvf(kd,m,q,x1,mj,f1)
5       do 10 it=1,100
           mj=mj+1
           x=x1-(x1-x0)/(1.0d0-f0/f1)
           call cvf(kd,m,q,x,mj,f)
           if (abs(1.0-x1/x).lt.eps.or.f.eq.0.0) go to 15
           x0=x1
           f0=f1
           x1=x
10         f1=f
15      a=x
        if (delta.gt.0.05) then
           a=ca
           if (iflag.lt.0) then
              iflag=-10
           endif
           return
        endif
        if (abs((a-ca)/ca).gt.0.05) then
           x0=ca
           delta=delta+0.005d0
           call cvf(kd,m,q,x0,mj,f0)
           x1=(1.0d0+delta)*ca
           call cvf(kd,m,q,x1,mj,f1)
           go to 5
        endif
        return
        end


        subroutine cvf(kd,m,q,a,mj,f)
c
c       ======================================================
c       purpose: compute the value of f for characteristic
c                equation of mathieu functions
c       input :  m --- order of mathieu functions
c                q --- parameter of mathieu functions
c                a --- characteristic value
c       output:  f --- value of f for characteristic equation
c       ======================================================
c
        implicit double precision (a-h,o-z)
        b=a
        ic=int(m/2)
        l=0
        l0=0
        j0=2
        jf=ic
        if (kd.eq.1) l0=2
        if (kd.eq.1) j0=3
        if (kd.eq.2.or.kd.eq.3) l=1
        if (kd.eq.4) jf=ic-1
        t1=0.0d0
        do 10 j=mj,ic+1,-1
10         t1=-q*q/((2.0d0*j+l)**2-b+t1)
        if (m.le.2) then
           t2=0.0d0
           if (kd.eq.1.and.m.eq.0) t1=t1+t1
           if (kd.eq.1.and.m.eq.2) t1=-2.0*q*q/(4.0-b+t1)-4.0
           if (kd.eq.2.and.m.eq.1) t1=t1+q
           if (kd.eq.3.and.m.eq.1) t1=t1-q
        else
           if (kd.eq.1) t0=4.0d0-b+2.0d0*q*q/b
           if (kd.eq.2) t0=1.0d0-b+q
           if (kd.eq.3) t0=1.0d0-b-q
           if (kd.eq.4) t0=4.0d0-b
           t2=-q*q/t0
           do 15 j=j0,jf
15            t2=-q*q/((2.0d0*j-l-l0)**2-b+t2)
        endif
        f=(2.0d0*ic+l)**2+t1+t2-b
        return
        end


        subroutine cv0(kd,m,q,a0)
c
c       =====================================================
c       purpose: compute the initial characteristic value of
c                mathieu functions for m ó 12  or q ó 300 or
c                q ò m*m
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c       output:  a0 --- characteristic value
c       routines called:
c             (1) cvqm for computing initial characteristic
c                 value for q ó 3*m
c             (2) cvql for computing initial characteristic
c                 value for q ò m*m
c       ====================================================
c
        implicit double precision (a-h,o-z)
        q2=q*q
        if (m.eq.0) then
           if (q.le.1.0) then
              a0=(((.0036392*q2-.0125868)*q2+.0546875)*q2-.5)*q2
           else if (q.le.10.0) then
              a0=((3.999267d-3*q-9.638957d-2)*q-.88297)*q
     &           +.5542818
           else
              call cvql(kd,m,q,a0)
           endif
        else if (m.eq.1) then
           if (q.le.1.0.and.kd.eq.2) then
              a0=(((-6.51e-4*q-.015625)*q-.125)*q+1.0)*q+1.0
           else if (q.le.1.0.and.kd.eq.3) then
              a0=(((-6.51e-4*q+.015625)*q-.125)*q-1.0)*q+1.0
           else if (q.le.10.0.and. kd.eq.2) then
              a0=(((-4.94603d-4*q+1.92917d-2)*q-.3089229)
     &           *q+1.33372)*q+.811752
           else if (q.le.10.0.and.kd.eq.3) then
              a0=((1.971096d-3*q-5.482465d-2)*q-1.152218)
     &           *q+1.10427
           else
              call cvql(kd,m,q,a0)
           endif
        else if (m.eq.2) then
           if (q.le.1.0.and.kd.eq.1) then
              a0=(((-.0036391*q2+.0125888)*q2-.0551939)*q2
     &           +.416667)*q2+4.0
           else if (q.le.1.0.and.kd.eq.4) then
              a0=(.0003617*q2-.0833333)*q2+4.0
           else if (q.le.15.and.kd.eq.1) then
              a0=(((3.200972d-4*q-8.667445d-3)*q
     &           -1.829032d-4)*q+.9919999)*q+3.3290504
           else if (q.le.10.0.and.kd.eq.4) then
              a0=((2.38446d-3*q-.08725329)*q-4.732542d-3)
     &           *q+4.00909
           else
              call cvql(kd,m,q,a0)
           endif
        else if (m.eq.3) then
           if (q.le.1.0.and.kd.eq.2) then
              a0=((6.348e-4*q+.015625)*q+.0625)*q2+9.0
           else if (q.le.1.0.and.kd.eq.3) then
              a0=((6.348e-4*q-.015625)*q+.0625)*q2+9.0
           else if (q.le.20.0.and.kd.eq.2) then
              a0=(((3.035731d-4*q-1.453021d-2)*q
     &           +.19069602)*q-.1039356)*q+8.9449274
           else if (q.le.15.0.and.kd.eq.3) then
              a0=((9.369364d-5*q-.03569325)*q+.2689874)*q
     &           +8.771735
           else
              call cvql(kd,m,q,a0)
           endif
        else if (m.eq.4) then
           if (q.le.1.0.and.kd.eq.1) then
              a0=((-2.1e-6*q2+5.012e-4)*q2+.0333333)*q2+16.0
           else if (q.le.1.0.and.kd.eq.4) then
              a0=((3.7e-6*q2-3.669e-4)*q2+.0333333)*q2+16.0
           else if (q.le.25.0.and.kd.eq.1) then
              a0=(((1.076676d-4*q-7.9684875d-3)*q
     &           +.17344854)*q-.5924058)*q+16.620847
           else if (q.le.20.0.and.kd.eq.4) then
              a0=((-7.08719d-4*q+3.8216144d-3)*q
     &           +.1907493)*q+15.744
           else
              call cvql(kd,m,q,a0)
           endif
        else if (m.eq.5) then
           if (q.le.1.0.and.kd.eq.2) then
              a0=((6.8e-6*q+1.42e-5)*q2+.0208333)*q2+25.0
           else if (q.le.1.0.and.kd.eq.3) then
              a0=((-6.8e-6*q+1.42e-5)*q2+.0208333)*q2+25.0
           else if (q.le.35.0.and.kd.eq.2) then
              a0=(((2.238231d-5*q-2.983416d-3)*q
     &           +.10706975)*q-.600205)*q+25.93515
           else if (q.le.25.0.and.kd.eq.3) then
              a0=((-7.425364d-4*q+2.18225d-2)*q
     &           +4.16399d-2)*q+24.897
           else
              call cvql(kd,m,q,a0)
           endif
        else if (m.eq.6) then
           if (q.le.1.0) then
              a0=(.4d-6*q2+.0142857)*q2+36.0
           else if (q.le.40.0.and.kd.eq.1) then
              a0=(((-1.66846d-5*q+4.80263d-4)*q
     &           +2.53998d-2)*q-.181233)*q+36.423
           else if (q.le.35.0.and.kd.eq.4) then
              a0=((-4.57146d-4*q+2.16609d-2)*q-2.349616d-2)*q
     &           +35.99251
           else
              call cvql(kd,m,q,a0)
           endif
        else if (m.eq.7) then
           if (q.le.10.0) then
              call cvqm(m,q,a0)
           else if (q.le.50.0.and.kd.eq.2) then
              a0=(((-1.411114d-5*q+9.730514d-4)*q
     &           -3.097887d-3)*q+3.533597d-2)*q+49.0547
           else if (q.le.40.0.and.kd.eq.3) then
              a0=((-3.043872d-4*q+2.05511d-2)*q
     &           -9.16292d-2)*q+49.19035
           else
              call cvql(kd,m,q,a0)
           endif
        else if (m.ge.8) then
           if (q.le.3.*m) then
              call cvqm(m,q,a0)
           else if (q.gt.m*m) then
              call cvql(kd,m,q,a0)
           else
              if (m.eq.8.and.kd.eq.1) then
                 a0=(((8.634308d-6*q-2.100289d-3)*q+.169072)*q
     &              -4.64336)*q+109.4211
              else if (m.eq.8.and.kd.eq.4) then
                 a0=((-6.7842d-5*q+2.2057d-3)*q+.48296)*q+56.59
              else if (m.eq.9.and.kd.eq.2) then
                 a0=(((2.906435d-6*q-1.019893d-3)*q+.1101965)*q
     &              -3.821851)*q+127.6098
              else if (m.eq.9.and.kd.eq.3) then
                 a0=((-9.577289d-5*q+.01043839)*q+.06588934)*q
     &              +78.0198
              else if (m.eq.10.and.kd.eq.1) then
                 a0=(((5.44927d-7*q-3.926119d-4)*q+.0612099)*q
     &              -2.600805)*q+138.1923
              else if (m.eq.10.and.kd.eq.4) then
                 a0=((-7.660143d-5*q+.01132506)*q-.09746023)*q
     &              +99.29494
              else if (m.eq.11.and.kd.eq.2) then
                 a0=(((-5.67615d-7*q+7.152722d-6)*q+.01920291)*q
     &              -1.081583)*q+140.88
              else if (m.eq.11.and.kd.eq.3) then
                 a0=((-6.310551d-5*q+.0119247)*q-.2681195)*q
     &              +123.667
              else if (m.eq.12.and.kd.eq.1) then
                 a0=(((-2.38351d-7*q-2.90139d-5)*q+.02023088)*q
     &              -1.289)*q+171.2723
              else if (m.eq.12.and.kd.eq.4) then
                 a0=(((3.08902d-7*q-1.577869d-4)*q+.0247911)*q
     &              -1.05454)*q+161.471
              endif
           endif
        endif
        return
        end


        subroutine cvql(kd,m,q,a0)
c
c       ========================================================
c       purpose: compute the characteristic value of mathieu
c                functions  for q ò 3m
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c       output:  a0 --- initial characteristic value
c       ========================================================
c
        implicit double precision (a-h,o-z)
        if (kd.eq.1.or.kd.eq.2) w=2.0d0*m+1.0d0
        if (kd.eq.3.or.kd.eq.4) w=2.0d0*m-1.0d0
        w2=w*w
        w3=w*w2
        w4=w2*w2
        w6=w2*w4
        d1=5.0+34.0/w2+9.0/w4
        d2=(33.0+410.0/w2+405.0/w4)/w
        d3=(63.0+1260.0/w2+2943.0/w4+486.0/w6)/w2
        d4=(527.0+15617.0/w2+69001.0/w4+41607.0/w6)/w3
        c1=128.0
        p2=q/w4
        p1=dsqrt(p2)
        cv1=-2.0*q+2.0*w*dsqrt(q)-(w2+1.0)/8.0
        cv2=(w+3.0/w)+d1/(32.0*p1)+d2/(8.0*c1*p2)
        cv2=cv2+d3/(64.0*c1*p1*p2)+d4/(16.0*c1*c1*p2*p2)
        a0=cv1-cv2/(c1*p1)
        return
        end


        subroutine cvqm(m,q,a0)
c
c       =====================================================
c       purpose: compute the characteristic value of mathieu
c                functions for q ó m*m
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c       output:  a0 --- initial characteristic value
c       =====================================================
c
        implicit double precision (a-h,o-z)
        hm1=.5*q/(m*m-1.0)
        hm3=.25*hm1**3/(m*m-4.0)
        hm5=hm1*hm3*q/((m*m-1.0)*(m*m-9.0))
        a0=m*m+q*(hm1+(5.0*m*m+7.0)*hm3
     &     +(9.0*m**4+58.0*m*m+29.0)*hm5)
        return
        end


