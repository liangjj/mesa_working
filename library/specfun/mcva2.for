        program mcva2
c
c       =============================================================
c       purpose: this program calculates a specific characteristic 
c                value of mathieu functions using subroutine cva2
c       input :  m  --- order of mathieu functions
c                q  --- parameter of mathieu functions
c                kd --- case code
c                       kd=1 for cem(x,q)  ( m = 0,2,4,...)
c                       kd=2 for cem(x,q)  ( m = 1,3,5,...)
c                       kd=3 for sem(x,q)  ( m = 1,3,5,...)
c                       kd=4 for sem(x,q)  ( m = 2,4,6,...)
c       output:  a  --- characteristic value
c       example: q = 25.0, m = 0,1,2,...,12
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
c       =============================================================
c
        implicit double precision (a-h,o-z)
        write(*,*)'please enter kd, m and q '
        read(*,*)kd,m,q
        call cva2(kd,m,q,a)
        if (kd.le.2) then
           write(*,20)
        else
           write(*,30)
        endif
        write(*,*)
        write(*,*)'  m         a (or b)'
        write(*,*)'------------------------'
        write(*,10)m,a
10      format(1x,i3,f18.8)
20      format(1x,'characteristic value of even mathieu function, a')
30      format(1x,'characteristic value of odd mathieu function, b')
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
c                 value using an iteration method
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
